#!/usr/bin/env python3
"""Bloomsize: size, analyse, simulate and gate Bloom filters from plain text.

A spec file has one filter per line, ``NAME key=value ...``, plus optional
``gate`` lines for CI. Each filter needs ``n`` (items) and exactly one of
``p`` (target false-positive rate), ``m`` (bit count) or ``bpi`` (bits per item).
Optional keys: ``k`` (hash count override) and ``counter`` (bits per counting
cell, used for overflow-risk analysis). Lines starting with ``#`` are comments.

Exit codes: 0 ok, 1 lint errors or failed gates, 2 usage or file errors.
"""

import argparse
import hashlib
import json
import math
import re
import sys
from dataclasses import asdict, dataclass
from typing import Optional

LN2 = math.log(2)
MODES = ("p", "m", "bpi")
KEYS = {"n", "p", "m", "k", "bpi", "counter"}
GATE_METRICS = ("fpr", "bytes", "overflow")
GATE_RE = re.compile(r"^gate\s+(\S+)\s+(\w+)\s*(<=|>=)\s*(\S+)$")
BPI_SWEEP = (2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24)


class SpecError(Exception):
    def __init__(self, line, message):
        super().__init__(f"line {line}: {message}")
        self.line = line
        self.message = message


@dataclass
class Filter:
    name: str
    line: int
    n: int
    p: Optional[float] = None
    m: Optional[int] = None
    k: Optional[int] = None
    bpi: Optional[float] = None
    counter: Optional[int] = None


@dataclass
class Gate:
    line: int
    name: str
    metric: str
    op: str
    value: float


@dataclass
class Spec:
    filters: list
    gates: list
    errors: list  # (line, message) pairs, sorted by line


@dataclass
class Result:
    name: str
    mode: str
    n: int
    m: int
    bytes: int
    bits_per_item: float
    k: int
    k_float: float
    k_optimal: int
    fpr: float
    fpr_optimal_k: float
    target_p: Optional[float]
    counter_bits: Optional[int]
    counter_overflow_p: Optional[float]
    expected_overflow_counters: Optional[float]
    any_overflow_p: Optional[float]


# ---------------------------------------------------------------- parsing

def _int(text, line, key, minimum=None, maximum=None):
    try:
        value = int(text)
    except ValueError:
        raise SpecError(line, f"{key} must be a whole number, got {text!r}")
    if minimum is not None and value < minimum:
        raise SpecError(line, f"{key} must be at least {minimum}, got {value}")
    if maximum is not None and value > maximum:
        raise SpecError(line, f"{key} must be at most {maximum}, got {value}")
    return value


def _float(text, line, key):
    try:
        value = float(text)
    except ValueError:
        raise SpecError(line, f"{key} must be a number, got {text!r}")
    if not math.isfinite(value):
        raise SpecError(line, f"{key} must be finite, got {text!r}")
    return value


def _parse_filter(line, lineno):
    tokens = line.split()
    name, pairs = tokens[0], tokens[1:]
    if "=" in name:
        raise SpecError(lineno, "filter line must start with a name")
    values = {}
    for token in pairs:
        if "=" not in token:
            raise SpecError(lineno, f"expected key=value, got {token!r}")
        key, _, text = token.partition("=")
        if key not in KEYS:
            raise SpecError(lineno, f"unknown key {key!r} (expected one of {', '.join(sorted(KEYS))})")
        if key in values:
            raise SpecError(lineno, f"duplicate key {key!r}")
        values[key] = text
    if "n" not in values:
        raise SpecError(lineno, "missing n= (number of items)")
    filt = Filter(name=name, line=lineno, n=_int(values["n"], lineno, "n", minimum=1))
    modes = [key for key in MODES if key in values]
    if len(modes) != 1:
        raise SpecError(lineno, "give exactly one of p=, m= or bpi=")
    mode = modes[0]
    if mode == "p":
        filt.p = _float(values["p"], lineno, "p")
        if not 0 < filt.p < 1:
            raise SpecError(lineno, f"p must be between 0 and 1 exclusive, got {values['p']}")
    elif mode == "m":
        filt.m = _int(values["m"], lineno, "m", minimum=1)
    else:
        filt.bpi = _float(values["bpi"], lineno, "bpi")
        if filt.bpi <= 0:
            raise SpecError(lineno, f"bpi must be positive, got {values['bpi']}")
    if "k" in values:
        filt.k = _int(values["k"], lineno, "k", minimum=1)
    if "counter" in values:
        filt.counter = _int(values["counter"], lineno, "counter", minimum=1, maximum=32)
    return filt


def _parse_gate(line, lineno):
    match = GATE_RE.match(line)
    if not match:
        raise SpecError(lineno, "gate must look like 'gate NAME fpr <= 0.01' (metrics: fpr, bytes, overflow)")
    name, metric, op, value = match.groups()
    if metric not in GATE_METRICS:
        raise SpecError(lineno, f"unknown gate metric {metric!r} (expected fpr, bytes or overflow)")
    return Gate(lineno, name, metric, op, _float(value, lineno, "gate value"))


def parse_spec(text):
    filters, gates, errors = [], [], []
    seen = set()
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        try:
            if line.split()[0] == "gate":
                gates.append(_parse_gate(line, lineno))
            else:
                filt = _parse_filter(line, lineno)
                if filt.name in seen:
                    raise SpecError(lineno, f"duplicate filter name {filt.name!r}")
                seen.add(filt.name)
                filters.append(filt)
        except SpecError as exc:
            errors.append((exc.line, exc.message))
    for gate in gates:
        if gate.name not in seen:
            errors.append((gate.line, f"gate refers to unknown filter {gate.name!r}"))
    return Spec(filters, gates, sorted(errors))


# ---------------------------------------------------------------- maths

def optimal_bits(n, p):
    return math.ceil(-n * math.log(p) / LN2 ** 2)


def false_positive_rate(m, n, k):
    return (1.0 - math.exp(-k * n / m)) ** k


def poisson_tail(lam, c):
    """P(X >= c) for X ~ Poisson(lam), summed in log space to avoid overflow."""
    if c < lam:
        cdf = sum(math.exp(i * math.log(lam) - lam - math.lgamma(i + 1)) for i in range(c))
        return max(0.0, 1.0 - cdf)
    total, i = 0.0, c
    while True:
        term = math.exp(i * math.log(lam) - lam - math.lgamma(i + 1))
        total += term
        if term <= total * 1e-16 or i - c > 5_000_000:
            return min(total, 1.0)
        i += 1


def analyse(filt):
    if filt.p is not None:
        m, mode = optimal_bits(filt.n, filt.p), "design"
    elif filt.m is not None:
        m, mode = filt.m, "geometry"
    else:
        m, mode = max(1, math.ceil(filt.bpi * filt.n)), "budget"
    k_float = m / filt.n * LN2
    k_opt = max(1, round(k_float))
    k = filt.k if filt.k is not None else k_opt
    overflow_p = expected = any_p = None
    if filt.counter is not None:
        overflow_p = poisson_tail(filt.n * k / m, 2 ** filt.counter)
        expected = m * overflow_p
        any_p = -math.expm1(m * math.log1p(-overflow_p)) if overflow_p < 1 else 1.0
    return Result(
        name=filt.name, mode=mode, n=filt.n, m=m, bytes=math.ceil(m / 8),
        bits_per_item=m / filt.n, k=k, k_float=k_float, k_optimal=k_opt,
        fpr=false_positive_rate(m, filt.n, k),
        fpr_optimal_k=false_positive_rate(m, filt.n, k_opt),
        target_p=filt.p, counter_bits=filt.counter, counter_overflow_p=overflow_p,
        expected_overflow_counters=expected, any_overflow_p=any_p,
    )


def lint_warnings(filt, res):
    warnings = []
    if res.bits_per_item < 2:
        warnings.append(f"only {res.bits_per_item:.2f} bits per item; FPR is {res.fpr:.3g}")
    if res.k > 2 * res.k_optimal or 2 * res.k < res.k_optimal:
        warnings.append(f"k={res.k} is far from the optimal k={res.k_optimal} for this size")
    if filt.p is not None and res.fpr > filt.p:
        excess = (res.fpr / filt.p - 1) * 100
        warnings.append(
            f"integer k={res.k} leaves FPR {res.fpr:.6g}, {excess:.2f}% above the target p={filt.p:g}"
        )
    if res.any_overflow_p is not None and res.any_overflow_p > 1e-6:
        warnings.append(
            f"{res.counter_bits}-bit counters may saturate: P(any) = {res.any_overflow_p:.3g}"
        )
    return warnings


# ---------------------------------------------------------------- filter simulation

def _positions(key, m, k):
    digest = hashlib.blake2b(key, digest_size=16).digest()
    h1 = int.from_bytes(digest[:8], "little")
    h2 = int.from_bytes(digest[8:], "little") | 1
    return [(h1 + i * h2) % m for i in range(k)]


def simulate_fpr(m, n, k, trials, seed):
    """Insert n members and probe with trials non-members; return false positives."""
    bits = bytearray((m + 7) // 8)
    for i in range(n):
        for pos in _positions(f"{seed}:member:{i}".encode(), m, k):
            bits[pos >> 3] |= 1 << (pos & 7)
    hits = 0
    for i in range(trials):
        for pos in _positions(f"{seed}:probe:{i}".encode(), m, k):
            if not bits[pos >> 3] & (1 << (pos & 7)):
                break
        else:
            hits += 1
    return hits


def wilson_interval(hits, trials, z=1.96):
    if trials == 0:
        return 0.0, 0.0
    p = hits / trials
    denom = 1 + z * z / trials
    centre = (p + z * z / (2 * trials)) / denom
    half = z * math.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


# ---------------------------------------------------------------- gates

def evaluate_gates(results, gates):
    rows = []
    for gate in gates:
        res = results[gate.name]
        if gate.metric == "fpr":
            actual = res.fpr
        elif gate.metric == "bytes":
            actual = res.bytes
        else:
            if res.counter_overflow_p is None:
                rows.append({"gate": gate, "actual": None, "passed": False,
                             "note": "overflow gates need counter= on the filter"})
                continue
            actual = res.counter_overflow_p
        passed = actual <= gate.value if gate.op == "<=" else actual >= gate.value
        rows.append({"gate": gate, "actual": actual, "passed": passed, "note": ""})
    return rows


# ---------------------------------------------------------------- rendering

def render_table(results):
    header = ["filter", "mode", "items", "bits", "bytes", "bits/item", "k", "k_opt", "FPR", "counter overflow"]
    rows = [header]
    for res in results:
        overflow = "-"
        if res.counter_overflow_p is not None:
            overflow = f"{res.counter_overflow_p:.3g}/cell, {res.expected_overflow_counters:.3g} cells"
        rows.append([
            res.name, res.mode, f"{res.n:,}", f"{res.m:,}", f"{res.bytes:,}",
            f"{res.bits_per_item:.2f}", str(res.k), str(res.k_optimal),
            f"{res.fpr:.3g}", overflow,
        ])
    widths = [max(len(row[i]) for row in rows) for i in range(len(header))]
    lines = ["  ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)).rstrip() for row in rows]
    return "\n".join(lines)


def _load(path):
    with open(path, encoding="utf-8") as handle:
        return parse_spec(handle.read())


def _print_errors(path, errors):
    for line, message in errors:
        print(f"{path}:{line}: error: {message}", file=sys.stderr)


def _analyse_all(spec):
    return [analyse(filt) for filt in spec.filters]


# ---------------------------------------------------------------- commands

def cmd_report(args):
    spec = _load(args.file)
    if spec.errors:
        _print_errors(args.file, spec.errors)
        return 1
    results = _analyse_all(spec)
    if args.json:
        print(json.dumps([asdict(res) for res in results], indent=2))
    else:
        print(render_table(results))
    return 0


def cmd_lint(args):
    spec = _load(args.file)
    _print_errors(args.file, spec.errors)
    warning_count = 0
    for filt in spec.filters:
        res = analyse(filt)
        for message in lint_warnings(filt, res):
            warning_count += 1
            print(f"{args.file}:{filt.line}: warning: {filt.name}: {message}")
    if spec.errors:
        print(f"{len(spec.errors)} error(s), {warning_count} warning(s)")
        return 1
    print(f"OK: {len(spec.filters)} filter(s), {warning_count} warning(s)")
    return 0


def cmd_check(args):
    spec = _load(args.file)
    if spec.errors:
        _print_errors(args.file, spec.errors)
        return 1
    if not spec.gates:
        print("no gates defined; nothing to check")
        return 0
    results = {res.name: res for res in _analyse_all(spec)}
    failed = 0
    for row in evaluate_gates(results, spec.gates):
        gate = row["gate"]
        label = f"gate {gate.name} {gate.metric} {gate.op} {gate.value:g}"
        if row["actual"] is None:
            print(f"FAIL {label}: {row['note']}")
            failed += 1
            continue
        status = "PASS" if row["passed"] else "FAIL"
        failed += 0 if row["passed"] else 1
        print(f"{status} {label}: actual {row['actual']:.4g}")
    print(f"{len(spec.gates) - failed} of {len(spec.gates)} gate(s) passed")
    return 0 if failed == 0 else 1


def cmd_simulate(args):
    spec = _load(args.file)
    if spec.errors:
        _print_errors(args.file, spec.errors)
        return 1
    filters = [f for f in spec.filters if args.name in (None, f.name)]
    if not filters:
        print(f"no filter named {args.name!r}", file=sys.stderr)
        return 2
    for filt in filters:
        res = analyse(filt)
        n_sim = min(res.n, args.max_items)
        hits = simulate_fpr(res.m, n_sim, res.k, args.trials, args.seed)
        low, high = wilson_interval(hits, args.trials)
        theory = false_positive_rate(res.m, n_sim, res.k)
        verdict = "consistent" if low <= theory <= high else "MISMATCH"
        print(f"{filt.name}: m={res.m:,} k={res.k} items={n_sim:,}"
              + (f" (capped from {res.n:,})" if n_sim < res.n else ""))
        print(f"  probes={args.trials:,} false_positives={hits:,} "
              f"empirical={hits / args.trials:.4g} (95% CI {low:.4g}-{high:.4g})")
        print(f"  theory={theory:.4g} -> {verdict}")
    return 0


def cmd_sweep(args):
    spec = _load(args.file)
    if spec.errors:
        _print_errors(args.file, spec.errors)
        return 1
    matches = [f for f in spec.filters if f.name == args.name]
    if not matches:
        print(f"no filter named {args.name!r}", file=sys.stderr)
        return 2
    n = matches[0].n
    rows = [["bits/item", "bits", "bytes", "k", "FPR"]]
    for bpi in BPI_SWEEP:
        m = math.ceil(bpi * n)
        k = max(1, round(m / n * LN2))
        rows.append([f"{bpi}", f"{m:,}", f"{math.ceil(m / 8):,}", str(k),
                     f"{false_positive_rate(m, n, k):.3g}"])
    widths = [max(len(row[i]) for row in rows) for i in range(len(rows[0]))]
    print(f"{args.name}: n={n:,}")
    for row in rows:
        print("  ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)).rstrip())
    return 0


def cmd_size(args):
    if (args.p is None) + (args.m is None) + (args.bpi is None) != 2:
        print("size: give exactly one of --p, --m or --bpi", file=sys.stderr)
        return 2
    parts = [f"n={args.n}"]
    for key, value in (("p", args.p), ("m", args.m), ("bpi", args.bpi), ("k", args.k),
                       ("counter", args.counter)):
        if value is not None:
            parts.append(f"{key}={value}")
    spec = parse_spec("design " + " ".join(parts))
    if spec.errors:
        for _, message in spec.errors:
            print(f"size: {message}", file=sys.stderr)
        return 2
    res = analyse(spec.filters[0])
    if args.json:
        print(json.dumps(asdict(res), indent=2))
    else:
        print(render_table([res]))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(prog="bloomsize", description="Size, analyse, simulate and gate Bloom filters.")
    sub = parser.add_subparsers(dest="command", required=True)

    report = sub.add_parser("report", help="print a sizing and FPR table for every filter")
    report.add_argument("file")
    report.add_argument("--json", action="store_true", help="emit JSON instead of a table")
    report.set_defaults(func=cmd_report)

    lint = sub.add_parser("lint", help="report line-numbered errors and warnings")
    lint.add_argument("file")
    lint.set_defaults(func=cmd_lint)

    check = sub.add_parser("check", help="evaluate gate lines for CI (exit 1 on failure)")
    check.add_argument("file")
    check.set_defaults(func=cmd_check)

    simulate = sub.add_parser("simulate", help="measure FPR with a real hash-based filter")
    simulate.add_argument("file")
    simulate.add_argument("--name", help="only simulate this filter")
    simulate.add_argument("--trials", type=int, default=20000, help="non-member probes (default 20000)")
    simulate.add_argument("--seed", default="bloomsize", help="seed for reproducible keys")
    simulate.add_argument("--max-items", type=int, default=200000,
                          help="cap inserted items to keep runs fast (default 200000)")
    simulate.set_defaults(func=cmd_simulate)

    sweep = sub.add_parser("sweep", help="show FPR across bits-per-item for one filter")
    sweep.add_argument("file")
    sweep.add_argument("name")
    sweep.set_defaults(func=cmd_sweep)

    size = sub.add_parser("size", help="analyse one filter without a spec file")
    size.add_argument("--n", type=int, required=True, help="number of items")
    size.add_argument("--p", type=float, help="target false-positive rate")
    size.add_argument("--m", type=int, help="bit count")
    size.add_argument("--bpi", type=float, help="bits per item")
    size.add_argument("--k", type=int, help="hash count override")
    size.add_argument("--counter", type=int, help="bits per counting cell")
    size.add_argument("--json", action="store_true", help="emit JSON")
    size.set_defaults(func=cmd_size)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except FileNotFoundError as exc:
        print(f"bloomsize: {exc.filename}: no such file", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
