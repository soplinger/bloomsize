# Bloomsize

**Size, analyse, simulate and gate Bloom filters from one plain-text file, offline.**

Bloomsize answers the questions you hit before shipping a Bloom filter: how many bits and hash functions do I need for my false-positive target, what FPR will a given memory budget actually give, will my counting filter's counters overflow, and does the real filter behave like the formula says. Specs are plain text, errors carry line numbers, and `check` exits non-zero in CI when a sizing gate is violated.

It uses only the Python 3 standard library. There are no servers, accounts, API keys or network calls.

## Features

1. **Design mode** – give `n=` and a target `p=` to get the optimal bit count and hash count.
2. **Geometry mode** – give `m=` (bits) and `n=` to see the FPR you will actually get, with an optional `k=` override.
3. **Budget mode** – give `bpi=` (bits per item) to size from a memory budget.
4. **Exact FPR** – `(1 - e^(-kn/m))^k` for the chosen hash count, and for the optimal integer `k` for comparison.
5. **Counting-filter overflow risk** – with `counter=` bits per cell, the per-cell overflow probability is computed from a Poisson tail, along with the expected number of saturated cells and the chance that any cell saturates.
6. **Line-numbered lint** – every parse error is reported with its file and line (duplicate names, unknown keys, bad numbers, conflicting modes, gates pointing at unknown filters), and warnings flag under-2-bit-per-item filters, hash counts far from optimal, rounding that pushes FPR above target, and counter overflow risk.
7. **CI gates** – `gate NAME fpr <= 0.001`, `gate NAME bytes <= 600000` or `gate NAME overflow <= 1e-6`. `check` prints PASS/FAIL and exits 1 if any gate fails.
8. **Simulation** – `simulate` builds a real double-hashed bit array (BLAKE2b, seeded and reproducible), probes it with non-members, and compares the empirical FPR against the formula using a 95% Wilson interval.
9. **Bits-per-item sweep** – `sweep` shows bits, bytes, k and FPR across a range of budgets for one filter.
10. **One-off sizing** – `size --n 1000000 --p 0.01` analyses a filter without writing a spec file.
11. **Table or JSON output** – a readable aligned table by default, and `report --json` for scripts.
12. **Deterministic results** – no randomness or clock is used, so the same spec always gives the same numbers.

## Quick start

The repo ships a spec at `examples/services.bloom`. It looks like this:

```text
seen_urls      n=2000000  p=0.001
blocklist      n=150000   bpi=12
session_ids    n=500000   m=4000000 k=5
recent_writes  n=80000    p=0.01 counter=4

gate seen_urls      fpr <= 0.0011
gate blocklist      fpr <= 0.004
gate session_ids    bytes <= 600000
gate recent_writes  overflow <= 1e-6
```

Then run:

```sh
python3 bloomsize.py report examples/services.bloom
```

This prints one row per filter with its mode, item count, bit count, bytes, bits per item, hash count, optimal hash count, FPR, and counter overflow figures. To check the sizing gates in CI:

```sh
python3 bloomsize.py check examples/services.bloom
```

It prints one PASS or FAIL line per gate and exits 0 only when all of them pass.

## Spec format

- One filter per line: `NAME key=value ...`.
- Keys: `n` (required), exactly one of `p`, `m` or `bpi`, and optionally `k` and `counter`.
- `gate NAME METRIC OP VALUE`, where METRIC is `fpr`, `bytes` or `overflow` and OP is `<=` or `>=`.
- `#` starts a comment, and blank lines are ignored.

## Commands

| Command | Purpose | Exit code |
| --- | --- | --- |
| `report FILE [--json]` | Sizing and FPR table, or JSON | 0, or 1 on spec errors |
| `lint FILE` | Errors with line numbers and warnings | 1 on errors |
| `check FILE` | Evaluate `gate` lines | 1 if any gate fails |
| `simulate FILE [--name N] [--trials T] [--seed S] [--max-items M]` | Empirical FPR versus theory | 0 |
| `sweep FILE NAME` | FPR across bits-per-item values | 0 |
| `size --n N (--p P \| --m M \| --bpi B) [--k K] [--counter C] [--json]` | One-off analysis | 0, or 2 on bad flags |

Exit code 2 is used for usage errors and missing files.

## Requirements

Python 3.8 or newer. No third-party packages are needed to run the tool. Tests need `pytest`.

## Running tests

```sh
pip install pytest
python3 -m pytest -q
```

The tests cover the formulas, parser error reporting with line numbers, counter-overflow tails, the simulation's determinism and accuracy, gates, lint warnings, JSON output, the one-off `size` command, `sweep`, and the usage and missing-file exit codes.

**Status note:** this build was produced in an environment where running Python and pytest was not permitted, so the test suite and the example commands above have not been executed yet. Run the commands above to confirm them.

## Notes on the maths

- Optimal bits: `m = ceil(-n ln p / (ln 2)^2)`.
- Optimal hashes: `k = round(m/n · ln 2)`. Because `k` is an integer, the achieved FPR can slightly exceed the target; lint warns when it does.
- Counter overflow: each item increments `k` cells, so a cell's count is approximated as Poisson with mean `nk/m`. A `b`-bit counter saturates when the count reaches `2^b`.
- Simulation uses double hashing (`h1 + i·h2 mod m`), so it can differ slightly from the independent-hash formula, which is why the check uses a confidence interval.
