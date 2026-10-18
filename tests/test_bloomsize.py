import json
import math
from pathlib import Path

import pytest

import bloomsize as bs

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "services.bloom"


def write_spec(tmp_path, text):
    path = tmp_path / "spec.bloom"
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_optimal_bits_matches_textbook_formula():
    # n=1000, p=1%: m = -n ln p / (ln 2)^2 ~= 9585.06 -> 9586
    assert bs.optimal_bits(1000, 0.01) == 9586


def test_optimal_k_rounds_m_over_n_ln2():
    m, n = 9586, 1000
    assert bs.analyse(bs.parse_spec("f n=1000 m=9586").filters[0]).k == round(m / n * math.log(2))


def test_fpr_formula_known_value():
    # m/n = 10, k = 7 gives roughly 0.82% (classic table value)
    assert bs.false_positive_rate(10_000, 1_000, 7) == pytest.approx(0.00819, abs=2e-5)


def test_design_mode_meets_target_p():
    res = bs.analyse(bs.parse_spec("f n=50000 p=0.001").filters[0])
    assert res.mode == "design"
    assert res.fpr <= 0.0012  # rounding k may nudge it slightly above p


def test_budget_mode_uses_bits_per_item():
    res = bs.analyse(bs.parse_spec("f n=1000 bpi=8.5").filters[0])
    assert res.m == math.ceil(8.5 * 1000)
    assert res.mode == "budget"


def test_k_override_is_used():
    res = bs.analyse(bs.parse_spec("f n=1000 m=10000 k=3").filters[0])
    assert res.k == 3
    assert res.k_optimal == 7


def test_parse_reports_every_error_with_line_numbers():
    text = "\n".join([
        "a n=10 p=0.01",
        "b n=10 p=1.5",
        "c n=ten m=100",
        "d n=10",
        "e n=10 m=100 bpi=4",
        "a n=10 m=100",
        "f n=10 m=100 zz=1",
        "gate nosuch fpr <= 0.1",
    ])
    spec = bs.parse_spec(text)
    lines = [line for line, _ in spec.errors]
    assert lines == [2, 3, 4, 5, 6, 7, 8]
    assert [f.name for f in spec.filters] == ["a"]


def test_gate_parsing_and_bad_metric():
    spec = bs.parse_spec("f n=10 m=100\ngate f latency <= 3")
    assert spec.errors == [(2, "unknown gate metric 'latency' (expected fpr, bytes or overflow)")]


def test_counter_overflow_tail_is_small_for_wide_counters():
    narrow = bs.analyse(bs.parse_spec("f n=100000 m=800000 k=4 counter=2").filters[0])
    wide = bs.analyse(bs.parse_spec("f n=100000 m=800000 k=4 counter=8").filters[0])
    assert narrow.counter_overflow_p > wide.counter_overflow_p
    assert wide.counter_overflow_p < 1e-12
    assert 0 < narrow.any_overflow_p <= 1


def test_poisson_tail_matches_direct_sum():
    lam, c = 3.0, 6
    direct = 1 - sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(c))
    assert bs.poisson_tail(lam, c) == pytest.approx(direct, rel=1e-9)


def test_simulation_is_deterministic_and_close_to_theory():
    m, n, k = 20_000, 2_000, 7
    hits_a = bs.simulate_fpr(m, n, k, 20_000, "seed-x")
    hits_b = bs.simulate_fpr(m, n, k, 20_000, "seed-x")
    assert hits_a == hits_b
    theory = bs.false_positive_rate(m, n, k)
    assert abs(hits_a / 20_000 - theory) < 0.3 * theory


def test_wilson_interval_bounds():
    low, high = bs.wilson_interval(10, 100)
    assert 0 <= low < 0.1 < high <= 1


def test_lint_warns_on_low_bits_and_far_k(tmp_path):
    text = "tight n=1000 m=1000 k=5\n"
    filt = bs.parse_spec(text).filters[0]
    warnings = bs.lint_warnings(filt, bs.analyse(filt))
    assert any("bits per item" in w for w in warnings)
    assert any("far from the optimal" in w for w in warnings)


def test_example_file_is_clean_and_gates_pass(capsys):
    assert bs.main(["lint", str(EXAMPLE)]) == 0
    assert bs.main(["check", str(EXAMPLE)]) == 0
    out = capsys.readouterr().out
    assert "4 of 4 gate(s) passed" in out


def test_check_fails_when_gate_is_violated(tmp_path, capsys):
    path = write_spec(tmp_path, "f n=1000 bpi=4\ngate f fpr <= 0.0001\n")
    assert bs.main(["check", path]) == 1
    assert "FAIL gate f fpr" in capsys.readouterr().out


def test_lint_returns_1_and_prints_line_errors(tmp_path, capsys):
    path = write_spec(tmp_path, "good n=10 m=100\nbad n=10\n")
    assert bs.main(["lint", path]) == 1
    err = capsys.readouterr().err
    assert f"{path}:2: error: give exactly one of p=, m= or bpi=" in err


def test_report_json_is_parseable(tmp_path, capsys):
    path = write_spec(tmp_path, "f n=1000 p=0.01\n")
    assert bs.main(["report", path, "--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data[0]["name"] == "f"
    assert data[0]["m"] == 9586


def test_size_one_off_and_bad_combination(capsys):
    assert bs.main(["size", "--n", "1000", "--p", "0.01"]) == 0
    assert "design" in capsys.readouterr().out
    assert bs.main(["size", "--n", "1000", "--p", "0.01", "--bpi", "4"]) == 2


def test_simulate_reports_consistency(tmp_path, capsys):
    path = write_spec(tmp_path, "f n=500 bpi=8\n")
    assert bs.main(["simulate", path, "--trials", "5000"]) == 0
    out = capsys.readouterr().out
    assert "probes=5,000" in out
    assert "-> consistent" in out or "-> MISMATCH" in out


def test_sweep_prints_each_bits_per_item_row(tmp_path, capsys):
    path = write_spec(tmp_path, "f n=1000 bpi=8\n")
    assert bs.main(["sweep", path, "f"]) == 0
    out = capsys.readouterr().out
    assert out.count("\n") >= len(bs.BPI_SWEEP) + 1


def test_missing_file_returns_usage_error(tmp_path, capsys):
    assert bs.main(["report", str(tmp_path / "nope.bloom")]) == 2
