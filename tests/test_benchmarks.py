"""Each plot is checked against a known answer, or against an independent
second construction of the same problem."""

import numpy as np
import pytest

import chemeng_plots as cp
from chemeng_plots.datasets import COUNTERCURRENT, DRYING, SYSTEM_1, SYSTEM_2


# ---------- pure helpers ----------

def test_ternary_vertices():
    assert np.allclose(cp.ternary_to_xy(1, 0, 0), (0.5, np.sqrt(3) / 2))
    assert np.allclose(cp.ternary_to_xy(0, 1, 0), (0, 0))
    assert np.allclose(cp.ternary_to_xy(0, 0, 1), (1, 0))
    # wt% and fractions give the same point
    assert np.allclose(cp.ternary_to_xy(20, 30, 50), cp.ternary_to_xy(0.2, 0.3, 0.5))


def test_lever_rule_recovers_mass_ratio():
    W = np.array(cp.ternary_to_xy(21.2, 7.5, 71.5))
    X = np.array(cp.ternary_to_xy(36.9, 58.6, 4.5))
    M = (50 * W + 100 * X) / 150          # 50 kg W + 100 kg X
    assert cp.lever_rule(M, W, X) == pytest.approx(0.5, abs=1e-9)


def test_falling_rate_integral_constant_rate():
    # R constant = 2 over X in [0.1, 0.3]  ->  ∫dX/R = 0.1
    assert cp.falling_rate_integral([0.3, 0.1, 0.2], [2, 2, 2]) == pytest.approx(0.1)


# ---------- plots ----------

@pytest.mark.parametrize("system", [SYSTEM_1, SYSTEM_2])
def test_ternary_data_sums_and_plait_inside(system):
    for p in system["raffinate"] + system["extract"]:
        assert sum(p) == pytest.approx(100, abs=1.0)
    _, _, res = cp.plot_ternary(system["raffinate"], system["extract"], labels=system["labels"])
    px, py = res["plait_xy"]
    assert 0 < py < np.sqrt(3) / 2 and 0 < px < 1
    assert len(res["coords"]["raffinate"]["rows"]) == len(system["raffinate"])
    assert "raffinate" in cp.coords_table(res)


def test_distribution_sides_of_diagonal():
    x1 = [p[0] for p in SYSTEM_1["raffinate"]]
    y1 = [p[0] for p in SYSTEM_1["extract"]]
    x2 = [p[0] for p in SYSTEM_2["raffinate"]]
    y2 = [p[0] for p in SYSTEM_2["extract"]]
    assert all(b > a for a, b in zip(x1, y1))   # K > 1
    assert all(b < a for a, b in zip(x2, y2))   # K < 1
    cp.plot_distribution(x1, y1)


def test_mccabe_thiele_three_stages_reach_known_raffinate():
    d = COUNTERCURRENT
    _, _, res = cp.plot_mccabe_thiele(d["eq_x"], d["eq_y"], d["raffinate_end"],
                                      d["feed_end"], n_stages=3)
    assert res["n_stages"] == 3
    # tutorial answer X_3 = 1.058e-3 (graphical accuracy, 25 %)
    assert res["x_final"] == pytest.approx(1.058e-3, rel=0.25)


def test_mccabe_thiele_design_mode_straight_lines():
    # equilibrium Y = 2X, operating line from (0.02, 0) to (0.12, 0.18): Kremser gives 4 stages
    _, _, res = cp.plot_mccabe_thiele([0, 0.12], [0, 0.24], (0.02, 0.0), (0.12, 0.18))
    assert res["reached_target"]
    assert res["n_stages"] == 4


def test_mccabe_thiele_pinch_warns():
    # operating line above equilibrium near the feed end -> never reaches target
    with pytest.warns(UserWarning, match="not reached"):
        cp.plot_mccabe_thiele([0, 1], [0, 1], (0.0, 0.0), (1.0, 1.0))


def test_drying_falling_rate_integral():
    d = DRYING
    _, _, res = cp.plot_drying_curve(d["X"], d["R"], X_c=d["X_c"], X_2=0.05, R_c=d["R_c"])
    # hand-computed trapezoid answer 0.1276 (15 %)
    assert res["falling_integral"] == pytest.approx(0.1276, rel=0.15)


def test_filtration_fit_recovers_line():
    rng = np.random.default_rng(0)
    V = np.linspace(0.1, 1.0, 8)
    tV = 2.0 * V + 5.0 + rng.normal(0, 0.02, len(V))
    _, _, res = cp.plot_filtration_linear(V, tV)
    assert res["slope"] == pytest.approx(2.0, abs=0.1)
    assert res["intercept"] == pytest.approx(5.0, abs=0.1)
    assert res["K_p"] == pytest.approx(2 * res["slope"])


def test_sle_matches_concentration_space_stepping():
    r = 0.4
    under_n = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    under_r = np.full_like(under_n, r)
    x0 = tuple(float(v) for v in cp.entrainment_to_xy(0.5, r))
    xn = tuple(float(v) for v in cp.entrainment_to_xy(0.05, r))
    _, _, res = cp.plot_sle_right_triangle(under_n, under_r, x0, (0.30, 0.70), xn,
                                           show_construction=True)
    # same problem in concentration space: 45° equilibrium line
    _, _, res2 = cp.plot_mccabe_thiele([0, 1], [0, 1], (0.05, 0.0), (0.5, 0.30))
    assert res["reached_target"]
    assert res["n_stages"] == res2["n_stages"]
    # underflow stage points lie on the constant-entrainment line x_A + x_S = r/(1+r)
    for a, b in res["coords"]["underflow_stages"]["rows"]:
        assert a + b == pytest.approx(r / (1 + r), abs=1e-6)


MS = 0.005  # mutual solubility of the near-immiscible synthetic system


def _near_immiscible(xs):
    raff = [(x, 1 - x - MS, MS) for x in xs]
    extr = [(2 * x, MS, 1 - 2 * x - MS) for x in xs]
    return raff, extr


def _lle_d(raff, extr):
    return cp.plot_lle_d_system(raff, extr, F=(0.12, 0.88, 0), S=(0, 0, 1),
                                RN=(0.02, 1 - 0.02 - MS, MS),
                                E1=(0.18, MS, 1 - 0.18 - MS))


def test_lle_difference_point_matches_immiscible_limit():
    # Reference: stage stepping in solute-free ratio coordinates, where the
    # operating line of a near-immiscible system is straight.
    x = np.linspace(0, 0.12, 200)
    X, Y = x / (1 - x - MS), 2 * x / (1 - 2 * x - MS)
    X_N, X_F, Y_1 = 0.02 / (1 - 0.02 - MS), 0.12 / 0.88, 0.18 / (1 - 0.18 - MS)
    _, _, ref = cp.plot_mccabe_thiele(X, Y, (X_N, 0), (X_F, Y_1))
    ref_xa = np.array(ref["stage_x"]) / (1 + np.array(ref["stage_x"]))

    raff, extr = _near_immiscible(np.r_[0.0, np.linspace(0.0005, 0.12, 40)])
    _, _, res = _lle_d(raff, extr)
    xa = [yy / (np.sqrt(3) / 2) for _, yy in res["coords"]["raffinate_stages"]["rows"]]

    assert res["reached_target"]
    assert res["n_stages"] == ref["n_stages"] == 5
    # stage-by-stage raffinate solute agrees with the reference to within 0.002
    assert np.allclose(xa[:4], ref_xa[:4], atol=2e-3)


def test_lle_warns_when_tie_lines_run_out():
    # tie lines stop at 2 % solute, but a 5th stage needs data below that
    raff, extr = _near_immiscible([0.02, 0.04, 0.06, 0.08, 0.10, 0.12])
    with pytest.warns(UserWarning, match="lower bound"):
        _, _, res = _lle_d(raff, extr)
    assert not res["reached_target"]
    assert res["n_stages"] == 4


def test_lle_solvent_ratio_mode_on_real_data():
    # E1 found from the solvent/feed ratio; target lies inside the tie-line data
    _, _, res = cp.plot_lle_d_system(SYSTEM_1["raffinate"], SYSTEM_1["extract"],
                                     F=(40, 60, 0), S=(0, 0, 100),
                                     RN=(17.04, 82.23, 0.73), S_over_F=0.5)
    assert res["reached_target"]
    assert res["n_stages"] == 2
    px, _ = res["difference_point"]
    assert px < 0          # P lies outside the triangle, on the feed side


def test_lle_requires_e1_or_ratio():
    with pytest.raises(ValueError):
        cp.plot_lle_d_system(SYSTEM_1["raffinate"], SYSTEM_1["extract"],
                             F=(40, 60, 0), S=(0, 0, 100), RN=(5.96, 93.52, 0.52))


def test_psychrometric_state_point():
    pytest.importorskip("psychrolib")
    _, _, res = cp.plot_psychrometric_chart(points=[{"Tdb": 30, "RH": 0.5}])
    Tdb, Twb, RH, Hr, h, Tdew = res["coords"]["states"]["rows"][0]
    # standard chart: 30 °C / 50 % RH -> H ≈ 13.3 g/kg, h ≈ 64 kJ/kg
    assert Hr == pytest.approx(0.0133, abs=0.001)
    assert 60 < h < 70


def test_language_switch():
    cp.set_language("zh")
    try:
        _, ax, _ = cp.plot_filtration_linear([0.1, 0.2, 0.3], [5.2, 5.4, 5.6])
        assert ax.get_title() == "恒压过滤 t/V vs V"
    finally:
        cp.set_language("en")
    _, ax, _ = cp.plot_filtration_linear([0.1, 0.2, 0.3], [5.2, 5.4, 5.6])
    assert ax.get_title().startswith("Constant-pressure")
    with pytest.raises(ValueError):
        cp.set_language("fr")
