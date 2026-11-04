"""Render the README gallery into docs/images/.

    python examples/gallery.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402

import chemeng_plots as cp  # noqa: E402
from chemeng_plots.datasets import COUNTERCURRENT, DRYING, SYSTEM_1  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "docs" / "images"
OUT.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.savefig(OUT / name, dpi=110)
    print("saved", name)


fig, _, _ = cp.plot_ternary(SYSTEM_1["raffinate"], SYSTEM_1["extract"],
                            labels=SYSTEM_1["labels"],
                            title="Acetone – water – trichloroethane")
save(fig, "ternary.png")

d = COUNTERCURRENT
fig, _, res = cp.plot_mccabe_thiele(d["eq_x"], d["eq_y"], d["raffinate_end"], d["feed_end"],
                                    n_stages=3, title="Countercurrent extraction, 3 stages")
save(fig, "mccabe_thiele.png")

r = 0.4
under_n = np.linspace(0, 0.6, 7)
fig, _, res = cp.plot_sle_right_triangle(
    under_n, np.full_like(under_n, r),
    feed_underflow=tuple(cp.entrainment_to_xy(0.5, r)),
    exit_overflow=(0.30, 0.70),
    final_underflow=tuple(cp.entrainment_to_xy(0.05, r)),
    show_construction=True, title="Countercurrent leaching, constant underflow")
save(fig, "leaching.png")

fig, _, res = cp.plot_lle_d_system(SYSTEM_1["raffinate"], SYSTEM_1["extract"],
                                   F=(40, 60, 0), S=(0, 0, 100), RN=(17.04, 82.23, 0.73),
                                   S_over_F=0.5, labels=SYSTEM_1["labels"],
                                   title="Difference-point stepping, S/F = 0.5")
save(fig, "difference_point.png")

fig, _, _ = cp.plot_psychrometric_chart(points=[{"Tdb": 30, "RH": 0.5}, {"Tdb": 40, "Twb": 25}])
save(fig, "psychrometric.png")

fig, _, _ = cp.plot_drying_curve(DRYING["X"], DRYING["R"], X_c=DRYING["X_c"], X_2=0.05,
                                 R_c=DRYING["R_c"])
save(fig, "drying.png")
