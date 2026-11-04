"""Matplotlib styles for engineering figures.

- :func:`apply_chemeng_style` — black-and-white engineering style (grey levels +
  line-style changes, never relies on colour). Used by every ``plot_*`` function.
- :func:`apply_pro_style` — colour-blind-safe Okabe-Ito colour style for report
  figures you draw yourself.
- :func:`check_axis_range` — warns when outliers squash the main trend.
"""

from __future__ import annotations

import warnings

import matplotlib as mpl
import numpy as np

#: Okabe-Ito colour-blind-safe palette.
PALETTE = {
    "blue": "#0072B2",
    "orange": "#D55E00",
    "green": "#009E73",
    "purple": "#CC79A7",
    "yellow": "#E69F00",
    "sky": "#56B4E9",
    "gray": "#555555",
    "black": "#000000",
}

# CJK fonts first so Chinese labels render; DejaVu Sans covers everything else.
_FONT_CHAIN = ["Microsoft YaHei", "SimHei", "PingFang SC", "Noto Sans CJK SC", "DejaVu Sans"]


def apply_chemeng_style() -> None:
    """Black-and-white engineering style (sets global ``rcParams``).

    - grey levels + line styles instead of colour (prints well, homework-friendly)
    - inward ticks with minor ticks
    - title 14 / label 12 / tick 10 pt
    - ``savefig``: 200 dpi, tight bounding box
    """
    mpl.rcParams["font.sans-serif"] = list(_FONT_CHAIN)
    mpl.rcParams["axes.unicode_minus"] = False

    mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
        color=["#000000", "#555555", "#888888", "#BBBBBB"]
    ) + mpl.cycler(linestyle=["-", "--", ":", "-."])

    mpl.rcParams["axes.titlesize"] = 14
    mpl.rcParams["axes.labelsize"] = 12
    mpl.rcParams["xtick.labelsize"] = 10
    mpl.rcParams["ytick.labelsize"] = 10
    mpl.rcParams["legend.fontsize"] = 10

    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams["ytick.direction"] = "in"
    mpl.rcParams["xtick.minor.visible"] = True
    mpl.rcParams["ytick.minor.visible"] = True

    mpl.rcParams["axes.linewidth"] = 1.0
    mpl.rcParams["grid.linestyle"] = ":"
    mpl.rcParams["grid.alpha"] = 0.4

    mpl.rcParams["figure.dpi"] = 100
    mpl.rcParams["savefig.dpi"] = 200
    mpl.rcParams["savefig.bbox"] = "tight"


def apply_pro_style() -> None:
    """Colour report style: Okabe-Ito palette, bold titles, light grid, 200 dpi."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"] + _FONT_CHAIN[:-1],
        "mathtext.fontset": "dejavusans",
        "axes.titlesize": 13, "axes.titleweight": "bold",
        "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "legend.fontsize": 9.5, "legend.frameon": True,
        "legend.framealpha": 0.9, "legend.edgecolor": "0.8",
        "axes.linewidth": 1.0, "axes.edgecolor": "0.3",
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.grid": True, "grid.linestyle": "-",
        "grid.linewidth": 0.5, "grid.alpha": 0.35,
        "axes.unicode_minus": False,
        "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
    })
    mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
        color=[PALETTE[k] for k in ("blue", "orange", "green", "purple", "yellow", "sky")]
    )


def check_axis_range(ax, data, threshold: float = 0.05) -> bool:
    """Warn when outliers squash the main trend of ``data`` on ``ax``.

    Compares the inter-quartile span of ``data`` with the current y-axis span.
    If the IQR covers less than ``threshold`` of the axis, emits a
    :class:`UserWarning` and returns ``False``.
    """
    arr = np.asarray(data, dtype=float).ravel()
    arr = arr[~np.isnan(arr)]
    if len(arr) < 4:
        return True

    ymin, ymax = ax.get_ylim()
    axis_span = ymax - ymin
    if axis_span <= 0:
        return True

    q25, q75 = np.percentile(arr, [25, 75])
    ratio = (q75 - q25) / axis_span
    if ratio < threshold:
        warnings.warn(
            f"check_axis_range: the data's IQR covers only {ratio:.1%} of the y axis; "
            "outliers are squashing the main trend. Consider ax.set_ylim(...) or "
            "ax.set_yscale('log').",
            stacklevel=2,
        )
        return False
    return True
