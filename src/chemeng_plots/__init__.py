"""chemeng-plots: tested graphical-method plots for chemical engineering separations.

Every ``plot_*`` function returns ``(fig, ax, result)`` — the figure plus the
numbers read off it (stage count, integral, slope...) and every plotted
coordinate in ``result["coords"]``.
"""

from ._i18n import get_language, set_language
from .plots import (
    coords_table,
    entrainment_to_xy,
    falling_rate_integral,
    lever_rule,
    plot_distribution,
    plot_drying_curve,
    plot_filtration_linear,
    plot_lle_d_system,
    plot_mccabe_thiele,
    plot_psychrometric_chart,
    plot_sle_right_triangle,
    plot_ternary,
    ternary_to_xy,
)
from .style import PALETTE, apply_chemeng_style, apply_pro_style, check_axis_range

__version__ = "0.1.0"

__all__ = [
    "plot_ternary",
    "plot_distribution",
    "plot_mccabe_thiele",
    "plot_drying_curve",
    "plot_filtration_linear",
    "plot_sle_right_triangle",
    "plot_psychrometric_chart",
    "plot_lle_d_system",
    "ternary_to_xy",
    "lever_rule",
    "falling_rate_integral",
    "entrainment_to_xy",
    "coords_table",
    "apply_chemeng_style",
    "apply_pro_style",
    "check_axis_range",
    "PALETTE",
    "set_language",
    "get_language",
]
