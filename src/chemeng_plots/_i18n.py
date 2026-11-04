"""Figure text in English and Chinese. Switch with :func:`set_language`."""

from __future__ import annotations

_STRINGS = {
    # shared ternary chrome
    "edge_left": {"en": "(left-edge scale, wt%)", "zh": "(左边刻度 ↗ wt%)"},
    "edge_bottom": {"en": "(bottom-edge scale, wt%)", "zh": "(底边刻度 wt%)"},
    "edge_right": {"en": "(right-edge scale, wt%)", "zh": "(右边刻度 wt%)"},
    "binodal": {"en": "Binodal curve", "zh": "Binodal 曲线"},
    "point": {"en": "Point", "zh": "点"},
    # ternary
    "ternary_title": {"en": "Ternary phase diagram", "zh": "三元相图"},
    "raffinate": {"en": "Raffinate (R)", "zh": "萃余相 (R)"},
    "extract": {"en": "Extract (E)", "zh": "萃取相 (E)"},
    "fs_line": {"en": "F–S line", "zh": "F-S 线"},
    # distribution
    "dist_title": {"en": "Distribution curve", "zh": "分配曲线"},
    "dist_x": {"en": "Solute in raffinate x (wt%)", "zh": "溶质 in 萃余相 x (wt%)"},
    "dist_y": {"en": "Solute in extract y (wt%)", "zh": "溶质 in 萃取相 y (wt%)"},
    "dist_label": {"en": "System", "zh": "体系"},
    "op_point": {"en": "Operating point", "zh": "工作点"},
    # McCabe-Thiele
    "mt_title": {"en": "Countercurrent extraction — stage stepping",
                 "zh": "逆流多级萃取 — 阶梯法"},
    "mt_x": {"en": "X = solute / carrier (kg/kg)", "zh": "X = 溶质/载体 (kg/kg)"},
    "mt_y": {"en": "Y = solute / solvent (kg/kg)", "zh": "Y = 溶质/溶剂 (kg/kg)"},
    "eq_curve": {"en": "Equilibrium curve", "zh": "平衡曲线"},
    "op_line": {"en": "Operating line", "zh": "操作线"},
    # drying
    "dry_title": {"en": "Drying rate curve R–X", "zh": "干燥速率曲线 R-X"},
    "dry_x": {"en": "Free moisture X (kg water / kg dry solid)",
              "zh": "自由含水率 X (kg/kg 干基)"},
    "dry_y": {"en": "Drying rate R (kg/(m²·h))", "zh": "干燥速率 R (kg/(m²·h))"},
    "dry_data": {"en": "R–X data", "zh": "R-X 数据"},
    "dry_const": {"en": "Constant-rate $R_c$ ≈ {:.2f}", "zh": "恒速段 $R_c$ ≈ {:.2f}"},
    "dry_fall": {"en": "Falling-rate integral", "zh": "降速段积分区"},
    # filtration
    "filt_title": {"en": "Constant-pressure filtration t/V vs V", "zh": "恒压过滤 t/V vs V"},
    "filt_x": {"en": "Filtrate volume V (m³)", "zh": "滤液体积 V (m³)"},
    "filt_y": {"en": "t/V (s/m³)", "zh": "t/V (s/m³)"},
    "filt_data": {"en": "Data", "zh": "实验数据"},
    "filt_fit": {"en": "Fit: slope={:.3g}, intercept={:.3g}",
                 "zh": "拟合：斜率={:.3g}, 截距={:.3g}"},
    # SLE
    "sle_title": {"en": "Leaching (SLE) right-triangle diagram — countercurrent",
                  "zh": "SLE 直角三角图 — 逆流多级"},
    "underflow": {"en": "Underflow line (entrainment)", "zh": "Underflow 线（夹带）"},
    "overflow": {"en": "Overflow line (hypotenuse)", "zh": "Overflow 线（斜边）"},
    "sle_x": {"en": "$x_A$ = {} mass fraction", "zh": "$x_A$ = {} 质量分数"},
    "sle_y": {"en": "$x_S$ = {} mass fraction", "zh": "$x_S$ = {} 质量分数"},
    # psychrometric
    "psy_title": {"en": "Psychrometric chart", "zh": "湿度图 Psychrometric Chart"},
    "psy_sat": {"en": "Saturation (RH = 100%)", "zh": "饱和线 (RH=100%)"},
    "psy_x": {"en": "Dry-bulb temperature (°C)", "zh": "干球温度 Dry-bulb T (°C)"},
    "psy_y": {"en": "Humidity ratio H (kg/kg dry air)",
              "zh": "湿度比 Humidity ratio H (kg/kg 干空气)"},
    # LLE difference point
    "lled_title": {"en": "Partially miscible countercurrent LLE — difference point",
                   "zh": "LLE [D] 部分互溶逆流 — 差点 P 图解"},
    "lled_offchart": {"en": "Difference point P is off-chart\nat (x, y) = ({:.2f}, {:.2f})",
                      "zh": "差点 P 在画外\n(x, y) = ({:.2f}, {:.2f})"},
}

_LANG = "en"


def set_language(lang: str) -> None:
    """Set the language of titles, axis labels and legends: ``"en"`` or ``"zh"``."""
    global _LANG
    if lang not in ("en", "zh"):
        raise ValueError(f"language must be 'en' or 'zh', got {lang!r}")
    _LANG = lang


def get_language() -> str:
    return _LANG


def t(key: str) -> str:
    return _STRINGS[key][_LANG]
