"""Graphical-method plots for separation processes / unit operations.

Every ``plot_*`` function returns ``(fig, ax, result)``. ``result`` is a dict
holding the computed numbers (stage count, integral, slope...) plus
``result["coords"]``: every plotted point, so you can copy the construction
onto graph paper or check it. ``print(coords_table(result))`` formats it.

Ternary convention (equilateral triangle): ``(a, b, c)`` = (solute at the top,
carrier at bottom-left, solvent at bottom-right). wt% or fractions both work —
points are normalised internally.
"""

from __future__ import annotations

import warnings

import matplotlib.pyplot as plt
import numpy as np

from ._i18n import t
from .style import apply_chemeng_style, check_axis_range

# numpy >= 2.0 renamed trapz -> trapezoid
_trapz = getattr(np, "trapezoid", None) or np.trapz

_MAX_STAGES = 50

# Equilateral-triangle vertices (Cartesian)
_TOP = np.array([0.5, np.sqrt(3) / 2])   # solute A
_LEFT = np.array([0.0, 0.0])             # carrier B
_RIGHT = np.array([1.0, 0.0])            # solvent C


# ============================================================
# Pure helpers (no plotting, unit-testable)
# ============================================================

def ternary_to_xy(a, b, c):
    """Ternary composition -> Cartesian ``(x, y)``.

    ``a`` = top vertex (solute), ``b`` = bottom-left (carrier),
    ``c`` = bottom-right (solvent). Inputs are normalised, so wt% is fine.
    Accepts scalars or equal-length arrays.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)
    total = a + b + c
    a, c = a / total, c / total
    x = 0.5 * a + c
    y = a * (np.sqrt(3) / 2)
    return x, y


def lever_rule(xy_M, xy_1, xy_2):
    """Lever rule for ``M = stream 1 + stream 2``: returns ``mass_1 / mass_2``.

    Masses are inversely proportional to their arm lengths from M:
    ``mass_1 / mass_2 = dist(M, 2) / dist(M, 1)``.
    """
    M = np.asarray(xy_M, dtype=float)
    p1 = np.asarray(xy_1, dtype=float)
    p2 = np.asarray(xy_2, dtype=float)
    return np.linalg.norm(M - p2) / np.linalg.norm(M - p1)


def falling_rate_integral(X, R):
    """Trapezoidal ``∫ dX / R`` over the falling-rate period.

    Multiply by ``L_s / A`` to get the falling-rate drying time.
    Points are sorted by ``X`` internally, so input order does not matter.
    """
    X = np.asarray(X, dtype=float)
    R = np.asarray(R, dtype=float)
    order = np.argsort(X)
    return float(_trapz(1.0 / R[order], X[order]))


def entrainment_to_xy(n, r):
    """Leaching entrainment data -> right-triangle coordinates.

    ``n`` = solute fraction in the solution (kg solute / kg solution),
    ``r`` = entrainment (kg solution / kg inert solid).
    On a 1 kg solid basis: solute = r·n, solvent = r·(1−n), total = 1 + r.
    Returns ``(x_A, x_S)``; solid fraction = 1 − x_A − x_S.
    """
    n = np.asarray(n, dtype=float)
    r = np.asarray(r, dtype=float)
    tot = 1.0 + r
    return r * n / tot, r * (1 - n) / tot


def _line_intersect(p1, p2, p3, p4):
    """Intersection of lines (p1p2) and (p3p4); ``None`` if parallel."""
    (x1, y1), (x2, y2), (x3, y3), (x4, y4) = p1, p2, p3, p4
    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(den) < 1e-15:
        return None
    a = x1 * y2 - y1 * x2
    b = x3 * y4 - y3 * x4
    px = (a * (x3 - x4) - (x1 - x2) * b) / den
    py = (a * (y3 - y4) - (y1 - y2) * b) / den
    return (px, py)


def _smooth_branch(xy_points):
    """Lightly smooth a (possibly non-function) curve with a parametric spline.

    Falls back to the raw points when there are fewer than 4 or scipy is missing.
    """
    pts = np.asarray(xy_points, dtype=float)
    if len(pts) < 4:
        return pts[:, 0], pts[:, 1]
    try:
        from scipy.interpolate import splev, splprep
        tck, _ = splprep([pts[:, 0], pts[:, 1]], s=1e-4, k=3)
        xs, ys = splev(np.linspace(0, 1, 200), tck)
        return np.asarray(xs), np.asarray(ys)
    except Exception:
        return pts[:, 0], pts[:, 1]


def _light_grid(ax):
    ax.grid(True, which="major", ls=":", lw=0.6, color="#999999", alpha=0.5)
    ax.grid(True, which="minor", ls=":", lw=0.4, color="#CCCCCC", alpha=0.4)


def _triangular_grid(ax, step=0.1):
    """Light triangular grid: lines parallel to all three edges every ``step``."""
    for k in np.arange(step, 1.0, step):
        segments = [
            (ternary_to_xy(k, 1 - k, 0), ternary_to_xy(k, 0, 1 - k)),   # a = const
            (ternary_to_xy(1 - k, k, 0), ternary_to_xy(0, k, 1 - k)),   # b = const
            (ternary_to_xy(1 - k, 0, k), ternary_to_xy(0, 1 - k, k)),   # c = const
        ]
        for (xa, ya), (xb, yb) in segments:
            ax.plot([xa, xb], [ya, yb], color="#CCCCCC", lw=0.5, ls="-", zorder=0)


def _triangular_ticks(ax, step=0.2):
    """Tick labels on the three edges, each on its own family of grid lines.

    left edge AB -> A scale; bottom edge BC -> B scale; right edge AC -> C scale.
    """
    for k in np.arange(0.0, 1.0 + 1e-9, step):
        kp = int(round(k * 100))
        xa, ya = ternary_to_xy(k, 1 - k, 0)
        ax.text(xa - 0.025, ya, f"{kp}", ha="right", va="center", fontsize=7, color="#555555")
        xb, yb = ternary_to_xy(0, k, 1 - k)
        ax.text(xb, yb - 0.028, f"{kp}", ha="center", va="top", fontsize=7, color="#555555")
        xc, yc = ternary_to_xy(1 - k, 0, k)
        ax.text(xc + 0.025, yc, f"{kp}", ha="left", va="center", fontsize=7, color="#555555")


def _triangle_frame(ax, labels, fontsize=10):
    """Triangle outline, grid, edge ticks and vertex labels."""
    _triangular_grid(ax, step=0.1)
    _triangular_ticks(ax, step=0.2)
    tri = np.array([_TOP, _LEFT, _RIGHT, _TOP])
    ax.plot(tri[:, 0], tri[:, 1], color="black", lw=1.2, ls="-")
    ax.text(_TOP[0], _TOP[1] + 0.05, f"{labels[0]}\n{t('edge_left')}",
            ha="center", va="bottom", fontsize=fontsize)
    ax.text(_LEFT[0] - 0.04, _LEFT[1] - 0.04, f"{labels[1]}\n{t('edge_bottom')}",
            ha="right", va="top", fontsize=fontsize)
    ax.text(_RIGHT[0] + 0.04, _RIGHT[1] - 0.04, f"{labels[2]}\n{t('edge_right')}",
            ha="left", va="top", fontsize=fontsize)


def _coord_block(columns, rows):
    return {"columns": list(columns), "rows": [list(r) for r in rows]}


def coords_table(result):
    """Format ``result["coords"]`` from any ``plot_*`` call as a printable table."""
    out = []
    for name, block in result.get("coords", {}).items():
        out.append(f"# {name}")
        out.append("  ".join(f"{str(c):>11}" for c in block["columns"]))
        for row in block["rows"]:
            out.append("  ".join(
                f"{v:>11.4g}" if isinstance(v, (int, float, np.floating)) else f"{str(v):>11}"
                for v in row
            ))
        out.append("")
    return "\n".join(out)


# ============================================================
# 1. Ternary phase diagram (LLE)
# ============================================================

def plot_ternary(raffinate, extract, labels=("Acetone (A)", "Water (B)", "Solvent (C)"),
                 mix=None, FS=None, title=None):
    """LLE ternary diagram: binodal curve, tie lines, plait point, optional M and F–S line.

    Parameters
    ----------
    raffinate : list of (a, b, c)
        Raffinate-phase (carrier-rich) equilibrium points, wt% or fractions.
    extract : list of (a, b, c)
        Extract-phase points; row *i* pairs with ``raffinate[i]`` as one tie line.
    labels : (top, bottom-left, bottom-right) vertex names.
    mix : (a, b, c), optional — mixing point M.
    FS : ((a, b, c)_F, (a, b, c)_S), optional — draws the F–S line for the lever rule.

    Returns ``(fig, ax, result)`` with ``result["plait_xy"]`` and ``result["coords"]``.
    """
    apply_chemeng_style()
    fig, ax = plt.subplots(figsize=(6.5, 6))

    raff = np.asarray(raffinate, dtype=float)
    extr = np.asarray(extract, dtype=float)

    _triangle_frame(ax, labels)

    xr, yr = ternary_to_xy(raff[:, 0], raff[:, 1], raff[:, 2])
    xe, ye = ternary_to_xy(extr[:, 0], extr[:, 1], extr[:, 2])

    for i in range(len(xr)):
        ax.plot([xr[i], xe[i]], [yr[i], ye[i]], color="#888888", lw=0.8, ls="-", zorder=1)

    # Join both branches into one dome: raffinate (solute rising) + extract (solute falling)
    raff_sorted = raff[np.argsort(raff[:, 0])]
    extr_sorted = extr[np.argsort(extr[:, 0])]
    rb_xy = np.column_stack(ternary_to_xy(raff_sorted[:, 0], raff_sorted[:, 1], raff_sorted[:, 2]))
    eb_xy = np.column_stack(ternary_to_xy(extr_sorted[:, 0], extr_sorted[:, 1], extr_sorted[:, 2]))
    bx, by = _smooth_branch(np.vstack([rb_xy, eb_xy[::-1]]))
    ax.plot(bx, by, color="black", lw=1.4, ls="--", label=t("binodal"), zorder=2)

    ax.scatter(xr, yr, s=28, facecolors="white", edgecolors="black", marker="o",
               label=t("raffinate"), zorder=3)
    ax.scatter(xe, ye, s=28, facecolors="black", edgecolors="black", marker="^",
               label=t("extract"), zorder=3)

    # Plait point estimate = top of the binodal
    pi = int(np.argmax(by))
    ax.scatter([bx[pi]], [by[pi]], s=60, marker="*", color="black", zorder=4)
    ax.annotate("Plait point P", (bx[pi], by[pi]), textcoords="offset points",
                xytext=(8, 6), fontsize=9)

    if mix is not None:
        mx, my = ternary_to_xy(*mix)
        ax.scatter([mx], [my], s=50, marker="x", color="black", zorder=5)
        ax.annotate("M", (mx, my), textcoords="offset points", xytext=(6, -10), fontsize=10)

    if FS is not None:
        fx, fy = ternary_to_xy(*FS[0])
        sx, sy = ternary_to_xy(*FS[1])
        ax.plot([fx, sx], [fy, sy], color="#555555", lw=1.0, ls="--", label=t("fs_line"))
        ax.annotate("F", (fx, fy), textcoords="offset points", xytext=(4, 4), fontsize=9)
        ax.annotate("S", (sx, sy), textcoords="offset points", xytext=(4, 4), fontsize=9)

    ax.set_title(title or t("ternary_title"))
    ax.set_aspect("equal")
    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.12, np.sqrt(3) / 2 + 0.12)
    ax.axis("off")
    ax.legend(loc="upper right", fontsize=9)

    def _rows(pts, xs, ys, tag):
        return [[f"{tag}{i + 1}", *(100 * p / p.sum()), xs[i], ys[i]] for i, p in enumerate(pts)]

    cols = [t("point"), f"{labels[0][:6]}%", f"{labels[1][:6]}%", f"{labels[2][:6]}%", "x", "y"]
    coords = {
        "raffinate": _coord_block(cols, _rows(raff, xr, yr, "R")),
        "extract": _coord_block(cols, _rows(extr, xe, ye, "E")),
        "plait_point": _coord_block(["x", "y"], [[bx[pi], by[pi]]]),
    }
    return fig, ax, {"coords": coords, "plait_xy": (float(bx[pi]), float(by[pi]))}


# ============================================================
# 2. Distribution curve
# ============================================================

def plot_distribution(x_raff, y_ext, label=None, op_point=None,
                      title=None, xlabel=None, ylabel=None):
    """Distribution curve y vs x with a ``y = x`` reference line.

    Above ``y = x`` means K > 1 (solvent favours the solute).

    Parameters
    ----------
    x_raff, y_ext : solute content at the raffinate / extract end of each tie line.
    op_point : (x, y), optional — marks an operating point.
    """
    apply_chemeng_style()
    fig, ax = plt.subplots(figsize=(6, 6))

    x = np.asarray(x_raff, dtype=float)
    y = np.asarray(y_ext, dtype=float)
    order = np.argsort(x)
    x, y = x[order], y[order]

    lim = max(x.max(), y.max()) * 1.05
    ax.plot([0, lim], [0, lim], color="#888888", ls=":", lw=1.0, label="y = x")
    ax.plot(x, y, color="black", marker="o", ms=4, lw=1.4, label=label or t("dist_label"))

    if op_point is not None:
        ax.scatter([op_point[0]], [op_point[1]], s=60, marker="s",
                   facecolors="none", edgecolors="black", zorder=5, label=t("op_point"))

    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.set_title(title or t("dist_title"))
    ax.set_xlabel(xlabel or t("dist_x"))
    ax.set_ylabel(ylabel or t("dist_y"))
    _light_grid(ax)
    ax.legend(loc="best")

    coords = {"data": _coord_block(["x", "y"], [[xi, yi] for xi, yi in zip(x, y)])}
    return fig, ax, {"coords": coords}


# ============================================================
# 3. Countercurrent stage stepping (McCabe-Thiele construction)
# ============================================================

def plot_mccabe_thiele(eq_x, eq_y, op_xy1, op_xy2, n_stages=None,
                       title=None, xlabel=None, ylabel=None):
    """Step off theoretical stages between the equilibrium curve and the operating line.

    Starts at the feed end (``op_xy2``, larger X): horizontal to the equilibrium
    curve = one stage, then vertical back to the operating line, until the
    raffinate end (``op_xy1``) is reached.

    Parameters
    ----------
    eq_x, eq_y : equilibrium data, both increasing.
    op_xy1 : (X_lo, Y_lo) raffinate end of the operating line.
    op_xy2 : (X_hi, Y_hi) feed end of the operating line.
    n_stages : int, optional
        Given: step exactly this many stages and report ``x_final`` (rating mode).
        ``None``: step until X <= X_lo and report the stage count (design mode).

    Returns ``(fig, ax, result)`` with ``n_stages``, ``x_final``, ``stage_x``,
    ``reached_target`` and ``coords``.
    """
    apply_chemeng_style()
    fig, ax = plt.subplots(figsize=(6.5, 6))

    ex = np.asarray(eq_x, dtype=float)
    ey = np.asarray(eq_y, dtype=float)

    x_lo, y_lo = op_xy1
    x_hi, y_hi = op_xy2
    slope = (y_hi - y_lo) / (x_hi - x_lo)

    def op_line(xq):
        return y_lo + slope * (xq - x_lo)

    def eq_inv(yq):
        return float(np.interp(yq, ey, ex))

    ax.plot(ex, ey, color="black", marker="o", ms=3, lw=1.4, label=t("eq_curve"))
    ax.plot([x_lo, x_hi], [y_lo, y_hi], color="#555555", ls="--", lw=1.2, label=t("op_line"))

    stage_x = []
    corners = [(x_hi, y_hi)]
    cur_x, cur_y = x_hi, y_hi
    count = 0
    while count < _MAX_STAGES:
        nx = eq_inv(cur_y)                       # horizontal to equilibrium = one stage
        ax.plot([cur_x, nx], [cur_y, cur_y], color="black", lw=0.8, ls="-")
        count += 1
        stage_x.append(nx)
        corners.append((nx, cur_y))
        ax.annotate(str(count), (nx, cur_y), textcoords="offset points",
                    xytext=(-10, 4), fontsize=8)
        cur_x = nx
        if n_stages is not None and count >= n_stages:
            break
        if n_stages is None and cur_x <= x_lo:
            break
        ny = op_line(cur_x)                      # vertical back to the operating line
        ax.plot([cur_x, cur_x], [cur_y, ny], color="black", lw=0.8, ls="-")
        corners.append((cur_x, ny))
        cur_y = ny

    reached = stage_x[-1] <= x_lo + 1e-12
    if n_stages is None and not reached:
        warnings.warn(
            f"plot_mccabe_thiele: target X={x_lo:g} not reached after {_MAX_STAGES} stages "
            "(pinch between operating line and equilibrium curve?)",
            stacklevel=2,
        )

    ax.set_title(title or t("mt_title"))
    ax.set_xlabel(xlabel or t("mt_x"))
    ax.set_ylabel(ylabel or t("mt_y"))
    _light_grid(ax)
    ax.legend(loc="best")

    coords = {
        "equilibrium": _coord_block(["X", "Y"], [[a, b] for a, b in zip(ex, ey)]),
        "operating_line": _coord_block(["X", "Y"], [[x_lo, y_lo], [x_hi, y_hi]]),
        "stage_corners": _coord_block(["X", "Y"], [[a, b] for a, b in corners]),
    }
    result = {"n_stages": count, "x_final": stage_x[-1], "stage_x": stage_x,
              "reached_target": bool(reached), "coords": coords}
    return fig, ax, result


# ============================================================
# 4. Drying rate curve R–X
# ============================================================

def plot_drying_curve(X, R, X_c, X_2=None, R_c=None, title=None, xlabel=None, ylabel=None):
    """Drying rate curve R vs X with critical moisture ``X_c`` and the falling-rate integral.

    Parameters
    ----------
    X, R : mean free moisture and drying rate for each interval.
    X_c : critical moisture content (end of the constant-rate period).
    X_2 : optional end of the falling-rate period; if given, shades and returns
        ``result["falling_integral"] = ∫ dX / R`` from X_2 to X_c.
    R_c : optional constant rate; defaults to R interpolated at X_c.
    """
    apply_chemeng_style()
    fig, ax = plt.subplots(figsize=(6.5, 5))

    X = np.asarray(X, dtype=float)
    R = np.asarray(R, dtype=float)
    order = np.argsort(X)
    Xs, Rs = X[order], R[order]

    ax.plot(Xs, Rs, color="black", marker="o", ms=4, lw=1.4, label=t("dry_data"))

    if R_c is None:
        R_c = float(np.interp(X_c, Xs, Rs))
    ax.axvline(X_c, color="#888888", ls=":", lw=1.0)
    ax.annotate(f"$X_c$ = {X_c:g}", (X_c, R_c), textcoords="offset points",
                xytext=(6, 6), fontsize=9)
    ax.plot([X_c, Xs.max()], [R_c, R_c], color="#555555", ls="--", lw=1.0,
            label=t("dry_const").format(R_c))

    result = {"R_c": R_c}

    if X_2 is not None:
        mask = (Xs >= X_2) & (Xs <= X_c)
        fx, fr = Xs[mask], Rs[mask]
        # pad the ends by linear interpolation so the integral spans exactly [X_2, X_c]
        if fx.size == 0 or fx.min() > X_2:
            fx = np.insert(fx, 0, X_2)
            fr = np.insert(fr, 0, float(np.interp(X_2, Xs, Rs)))
        if fx.max() < X_c:
            fx = np.append(fx, X_c)
            fr = np.append(fr, R_c)
        ax.fill_between(fx, fr, color="#BBBBBB", alpha=0.5, label=t("dry_fall"))
        result["falling_integral"] = falling_rate_integral(fx, fr)

    ax.set_ylim(0, Rs.max() * 1.15)
    ax.set_title(title or t("dry_title"))
    ax.set_xlabel(xlabel or t("dry_x"))
    ax.set_ylabel(ylabel or t("dry_y"))
    _light_grid(ax)
    ax.legend(loc="best")
    check_axis_range(ax, Rs)

    result["coords"] = {
        "data": _coord_block(["X", "R"], [[a, b] for a, b in zip(Xs, Rs)]),
        "critical_point": _coord_block(["X_c", "R_c"], [[X_c, R_c]]),
    }
    return fig, ax, result


# ============================================================
# 5. Constant-pressure filtration t/V vs V
# ============================================================

def plot_filtration_linear(V, t_over_V, title=None, xlabel=None, ylabel=None):
    """Linearised constant-pressure filtration ``t/V = (K_p/2)·V + B`` with a least-squares fit.

    Returns ``result`` with ``slope`` (= K_p/2), ``intercept`` (= B) and ``K_p``.
    """
    apply_chemeng_style()
    fig, ax = plt.subplots(figsize=(6, 5))

    V = np.asarray(V, dtype=float)
    y = np.asarray(t_over_V, dtype=float)
    slope, intercept = np.polyfit(V, y, 1)

    ax.scatter(V, y, s=30, facecolors="white", edgecolors="black", marker="o",
               label=t("filt_data"), zorder=3)
    xx = np.linspace(0, V.max() * 1.05, 50)
    ax.plot(xx, slope * xx + intercept, color="black", lw=1.4,
            label=t("filt_fit").format(slope, intercept))

    ax.set_title(title or t("filt_title"))
    ax.set_xlabel(xlabel or t("filt_x"))
    ax.set_ylabel(ylabel or t("filt_y"))
    _light_grid(ax)
    ax.legend(loc="best")

    coords = {"data": _coord_block(["V", "t/V"], [[vi, yi] for vi, yi in zip(V, y)])}
    return fig, ax, {"slope": float(slope), "intercept": float(intercept),
                     "K_p": float(2 * slope), "coords": coords}


# ============================================================
# 6. Leaching (SLE) right-triangle diagram with net-flow point
# ============================================================

def plot_sle_right_triangle(under_n, under_r, feed_underflow, exit_overflow,
                            final_underflow, fresh_solvent=(0.0, 1.0),
                            labels=("Solute A", "Solvent S", "Solid B"),
                            title=None, show_construction=False):
    """Countercurrent leaching on a right-triangle diagram: net-flow point Δ + stage stepping.

    Unlike the equilateral LLE diagram: x axis = solute fraction ``x_A``,
    y axis = solvent fraction ``x_S``, solid = 1 − x_A − x_S (hypotenuse = no solid
    = overflow line). Stages alternate between the underflow line and the hypotenuse.

    Parameters
    ----------
    under_n : solute fraction in the entrained solution, increasing.
    under_r : entrainment (kg solution / kg solid) at each ``under_n``;
        pass a constant array for constant underflow.
    feed_underflow, exit_overflow, final_underflow, fresh_solvent : ``(x_A, x_S)``
        terminal compositions x_0, y_1, x_n and y_{n+1}.
    show_construction : draw every construction ray (the two terminal lines that
        fix Δ, each equilibrium ray through the solid vertex, each operating ray to Δ).

    Returns ``result`` with ``n_stages``, ``net_flow_point`` and ``coords``.
    """
    apply_chemeng_style()
    fig, ax = plt.subplots(figsize=(6.5, 6))

    under_n = np.asarray(under_n, dtype=float)
    under_r = np.asarray(under_r, dtype=float)

    V_A = (1.0, 0.0)   # 100 % solute
    V_S = (0.0, 1.0)   # 100 % solvent
    ax.plot([0, 1], [0, 0], color="black", lw=1.2, ls="-")
    ax.plot([0, 0], [0, 1], color="black", lw=1.2, ls="-")
    ax.plot([1, 0], [0, 1], color="black", lw=1.2, ls="-")
    _light_grid(ax)

    def U(n):
        r = float(np.interp(n, under_n, under_r))
        xa, xs = entrainment_to_xy(n, r)
        return (float(xa), float(xs))

    uxy = np.array([U(n) for n in np.linspace(under_n.min(), under_n.max(), 100)])
    ax.plot(uxy[:, 0], uxy[:, 1], color="black", lw=1.4, ls="--", label=t("underflow"))
    ax.plot([V_A[0], V_S[0]], [V_A[1], V_S[1]], color="#555555", lw=1.0, ls=":",
            label=t("overflow"))

    delta = _line_intersect(feed_underflow, exit_overflow, final_underflow, fresh_solvent)

    for p, name in [(feed_underflow, "$x_0$"), (exit_overflow, "$y_1$"),
                    (final_underflow, "$x_n$"), (fresh_solvent, "$y_{n+1}$")]:
        ax.scatter([p[0]], [p[1]], s=40, facecolors="white", edgecolors="black", zorder=4)
        ax.annotate(name, p, textcoords="offset points", xytext=(5, 5), fontsize=9)

    def conc(p):
        s = p[0] + p[1]
        return p[0] / s if s > 1e-12 else 0.0

    n_target = conc(final_underflow)

    if show_construction and delta is not None:
        for a, b in [(feed_underflow, exit_overflow), (final_underflow, fresh_solvent)]:
            ax.plot([a[0], delta[0]], [a[1], delta[1]], color="#999999", lw=0.6, ls="--", zorder=1)

    corners_under, corners_over = [], [exit_overflow]
    y = exit_overflow
    count = 0
    reached = False
    while count < _MAX_STAGES:
        n_y = conc(y)
        x = U(n_y)                      # equilibrium: underflow at the same concentration
        count += 1
        corners_under.append(x)
        if show_construction:
            ax.plot([0, y[0]], [0, y[1]], color="#BBBBBB", lw=0.5, ls=":", zorder=1)
        ax.plot([y[0], x[0]], [y[1], x[1]], color="black", lw=0.9, ls="-")
        ax.annotate(str(count), x, textcoords="offset points", xytext=(-12, 2), fontsize=8)
        if n_y <= n_target + 1e-9:
            reached = True
            break
        if delta is None:
            break
        ny = _line_intersect(delta, x, V_A, V_S)   # operating ray through Δ ∩ overflow line
        if ny is None:
            break
        if show_construction:
            ax.plot([x[0], delta[0]], [x[1], delta[1]], color="#999999", lw=0.5, ls="--", zorder=1)
        ax.plot([x[0], ny[0]], [x[1], ny[1]], color="black", lw=0.9, ls="-")
        corners_over.append(ny)
        y = ny

    if not reached:
        warnings.warn("plot_sle_right_triangle: stepping stopped before reaching the "
                      "final underflow concentration", stacklevel=2)

    if delta is not None:
        ax.scatter([delta[0]], [delta[1]], s=50, marker="x", color="black", zorder=5)
        ax.annotate("Δ", delta, textcoords="offset points", xytext=(6, 4), fontsize=10)

    # Expand the view to include Δ when it falls outside the triangle
    x_lo, x_hi, y_lo, y_hi = -0.05, 1.05, -0.05, 1.05
    if delta is not None and np.all(np.isfinite(delta)):
        x_lo = min(x_lo, delta[0] - 0.08)
        x_hi = max(x_hi, delta[0] + 0.08)
        y_lo = min(y_lo, delta[1] - 0.12)
        y_hi = max(y_hi, delta[1] + 0.12)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_aspect("equal")
    xspan, yspan = x_hi - x_lo, y_hi - y_lo
    if yspan > 1.6 * xspan:
        fig.set_size_inches(6.5, min(12.0, 6.0 * yspan / xspan))
    elif xspan > 1.6 * yspan:
        fig.set_size_inches(min(12.0, 6.5 * xspan / yspan), 6.0)
    ax.set_title(title or t("sle_title"))
    ax.set_xlabel(t("sle_x").format(labels[0]))
    ax.set_ylabel(t("sle_y").format(labels[1]))
    ax.legend(loc="upper right", fontsize=9)

    coords = {
        "underflow_stages": _coord_block(["x_A", "x_S"], corners_under),
        "overflow_stages": _coord_block(["x_A", "x_S"], corners_over),
        "net_flow_point": _coord_block(["x_A", "x_S"], [delta] if delta is not None else []),
    }
    return fig, ax, {"n_stages": count, "net_flow_point": delta,
                     "reached_target": reached, "coords": coords}


# ============================================================
# 7. Psychrometric chart
# ============================================================

def plot_psychrometric_chart(points=None, P=101325.0, t_range=(0.0, 50.0), title=None):
    """Psychrometric chart: saturation line, RH family, constant wet-bulb lines, state points.

    Properties come from ``psychrolib`` (ASHRAE formulations, SI units). Install
    with ``pip install chemeng-plots[psychro]``.

    Parameters
    ----------
    points : list of dict, optional
        Each gives ``Tdb`` (°C) plus one of ``RH`` (0–1), ``Twb`` (°C) or
        ``Hr`` (kg/kg). Each point is marked and annotated with Tdb, Twb, RH, H, h.
    P : total pressure in Pa.
    t_range : dry-bulb range in °C.

    Returns ``result["coords"]["states"]`` with
    ``[Tdb, Twb, RH, H, h (kJ/kg), Tdew]`` for each point.
    """
    try:
        import psychrolib
    except ImportError as e:  # pragma: no cover
        raise ImportError("plot_psychrometric_chart needs psychrolib: "
                          "pip install chemeng-plots[psychro]") from e
    psychrolib.SetUnitSystem(psychrolib.SI)

    apply_chemeng_style()
    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    T = np.linspace(t_range[0], t_range[1], 200)

    Hr_sat = np.array([psychrolib.GetHumRatioFromRelHum(tt, 1.0, P) for tt in T])
    ax.plot(T, Hr_sat, color="black", lw=1.6, ls="-", label=t("psy_sat"))

    for rh in [0.1, 0.2, 0.4, 0.6, 0.8]:
        Hr = np.array([psychrolib.GetHumRatioFromRelHum(tt, rh, P) for tt in T])
        ax.plot(T, Hr, color="#888888", lw=0.7, ls="-")
        ax.annotate(f"{int(rh * 100)}%", (T[-1], Hr[-1]), fontsize=7, color="#888888",
                    textcoords="offset points", xytext=(2, 0))

    for twb in np.arange(np.ceil(t_range[0] / 5) * 5, t_range[1], 5):
        Tline = np.linspace(twb, t_range[1], 30)
        Hr_wb = []
        for tt in Tline:
            try:
                Hr_wb.append(psychrolib.GetHumRatioFromTWetBulb(tt, twb, P))
            except Exception:
                Hr_wb.append(np.nan)
        ax.plot(Tline, Hr_wb, color="#BBBBBB", lw=0.5, ls="--")

    rows = []
    for pt in points or []:
        Tdb = pt["Tdb"]
        if "RH" in pt:
            Hr = psychrolib.GetHumRatioFromRelHum(Tdb, pt["RH"], P)
        elif "Twb" in pt:
            Hr = psychrolib.GetHumRatioFromTWetBulb(Tdb, pt["Twb"], P)
        elif "Hr" in pt:
            Hr = pt["Hr"]
        else:
            raise ValueError(f"state point {pt!r} needs one of RH, Twb or Hr")
        RH = psychrolib.GetRelHumFromHumRatio(Tdb, Hr, P)
        Twb = psychrolib.GetTWetBulbFromHumRatio(Tdb, Hr, P)
        h = psychrolib.GetMoistAirEnthalpy(Tdb, Hr) / 1000.0
        Tdew = psychrolib.GetTDewPointFromHumRatio(Tdb, Hr, P)
        ax.scatter([Tdb], [Hr], s=55, marker="o", color="black", zorder=5)
        ax.annotate(f"Tdb={Tdb}°C\nTwb={Twb:.1f}°C\nRH={RH * 100:.0f}%\n"
                    f"H={Hr * 1000:.1f} g/kg\nh={h:.1f} kJ/kg",
                    (Tdb, Hr), fontsize=7, textcoords="offset points", xytext=(6, 6),
                    bbox=dict(boxstyle="round", fc="white", ec="#888888", lw=0.5))
        rows.append([Tdb, round(Twb, 2), round(RH, 3), round(Hr, 5), round(h, 2), round(Tdew, 2)])

    ax.set_xlim(*t_range)
    ax.set_ylim(0, Hr_sat.max() * 1.05)
    ax.set_title(title or t("psy_title"))
    ax.set_xlabel(t("psy_x"))
    ax.set_ylabel(t("psy_y"))
    ax.legend(loc="upper left", fontsize=9)

    coords = {"states": _coord_block(
        ["Tdb (°C)", "Twb (°C)", "RH", "H (kg/kg)", "h (kJ/kg)", "Tdew (°C)"], rows)}
    return fig, ax, {"coords": coords}


# ============================================================
# 8. Partially miscible countercurrent LLE: difference-point stepping
# ============================================================

def _cross2d(u, v):
    return u[0] * v[1] - u[1] * v[0]


def plot_lle_d_system(raffinate, extract, F, S, RN, E1=None, S_over_F=None,
                      labels=("Solute A", "Carrier B", "Solvent C"),
                      n_max=30, title=None):
    """Countercurrent extraction with partially miscible solvent, on the equilateral diagram.

    Difference point P = line(E1, F) ∩ line(RN, S). Starting from E1, alternate a
    tie line (equilibrium) and a ray through P (operating) until RN is reached.

    Parameters
    ----------
    raffinate, extract : paired tie-line data, each point (A, B, C).
    F : feed composition. S : solvent composition. RN : final raffinate (target).
    E1 : extract leaving stage 1; if omitted it is found from ``S_over_F``
        (solvent / feed mass ratio) via the mixing point M.

    Returns ``result`` with ``n_stages``, ``difference_point`` and ``coords``.
    """
    if E1 is None and S_over_F is None:
        raise ValueError("give either E1 or S_over_F")

    apply_chemeng_style()
    fig, ax = plt.subplots(figsize=(6.5, 6))

    def norm(p):
        p = np.asarray(p, dtype=float)
        return p / p.sum()

    raff = np.array([norm(p) for p in raffinate])
    extr = np.array([norm(p) for p in extract])

    _triangle_frame(ax, labels, fontsize=9)

    xr, yr = ternary_to_xy(raff[:, 0], raff[:, 1], raff[:, 2])
    xe, ye = ternary_to_xy(extr[:, 0], extr[:, 1], extr[:, 2])
    for i in range(len(xr)):
        ax.plot([xr[i], xe[i]], [yr[i], ye[i]], color="#CCCCCC", lw=0.6, ls="-", zorder=1)
    dome = np.vstack([np.column_stack([xr, yr])[np.argsort(raff[:, 0])],
                      np.column_stack([xe, ye])[np.argsort(extr[:, 0])][::-1]])
    bx, by = _smooth_branch(dome)
    ax.plot(bx, by, color="black", lw=1.2, ls="--", label=t("binodal"), zorder=2)

    # Branches parameterised by solute fraction A
    raA, exA = raff[:, 0], extr[:, 0]
    o_r, o_e = np.argsort(raA), np.argsort(exA)
    raA_s, rxy_s = raA[o_r], np.column_stack([xr, yr])[o_r]
    exA_s, exy_s = exA[o_e], np.column_stack([xe, ye])[o_e]

    def raff_xy(xa):
        return np.array([np.interp(xa, raA_s, rxy_s[:, 0]), np.interp(xa, raA_s, rxy_s[:, 1])])

    def extr_xy(ya):
        return np.array([np.interp(ya, exA_s, exy_s[:, 0]), np.interp(ya, exA_s, exy_s[:, 1])])

    def conj_inv(ya):  # extract solute fraction -> conjugate raffinate solute fraction
        return float(np.interp(ya, exA[o_e], raA[o_e]))

    ya_grid = np.linspace(exA_s.min(), exA_s.max(), 400)

    def _ray_hit_extract(A_pt, B_pt, upper=None):
        """Where the line through A_pt and B_pt crosses the extract branch (sign change)."""
        gs = np.array([_cross2d(extr_xy(y) - A_pt, B_pt - A_pt) for y in ya_grid])
        mask = ya_grid < upper if upper is not None else np.ones_like(ya_grid, bool)
        yg, gg = ya_grid[mask], gs[mask]
        if len(yg) < 2:
            return None
        idx = np.where(np.diff(np.sign(gg)) != 0)[0]
        if len(idx) == 0:
            return None
        i0 = idx[-1]
        return yg[i0] - gg[i0] * (yg[i0 + 1] - yg[i0]) / (gg[i0 + 1] - gg[i0])

    F_xy = np.array(ternary_to_xy(*norm(F)))
    S_xy = np.array(ternary_to_xy(*norm(S)))
    RN_xy = np.array(ternary_to_xy(*norm(RN)))

    if E1 is None:
        M = (norm(F) + S_over_F * norm(S)) / (1 + S_over_F)
        M_xy = np.array(ternary_to_xy(*M))
        ya_E1 = _ray_hit_extract(RN_xy, M_xy)
        if ya_E1 is None:
            plt.close(fig)
            raise ValueError("line RN–M does not cross the extract branch; "
                             "check S_over_F and the tie-line data")
        E1_xy = extr_xy(ya_E1)
        ax.scatter([M_xy[0]], [M_xy[1]], s=30, marker="s", color="#555555", zorder=4)
        ax.annotate("M", M_xy, textcoords="offset points", xytext=(4, -10), fontsize=8)
    else:
        ya_E1 = norm(E1)[0]
        E1_xy = np.array(ternary_to_xy(*norm(E1)))

    P = _line_intersect(E1_xy, F_xy, RN_xy, S_xy)
    P = np.array(P) if P is not None else None

    for pt, name in [(F_xy, "F"), (S_xy, "S"), (E1_xy, "$E_1$"), (RN_xy, "$R_N$")]:
        ax.scatter([pt[0]], [pt[1]], s=40, facecolors="white", edgecolors="black", zorder=4)
        ax.annotate(name, pt, textcoords="offset points", xytext=(5, 4), fontsize=9)

    ya = ya_E1
    xa_target = norm(RN)[0]
    under_pts = []

    xa = conj_inv(ya)
    stages = 1
    ax.plot([extr_xy(ya)[0], raff_xy(xa)[0]], [extr_xy(ya)[1], raff_xy(xa)[1]],
            color="black", lw=0.8, ls="-")
    ax.annotate("1", raff_xy(xa), textcoords="offset points", xytext=(-12, 2), fontsize=8)
    under_pts.append(tuple(raff_xy(xa)))

    while xa > xa_target + 1e-9 and stages < n_max and P is not None:
        Rpt = raff_xy(xa)
        ya_next = _ray_hit_extract(P, Rpt, upper=ya)   # operating ray through P
        if ya_next is None:
            break
        E_next = extr_xy(ya_next)
        ax.plot([Rpt[0], E_next[0]], [Rpt[1], E_next[1]], color="black", lw=0.8, ls="-")
        ya = ya_next
        xa = conj_inv(ya)
        stages += 1
        ax.plot([E_next[0], raff_xy(xa)[0]], [E_next[1], raff_xy(xa)[1]],
                color="black", lw=0.8, ls="-")
        ax.annotate(str(stages), raff_xy(xa), textcoords="offset points",
                    xytext=(-12, 2), fontsize=8)
        under_pts.append(tuple(raff_xy(xa)))

    reached = xa <= xa_target + 1e-9
    if not reached:
        warnings.warn(
            "plot_lle_d_system: stepping stopped before reaching RN "
            f"(solute {xa:.4g} > target {xa_target:.4g}) — usually the tie-line data does not "
            "extend to low enough solute content. n_stages is then a lower bound.",
            stacklevel=2,
        )

    if P is not None:
        ax.scatter([P[0]], [P[1]], s=50, marker="x", color="black", zorder=5)
        ax.annotate("P", P, textcoords="offset points", xytext=(6, 4), fontsize=10)

    # Keep the view on the triangle; P often lies far outside.
    ax.set_xlim(-0.14, 1.14)
    ax.set_ylim(-0.14, np.sqrt(3) / 2 + 0.12)
    if P is not None and not (-0.08 <= P[0] <= 1.08 and -0.08 <= P[1] <= 0.95):
        ax.text(0.0, 1.0, t("lled_offchart").format(P[0], P[1]), transform=ax.transAxes,
                ha="left", va="top", fontsize=8, color="#555555")
    ax.set_title(title or t("lled_title"))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(loc="upper right", fontsize=9)

    coords = {"raffinate_stages": _coord_block(["x", "y"], under_pts),
              "difference_point": _coord_block(["x", "y"], [P] if P is not None else [])}
    return fig, ax, {"n_stages": stages,
                     "difference_point": tuple(float(v) for v in P) if P is not None else None,
                     "reached_target": bool(reached), "coords": coords}
