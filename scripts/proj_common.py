# -*- coding: utf-8 -*-
"""Shared helpers for the ternary liquidus-projection panels.

Two jobs the raw tricontourf can't do on its own:
- fill the greyscale primary-phase fields all the way to the triangle frame
  (tricontourf only fills the convex hull of the data, and the edge samples sit a
  grid-step inside, so a thin white strip is left against the frame),
- place each field label in the fattest part of its own field, not at the raw
  centroid, which for a thin corner band lands on a contour or the frame.
"""
import numpy as np


def tri_xy(a, b, c):
    """Barycentric (a bottom-left, b top, c bottom-right) -> cartesian."""
    s = a + b + c
    return 0.5 * (2 * c + b) / s, (np.sqrt(3) / 2) * b / s


def augment_edges(A, B, C, prim, n=120):
    """Add edge and corner samples carrying the nearest interior field code, so a
    triangulation of the result spans the whole triangle and the fill reaches the
    frame. Returns (Xaug, Yaug, prim_aug) for the FILL only; keep isotherms on the
    original data so their positions are not pulled to the edge."""
    X, Y = tri_xy(A, B, C)
    edges = [((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),   # c = 0  (bottom-left -> top)
             ((0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),   # a = 0  (top -> bottom-right)
             ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0))]   # b = 0  (bottom-right -> bottom-left)
    ea, eb, ec = [], [], []
    for (p0, p1) in edges:
        for t in np.linspace(0.0, 1.0, n):
            ea.append(p0[0] + (p1[0] - p0[0]) * t)
            eb.append(p0[1] + (p1[1] - p0[1]) * t)
            ec.append(p0[2] + (p1[2] - p0[2]) * t)
    ea, eb, ec = np.array(ea), np.array(eb), np.array(ec)
    ex, ey = tri_xy(ea, eb, ec)
    # nearest interior sample for each edge/corner point
    idx = np.array([np.argmin((X - x) ** 2 + (Y - y) ** 2) for x, y in zip(ex, ey)])
    Xaug = np.concatenate([X, ex])
    Yaug = np.concatenate([Y, ey])
    prim_aug = np.concatenate([np.asarray(prim), np.asarray(prim)[idx]])
    return Xaug, Yaug, prim_aug


_CORNERS = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]])


def _seg_dist(px, py, a, b):
    """Distance from points (px,py) to segment a-b."""
    ax, ay = a; bx, by = b
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = np.clip(((px - ax) * dx + (py - ay) * dy) / L2, 0.0, 1.0)
    return np.hypot(px - (ax + t * dx), py - (ay + t * dy))


def label_anchor(X, Y, prim, p):
    """Anchor a field label at the point of field p that is farthest from the
    nearest obstacle, where an obstacle is any other field's sample, a triangle
    edge, or a triangle corner. This keeps a large field's label near its interior
    (the incenter), pulls a thin band's label into its widest stretch, and keeps
    every label off the frame and away from the corner element labels."""
    inside = prim == p
    xi, yi = X[inside], Y[inside]
    if xi.size == 0:
        return 0.0, 0.0
    xo, yo = X[~inside], Y[~inside]
    clear = np.full(xi.shape, np.inf)
    if xo.size:
        clear = np.array([np.min(np.hypot(xo - x, yo - y)) for x, y in zip(xi, yi)])
    edges = [(_CORNERS[0], _CORNERS[1]), (_CORNERS[1], _CORNERS[2]), (_CORNERS[2], _CORNERS[0])]
    for a, b in edges:
        clear = np.minimum(clear, _seg_dist(xi, yi, a, b))
    for cx, cy in _CORNERS:
        clear = np.minimum(clear, np.hypot(xi - cx, yi - cy))
    j = int(np.argmax(clear))
    return float(xi[j]), float(yi[j])


# triangle edges as (P0, P1) with their outward unit normals
_EDGES = [
    (_CORNERS[0], _CORNERS[1], np.array([0.0, -1.0])),                       # bottom
    (_CORNERS[1], _CORNERS[2], np.array([np.sqrt(3) / 2, 0.5])),             # right
    (_CORNERS[2], _CORNERS[0], np.array([-np.sqrt(3) / 2, 0.5])),            # left
]


def _proj_t(a, P0, P1):
    d = P1 - P0
    return float(np.clip(np.dot(a - P0, d) / np.dot(d, d), 0.0, 1.0))


def draw_field_labels(ax, X, Y, prim, present, lab, small_frac=0.06,
                      off=0.085, dmin=0.15, force_leader=(), nudge=None):
    """Label each primary-phase field. Large fields get an interior label at their
    obstacle-aware anchor. Small edge-hugging fields (fraction < small_frac), plus
    any field named in force_leader, get a leader line to a label in the outer
    margin, spread along the nearest edge so they neither stack nor clip the frame."""
    force_leader = set(force_leader)
    leaders = {0: [], 1: [], 2: []}  # per edge, entries (t, name, anchor)
    for p in present:
        sel = prim == p
        if sel.sum() < 3:
            continue
        ax_, ay_ = label_anchor(X, Y, prim, p)
        if nudge and p in nudge:            # hand offset to clear a contour label or data point
            ax_, ay_ = ax_ + nudge[p][0], ay_ + nudge[p][1]
        if sel.mean() >= small_frac and p not in force_leader:
            ax.text(ax_, ay_, lab[p], ha="center", va="center", fontsize=12, fontweight="bold",
                    bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none", alpha=0.85))
            continue
        a = np.array([ax_, ay_])
        ei = int(np.argmin([_seg_dist(np.array([ax_]), np.array([ay_]), P0, P1)[0]
                            for P0, P1, _ in _EDGES]))
        P0, P1, _ = _EDGES[ei]
        leaders[ei].append([_proj_t(a, P0, P1), lab[p], a])
    for ei, items in leaders.items():
        if not items:
            continue
        P0, P1, n = _EDGES[ei]
        items.sort(key=lambda r: r[0])
        ts = [max(0.12, min(0.88, it[0])) for it in items]
        for i in range(1, len(ts)):            # spread apart along the edge
            if ts[i] - ts[i - 1] < dmin:
                ts[i] = ts[i - 1] + dmin
        shift = max(0.0, ts[-1] - 0.88)
        ts = [t - shift for t in ts]
        ha = "center" if abs(n[0]) < 0.3 else ("left" if n[0] > 0 else "right")
        va = "center" if abs(n[1]) < 0.3 else ("bottom" if n[1] > 0 else "top")
        for t, (_, name, a) in zip(ts, items):
            edge_pt = P0 + (P1 - P0) * t
            tp = edge_pt + n * off
            ax.annotate(name, xy=a, xytext=tp, ha=ha, va=va,
                        fontsize=11, fontweight="bold",
                        arrowprops=dict(arrowstyle="-", lw=0.6, color="0.35",
                                        shrinkA=1.0, shrinkB=1.0))
