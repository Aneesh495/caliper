"""
Caliper Native Bridge & Core Mathematical Acceleration.
Module: caliper.core

Binds Python execution to compiled Futhark multicore C shared libraries via ctypes,
with high-performance pure-Python vectorized fallbacks for zero-setup execution.
"""

import os
import sys
import ctypes
import math
from typing import Tuple, List, Optional, Dict, Any

_NATIVE_LIB = None

def _try_load_native_library():
    global _NATIVE_LIB
    if _NATIVE_LIB is not None:
        return _NATIVE_LIB

    candidate_paths = [
        os.path.join(os.path.dirname(__file__), "..", "_build", "libcaliper.dylib"),
        os.path.join(os.path.dirname(__file__), "..", "_build", "libcaliper.so"),
        os.path.join(os.path.dirname(__file__), "..", "bin", "libcaliper.dylib"),
        os.path.join(os.path.dirname(__file__), "..", "bin", "libcaliper.so"),
    ]
    for path in candidate_paths:
        if os.path.exists(path):
            try:
                _NATIVE_LIB = ctypes.CDLL(path)
                return _NATIVE_LIB
            except Exception:
                continue
    return None

def point_triangle_closest_pt(p: Tuple[float, float, float],
                              tri: Tuple[Tuple[float, float, float],
                                         Tuple[float, float, float],
                                         Tuple[float, float, float]]) -> Tuple[Tuple[float, float, float], float]:
    """
    Computes exact closest point on triangle (a, b, c) to point p
    using 7-region Voronoi feature classification.
    """
    a, b, c = tri
    ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    ac = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    ap = (p[0] - a[0], p[1] - a[1], p[2] - a[2])

    d1 = ab[0]*ap[0] + ab[1]*ap[1] + ab[2]*ap[2]
    d2 = ac[0]*ap[0] + ac[1]*ap[1] + ac[2]*ap[2]
    if d1 <= 0.0 and d2 <= 0.0:
        dist = math.sqrt((p[0] - a[0])**2 + (p[1] - a[1])**2 + (p[2] - a[2])**2)
        return a, dist

    bp = (p[0] - b[0], p[1] - b[1], p[2] - b[2])
    d3 = ab[0]*bp[0] + ab[1]*bp[1] + ab[2]*bp[2]
    d4 = ac[0]*bp[0] + ac[1]*bp[1] + ac[2]*bp[2]
    if d3 >= 0.0 and d4 <= d3:
        dist = math.sqrt((p[0] - b[0])**2 + (p[1] - b[1])**2 + (p[2] - b[2])**2)
        return b, dist

    vc = d1*d4 - d3*d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        v = d1 / (d1 - d3)
        q = (a[0] + v*ab[0], a[1] + v*ab[1], a[2] + v*ab[2])
        dist = math.sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2 + (p[2] - q[2])**2)
        return q, dist

    cp = (p[0] - c[0], p[1] - c[1], p[2] - c[2])
    d5 = ab[0]*cp[0] + ab[1]*cp[1] + ab[2]*cp[2]
    d6 = ac[0]*cp[0] + ac[1]*cp[1] + ac[2]*cp[2]
    if d6 >= 0.0 and d5 <= d6:
        dist = math.sqrt((p[0] - c[0])**2 + (p[1] - c[1])**2 + (p[2] - c[2])**2)
        return c, dist

    vb = d5*d2 - d1*d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        w = d2 / (d2 - d6)
        q = (a[0] + w*ac[0], a[1] + w*ac[1], a[2] + w*ac[2])
        dist = math.sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2 + (p[2] - q[2])**2)
        return q, dist

    va = d3*d6 - d5*d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        bc = (c[0] - b[0], c[1] - b[1], c[2] - b[2])
        q = (b[0] + w*bc[0], b[1] + w*bc[1], b[2] + w*bc[2])
        dist = math.sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2 + (p[2] - q[2])**2)
        return q, dist

    denom = 1.0 / (va + vb + vc)
    v = vb * denom
    w = vc * denom
    q = (a[0] + v*ab[0] + w*ac[0], a[1] + v*ab[1] + w*ac[1], a[2] + v*ab[2] + w*ac[2])
    dist = math.sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2 + (p[2] - q[2])**2)
    return q, dist

def deviation_to_color_rgb(delta: float, tolerance: float) -> Tuple[float, float, float]:
    """
    Maps signed deviation delta to false-color RGB tuple [0.0, 1.0]:
    -tol (Blue) -> -tol/2 (Cyan) -> 0 (Green) -> +tol/2 (Yellow) -> +tol (Red).
    """
    t = max(-1.0, min(1.0, delta / tolerance))
    s = (t + 1.0) * 0.5
    if s < 0.25:
        w = s / 0.25
        return (0.0, w, 1.0)
    elif s < 0.50:
        w = (s - 0.25) / 0.25
        return (0.0, 1.0, 1.0 - w)
    elif s < 0.75:
        w = (s - 0.50) / 0.25
        return (w, 1.0, 0.0)
    else:
        w = (s - 0.75) / 0.25
        return (1.0, 1.0 - w, 0.0)
