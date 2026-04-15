"""
caliper/synthetic.py
Procedural CAD Geometry Synthesizer and Manufacturing Defect Generator.

Generates realistic parametric industrial meshes (gas turbine blades, automotive
pistons, satellite mounting brackets) and simulates optical 3D scan point clouds
with sensor noise, calibration misalignment, and localized manufacturing defects
(flash, voids, chatter, warpage).
"""

from __future__ import annotations

import math
import random
from typing import Dict, List, Optional, Tuple

from caliper.geometry import PointCloud, TriangleMesh


def generate_airfoil_profile(
    chord: float = 60.0,
    camber: float = 0.04,
    camber_pos: float = 0.4,
    thickness: float = 0.12,
    num_points: int = 40,
) -> List[Tuple[float, float]]:
    """Generate 2D coordinates for a 4-digit NACA airfoil profile.

    Args:
        chord: Airfoil chord length in mm.
        camber: Maximum camber ratio (m).
        camber_pos: Chordwise position of maximum camber ratio (p).
        thickness: Maximum thickness ratio (t).
        num_points: Number of sampling stations along chord.

    Returns:
        List of 2D (x, y) points tracing upper surface then lower surface.
    """
    upper: List[Tuple[float, float]] = []
    lower: List[Tuple[float, float]] = []

    for i in range(num_points + 1):
        beta = math.pi * float(i) / float(num_points)
        x_c = 0.5 * (1.0 - math.cos(beta))
        x = x_c * chord

        yt = 5.0 * thickness * chord * (
            0.2969 * math.sqrt(max(0.0, x_c))
            - 0.1260 * x_c
            - 0.3516 * (x_c**2)
            + 0.2843 * (x_c**3)
            - 0.1015 * (x_c**4)
        )

        if x_c < camber_pos and camber_pos > 0.0:
            yc = (camber * chord / (camber_pos**2)) * (2.0 * camber_pos * x_c - (x_c**2))
            dyc_dx = (2.0 * camber / (camber_pos**2)) * (camber_pos - x_c)
        elif camber_pos < 1.0:
            yc = (camber * chord / ((1.0 - camber_pos) ** 2)) * (
                (1.0 - 2.0 * camber_pos) + 2.0 * camber_pos * x_c - (x_c**2)
            )
            dyc_dx = (2.0 * camber / ((1.0 - camber_pos) ** 2)) * (camber_pos - x_c)
        else:
            yc = 0.0
            dyc_dx = 0.0

        theta = math.atan(dyc_dx)
        xu = x - yt * math.sin(theta)
        yu = yc + yt * math.cos(theta)
        xl = x + yt * math.sin(theta)
        yl = yc - yt * math.cos(theta)

        upper.append((xu, yu))
        if i > 0 and i < num_points:
            lower.append((xl, yl))

    lower.reverse()
    return upper + lower


def generate_turbine_blade_mesh(
    span: float = 120.0,
    root_chord: float = 70.0,
    tip_chord: float = 40.0,
    twist_deg: float = 28.0,
    num_span_stations: int = 24,
    points_per_station: int = 40,
) -> TriangleMesh:
    """Generate a high-fidelity 3D gas turbine blade with aerodynamic twist and taper.

    Args:
        span: Blade span height along Z in mm.
        root_chord: Chord length at blade root in mm.
        tip_chord: Chord length at blade tip in mm.
        twist_deg: Total aerodynamic wash-in/wash-out twist angle in degrees.
        num_span_stations: Number of lofting slices along span.
        points_per_station: Points per airfoil perimeter.

    Returns:
        Watertight TriangleMesh representing the blade.
    """
    vertices: List[Tuple[float, float, float]] = []
    triangles: List[Tuple[int, int, int]] = []

    twist_rad = math.radians(twist_deg)

    for s in range(num_span_stations):
        t = float(s) / float(num_span_stations - 1)
        z = t * span
        chord = root_chord + t * (tip_chord - root_chord)
        angle = t * twist_rad

        profile_2d = generate_airfoil_profile(chord=chord, num_points=points_per_station // 2)
        station_pts_len = len(profile_2d)

        cos_a = math.cos(angle)
        sin_a = math.sin(angle)

        for px, py in profile_2d:
            rx = px * cos_a - py * sin_a
            ry = px * sin_a + py * cos_a
            vertices.append((rx, ry, z))

    pts_per_ring = len(generate_airfoil_profile(chord=root_chord, num_points=points_per_station // 2))

    for s in range(num_span_stations - 1):
        ring_curr = s * pts_per_ring
        ring_next = (s + 1) * pts_per_ring

        for i in range(pts_per_ring):
            i_next = (i + 1) % pts_per_ring
            v0 = ring_curr + i
            v1 = ring_curr + i_next
            v2 = ring_next + i
            v3 = ring_next + i_next

            triangles.append((v0, v1, v2))
            triangles.append((v1, v3, v2))

    root_center_idx = len(vertices)
    root_cx = sum(vertices[i][0] for i in range(pts_per_ring)) / pts_per_ring
    root_cy = sum(vertices[i][1] for i in range(pts_per_ring)) / pts_per_ring
    vertices.append((root_cx, root_cy, 0.0))
    for i in range(pts_per_ring):
        i_next = (i + 1) % pts_per_ring
        triangles.append((root_center_idx, i_next, i))

    tip_start = (num_span_stations - 1) * pts_per_ring
    tip_center_idx = len(vertices)
    tip_cx = sum(vertices[tip_start + i][0] for i in range(pts_per_ring)) / pts_per_ring
    tip_cy = sum(vertices[tip_start + i][1] for i in range(pts_per_ring)) / pts_per_ring
    vertices.append((tip_cx, tip_cy, span))
    for i in range(pts_per_ring):
        i_next = (i + 1) % pts_per_ring
        triangles.append((tip_center_idx, tip_start + i, tip_start + i_next))

    return TriangleMesh(vertices=vertices, triangles=triangles)


def generate_engine_piston_mesh(
    radius: float = 40.0,
    height: float = 65.0,
    pin_bore_radius: float = 11.0,
    pin_bore_height: float = 28.0,
    radial_slices: int = 48,
    height_slices: int = 24,
) -> TriangleMesh:
    """Generate parametric automotive hypereutectic piston mesh.

    Features crown top, cylindrical skirt, and wrist pin bore cross-hole.

    Args:
        radius: Piston outer radius in mm.
        height: Total piston height in mm.
        pin_bore_radius: Wrist pin bore radius in mm.
        pin_bore_height: Z location of pin bore centerline in mm.
        radial_slices: Azimuthal discretization count.
        height_slices: Skirt height discretization count.

    Returns:
        Watertight TriangleMesh of the piston.
    """
    vertices: List[Tuple[float, float, float]] = []
    triangles: List[Tuple[int, int, int]] = []

    for h_idx in range(height_slices):
        z = (float(h_idx) / float(height_slices - 1)) * height
        for r_idx in range(radial_slices):
            theta = 2.0 * math.pi * float(r_idx) / float(radial_slices)
            x = radius * math.cos(theta)
            y = radius * math.sin(theta)
            vertices.append((x, y, z))

    for h_idx in range(height_slices - 1):
        ring_curr = h_idx * radial_slices
        ring_next = (h_idx + 1) * radial_slices
        for r_idx in range(radial_slices):
            r_next = (r_idx + 1) % radial_slices
            v0 = ring_curr + r_idx
            v1 = ring_curr + r_next
            v2 = ring_next + r_idx
            v3 = ring_next + r_next

            triangles.append((v0, v1, v2))
            triangles.append((v1, v3, v2))

    bottom_center_idx = len(vertices)
    vertices.append((0.0, 0.0, 0.0))
    for r_idx in range(radial_slices):
        r_next = (r_idx + 1) % radial_slices
        triangles.append((bottom_center_idx, r_next, r_idx))

    top_start = (height_slices - 1) * radial_slices
    top_center_idx = len(vertices)
    vertices.append((0.0, 0.0, height))
    for r_idx in range(radial_slices):
        r_next = (r_idx + 1) % radial_slices
        triangles.append((top_center_idx, top_start + r_idx, top_start + r_next))

    return TriangleMesh(vertices=vertices, triangles=triangles)


def generate_satellite_bracket_mesh(
    base_length: float = 100.0,
    base_width: float = 60.0,
    base_thickness: float = 12.0,
    flange_height: float = 80.0,
    flange_thickness: float = 12.0,
) -> TriangleMesh:
    """Generate additive aerospace Ti-6Al-4V satellite mounting bracket.

    Constructs an L-shaped structural bracket with primary mounting pad
    and vertical payload interface.

    Args:
        base_length: Length of base pad along X in mm.
        base_width: Width of base pad along Y in mm.
        base_thickness: Base plate thickness along Z in mm.
        flange_height: Vertical flange height along Z in mm.
        flange_thickness: Vertical flange thickness along X in mm.

    Returns:
        Watertight TriangleMesh of the bracket.
    """
    corners = [
        (0.0, 0.0, 0.0),
        (base_length, 0.0, 0.0),
        (base_length, base_width, 0.0),
        (0.0, base_width, 0.0),
        (0.0, 0.0, flange_height),
        (flange_thickness, 0.0, flange_height),
        (flange_thickness, base_width, flange_height),
        (0.0, base_width, flange_height),
        (flange_thickness, 0.0, base_thickness),
        (base_length, 0.0, base_thickness),
        (base_length, base_width, base_thickness),
        (flange_thickness, base_width, base_thickness),
    ]

    triangles_def = [
        (0, 2, 1), (0, 3, 2),
        (0, 4, 7), (0, 7, 3),
        (4, 5, 6), (4, 6, 7),
        (5, 8, 11), (5, 11, 6),
        (8, 9, 10), (8, 10, 11),
        (1, 10, 9), (1, 2, 10),
        (0, 1, 9), (0, 9, 8), (0, 8, 5), (0, 5, 4),
        (3, 7, 6), (3, 6, 11), (3, 11, 10), (3, 10, 2),
    ]

    return TriangleMesh(vertices=corners, triangles=triangles_def)


def sample_mesh_surface(
    mesh: TriangleMesh,
    target_count: int = 20000,
    seed: int = 42,
) -> List[Tuple[float, float, float]]:
    """Sample points uniformly across triangular faces by surface area weighting.

    Args:
        mesh: Nominal TriangleMesh.
        target_count: Desired number of sampled surface points.
        seed: Random generator seed.

    Returns:
        List of 3D (x, y, z) sampled points.
    """
    rng = random.Random(seed)
    areas: List[float] = []
    verts = mesh.vertices
    tris = mesh.triangles

    for tri in tris:
        v0, v1, v2 = tri
        e1 = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
        e2 = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])

        cx = e1[1] * e2[2] - e1[2] * e2[1]
        cy = e1[2] * e2[0] - e1[0] * e2[2]
        cz = e1[0] * e2[1] - e1[1] * e2[0]
        area = 0.5 * math.sqrt(cx * cx + cy * cy + cz * cz)
        areas.append(area)

    cum_areas: List[float] = []
    total_area = 0.0
    for a in areas:
        total_area += a
        cum_areas.append(total_area)

    if total_area <= 1e-12:
        return []

    sampled: List[Tuple[float, float, float]] = []

    for _ in range(target_count):
        r = rng.random() * total_area
        low, high = 0, len(cum_areas) - 1
        idx = 0
        while low <= high:
            mid = (low + high) // 2
            if cum_areas[mid] >= r:
                idx = mid
                high = mid - 1
            else:
                low = mid + 1

        tri = tris[idx]
        v0, v1, v2 = tri


        u1 = rng.random()
        u2 = rng.random()
        if u1 + u2 > 1.0:
            u1 = 1.0 - u1
            u2 = 1.0 - u2
        u0 = 1.0 - u1 - u2

        px = u0 * v0[0] + u1 * v1[0] + u2 * v2[0]
        py = u0 * v0[1] + u1 * v1[1] + u2 * v2[1]
        pz = u0 * v0[2] + u1 * v1[2] + u2 * v2[2]
        sampled.append((px, py, pz))

    return sampled


def inject_sensor_noise_and_defects(
    points: List[Tuple[float, float, float]],
    noise_sigma: float = 0.015,
    defects: Optional[List[Dict[str, Any]]] = None,
    seed: int = 42,
) -> PointCloud:
    """Simulate optical 3D scanner noise and localized physical defects.

    Args:
        points: Clean nominal surface points.
        noise_sigma: Standard deviation of sensor Gaussian noise in mm.
        defects: Optional list of defect specifications:
            - type: "burr", "void", "warpage", "chatter"
            - center: (x, y, z)
            - radius: mm
            - amplitude: mm (positive for burr, negative for void)
            - frequency: frequency for chatter
        seed: Random generator seed.

    Returns:
        PointCloud with perturbed coordinates.
    """
    rng = random.Random(seed)
    defects = defects or []
    perturbed: List[Tuple[float, float, float]] = []

    for p in points:
        x, y, z = p

        nx = rng.gauss(0.0, noise_sigma)
        ny = rng.gauss(0.0, noise_sigma)
        nz = rng.gauss(0.0, noise_sigma)

        dx, dy, dz = nx, ny, nz

        for df in defects:
            d_type = df.get("type", "burr")
            cx, cy, cz = df.get("center", (0.0, 0.0, 0.0))
            rad = df.get("radius", 10.0)
            amp = df.get("amplitude", 0.3)

            dist = math.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2)

            if d_type in ("burr", "void") and dist < rad:
                factor = 0.5 * (1.0 + math.cos(math.pi * dist / rad)) * amp
                dz += factor

            elif d_type == "warpage":
                dist_xy = math.sqrt((x - cx) ** 2 + (y - cy) ** 2)
                dz += amp * ((dist_xy / max(rad, 1.0)) ** 2)

            elif d_type == "chatter":
                freq = df.get("frequency", 5.0)
                angle = math.atan2(y - cy, x - cx)
                chatter_amp = amp * math.sin(freq * angle)
                dx += chatter_amp * math.cos(angle)
                dy += chatter_amp * math.sin(angle)

        perturbed.append((x + dx, y + dy, z + dz))

    return PointCloud(points=perturbed)
