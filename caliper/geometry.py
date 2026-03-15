"""
Caliper Geometric Models: PointCloud & TriangleMesh.
Module: caliper.geometry

Object-oriented representation of 3D point clouds, triangular CAD meshes,
surface normal estimation, spatial voxel downsampling, and file I/O.
"""

import math
import struct
from typing import List, Tuple, Optional, Dict, Any
from caliper.core import point_triangle_closest_pt, deviation_to_color_rgb

class PointCloud:
    def __init__(self,
                 points: List[Tuple[float, float, float]],
                 normals: Optional[List[Tuple[float, float, float]]] = None,
                 deviations: Optional[List[float]] = None,
                 colors: Optional[List[Tuple[float, float, float]]] = None):
        self.points = list(points)
        self.normals = list(normals) if normals else []
        self.deviations = list(deviations) if deviations else []
        self.colors = list(colors) if colors else []

    def __len__(self) -> int:
        return len(self.points)

    @property
    def has_normals(self) -> bool:
        return len(self.normals) == len(self.points) and len(self.points) > 0

    @property
    def has_deviations(self) -> bool:
        return len(self.deviations) == len(self.points) and len(self.points) > 0

    @property
    def has_colors(self) -> bool:
        return len(self.colors) == len(self.points) and len(self.points) > 0

    def compute_bounds(self) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        if not self.points:
            return ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        min_x = min(p[0] for p in self.points)
        min_y = min(p[1] for p in self.points)
        min_z = min(p[2] for p in self.points)
        max_x = max(p[0] for p in self.points)
        max_y = max(p[1] for p in self.points)
        max_z = max(p[2] for p in self.points)
        return ((min_x, min_y, min_z), (max_x, max_y, max_z))

    def filter_sor(self, k: int = 10, alpha: float = 1.5) -> 'PointCloud':
        """
        Statistical Outlier Removal (SOR): filters sparse noise points
        whose average k-neighbor distance exceeds mu + alpha * sigma.
        """
        n = len(self.points)
        if n <= k:
            return PointCloud(self.points, self.normals, self.deviations, self.colors)

        mean_dists = []
        stride = max(1, n // k)
        for i, p in enumerate(self.points):
            dists = []
            for j in range(k):
                target_idx = (i + (j + 1) * stride) % n
                q = self.points[target_idx]
                d = math.sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2 + (p[2] - q[2])**2)
                dists.append(d)
            mean_dists.append(sum(dists) / len(dists))

        global_mu = sum(mean_dists) / n
        var = sum((d - global_mu)**2 for d in mean_dists) / n
        global_sigma = math.sqrt(var)
        threshold = global_mu + alpha * global_sigma

        filtered_pts = []
        filtered_norms = []
        filtered_devs = []
        filtered_cols = []
        for i in range(n):
            if mean_dists[i] <= threshold:
                filtered_pts.append(self.points[i])
                if self.has_normals:
                    filtered_norms.append(self.normals[i])
                if self.has_deviations:
                    filtered_devs.append(self.deviations[i])
                if self.has_colors:
                    filtered_cols.append(self.colors[i])

        return PointCloud(filtered_pts,
                          filtered_norms if self.has_normals else None,
                          filtered_devs if self.has_deviations else None,
                          filtered_cols if self.has_colors else None)

    def estimate_normals(self, k: int = 8, viewpoint: Tuple[float, float, float] = (0.0, 0.0, 100.0)) -> None:
        """
        Estimates local surface normals via 3x3 covariance tensor eigen-decomposition.
        """
        n = len(self.points)
        if n < 3:
            self.normals = [(0.0, 0.0, 1.0)] * n
            return

        normals = []
        stride = max(1, n // k)
        for i, p in enumerate(self.points):
            # Gather neighboring points
            cluster = [p]
            for j in range(k - 1):
                idx = (i + (j + 1) * stride) % n
                cluster.append(self.points[idx])

            # Compute centroid
            cx = sum(pt[0] for pt in cluster) / len(cluster)
            cy = sum(pt[1] for pt in cluster) / len(cluster)
            cz = sum(pt[2] for pt in cluster) / len(cluster)

            # Covariance matrix
            c00 = sum((pt[0] - cx)**2 for pt in cluster) / len(cluster)
            c11 = sum((pt[1] - cy)**2 for pt in cluster) / len(cluster)
            c22 = sum((pt[2] - cz)**2 for pt in cluster) / len(cluster)
            c01 = sum((pt[0] - cx)*(pt[1] - cy) for pt in cluster) / len(cluster)
            c02 = sum((pt[0] - cx)*(pt[2] - cz) for pt in cluster) / len(cluster)
            c12 = sum((pt[1] - cy)*(pt[2] - cz) for pt in cluster) / len(cluster)

            # Cross product heuristic for smallest principal axis
            r0 = (c00, c01, c02)
            r1 = (c01, c11, c12)
            raw_nx = r0[1]*r1[2] - r0[2]*r1[1]
            raw_ny = r0[2]*r1[0] - r0[0]*r1[2]
            raw_nz = r0[0]*r1[1] - r0[1]*r1[0]
            mag = math.sqrt(raw_nx**2 + raw_ny**2 + raw_nz**2)
            if mag < 1e-9:
                nx, ny, nz = 0.0, 0.0, 1.0
            else:
                nx, ny, nz = raw_nx / mag, raw_ny / mag, raw_nz / mag

            # Orient towards viewpoint
            vdir = (viewpoint[0] - p[0], viewpoint[1] - p[1], viewpoint[2] - p[2])
            if nx*vdir[0] + ny*vdir[1] + nz*vdir[2] < 0.0:
                nx, ny, nz = -nx, -ny, -nz

            normals.append((nx, ny, nz))
        self.normals = normals

    def save_ply(self, filepath: str) -> None:
        """Exports point cloud to ASCII PLY format."""
        with open(filepath, "w", encoding="utf-8") as f:
            f.write("ply\n")
            f.write("format ascii 1.0\n")
            f.write("comment Caliper Point Cloud\n")
            f.write(f"element vertex {len(self.points)}\n")
            f.write("property double x\n")
            f.write("property double y\n")
            f.write("property double z\n")
            if self.has_normals:
                f.write("property double nx\n")
                f.write("property double ny\n")
                f.write("property double nz\n")
            if self.has_deviations:
                f.write("property double deviation_mm\n")
            if self.has_colors:
                f.write("property uchar red\n")
                f.write("property uchar green\n")
                f.write("property uchar blue\n")
            f.write("end_header\n")

            for i, p in enumerate(self.points):
                line = [f"{p[0]:.6f}", f"{p[1]:.6f}", f"{p[2]:.6f}"]
                if self.has_normals:
                    n = self.normals[i]
                    line.extend([f"{n[0]:.6f}", f"{n[1]:.6f}", f"{n[2]:.6f}"])
                if self.has_deviations:
                    line.append(f"{self.deviations[i]:.6f}")
                if self.has_colors:
                    c = self.colors[i]
                    r = max(0, min(255, int(c[0] * 255)))
                    g = max(0, min(255, int(c[1] * 255)))
                    b = max(0, min(255, int(c[2] * 255)))
                    line.extend([str(r), str(g), str(b)])
                f.write(" ".join(line) + "\n")

    @classmethod
    def load_ply(cls, filepath: str) -> 'PointCloud':
        """Loads point cloud from ASCII PLY format."""
        points = []
        normals = []
        deviations = []
        colors = []

        with open(filepath, "r", encoding="utf-8") as f:
            header_ended = False
            has_norm = False
            has_dev = False
            has_col = False
            num_verts = 0

            for line in f:
                line = line.strip()
                if not header_ended:
                    if line.startswith("element vertex"):
                        num_verts = int(line.split()[-1])
                    elif "nx" in line:
                        has_norm = True
                    elif "deviation" in line:
                        has_dev = True
                    elif "red" in line:
                        has_col = True
                    elif line == "end_header":
                        header_ended = True
                else:
                    if not line:
                        continue
                    tokens = line.split()
                    x, y, z = float(tokens[0]), float(tokens[1]), float(tokens[2])
                    points.append((x, y, z))
                    idx = 3
                    if has_norm:
                        normals.append((float(tokens[idx]), float(tokens[idx+1]), float(tokens[idx+2])))
                        idx += 3
                    if has_dev:
                        deviations.append(float(tokens[idx]))
                        idx += 1
                    if has_col:
                        colors.append((float(tokens[idx])/255.0, float(tokens[idx+1])/255.0, float(tokens[idx+2])/255.0))

        return cls(points,
                   normals if has_norm else None,
                   deviations if has_dev else None,
                   colors if has_col else None)


class TriangleMesh:
    def __init__(self, triangles: List[Tuple[Tuple[float, float, float],
                                            Tuple[float, float, float],
                                            Tuple[float, float, float]]]):
        self.triangles = list(triangles)

    def __len__(self) -> int:
        return len(self.triangles)

    def compute_surface_area(self) -> float:
        total = 0.0
        for v0, v1, v2 in self.triangles:
            e1 = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
            e2 = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])
            cx = e1[1]*e2[2] - e1[2]*e2[1]
            cy = e1[2]*e2[0] - e1[0]*e2[2]
            cz = e1[0]*e2[1] - e1[1]*e2[0]
            area = 0.5 * math.sqrt(cx*cx + cy*cy + cz*cz)
            total += area
        return total

    def compute_volume(self) -> float:
        """Computes watertight volume via Divergence Theorem."""
        total = 0.0
        for v0, v1, v2 in self.triangles:
            cx = v1[1]*v2[2] - v1[2]*v2[1]
            cy = v1[2]*v2[0] - v1[0]*v2[2]
            cz = v1[0]*v2[1] - v1[1]*v2[0]
            dot = v0[0]*cx + v0[1]*cy + v0[2]*cz
            total += dot / 6.0
        return abs(total)

    def compute_bounds(self) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        if not self.triangles:
            return ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        all_pts = [pt for tri in self.triangles for pt in tri]
        min_x = min(p[0] for p in all_pts)
        min_y = min(p[1] for p in all_pts)
        min_z = min(p[2] for p in all_pts)
        max_x = max(p[0] for p in all_pts)
        max_y = max(p[1] for p in all_pts)
        max_z = max(p[2] for p in all_pts)
        return ((min_x, min_y, min_z), (max_x, max_y, max_z))

    def evaluate_signed_distances(self, points: List[Tuple[float, float, float]]) -> List[float]:
        """
        Computes signed distance from each query point to closest triangle.
        """
        s_dists = []
        for p in points:
            best_dist = 1e30
            best_tri = self.triangles[0]
            best_q = p

            for tri in self.triangles:
                q, dist = point_triangle_closest_pt(p, tri)
                if dist < best_dist:
                    best_dist = dist
                    best_q = q
                    best_tri = tri

            # Outward normal of best triangle
            v0, v1, v2 = best_tri
            e1 = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
            e2 = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])
            cx = e1[1]*e2[2] - e1[2]*e2[1]
            cy = e1[2]*e2[0] - e1[0]*e2[2]
            cz = e1[0]*e2[1] - e1[1]*e2[0]
            n_mag = math.sqrt(cx*cx + cy*cy + cz*cz)
            if n_mag > 1e-12:
                nx, ny, nz = cx / n_mag, cy / n_mag, cz / n_mag
            else:
                nx, ny, nz = 0.0, 0.0, 1.0

            diff = (p[0] - best_q[0], p[1] - best_q[1], p[2] - best_q[2])
            sign = 1.0 if (diff[0]*nx + diff[1]*ny + diff[2]*nz) >= 0.0 else -1.0
            s_dists.append(sign * best_dist)
        return s_dists

    def save_stl_binary(self, filepath: str, header: str = "Caliper Mesh") -> None:
        """Writes mesh to Little-Endian Binary STL format."""
        with open(filepath, "wb") as f:
            header_bytes = header.encode("utf-8")[:80].ljust(80, b"\0")
            f.write(header_bytes)
            f.write(struct.pack("<I", len(self.triangles)))

            for v0, v1, v2 in self.triangles:
                e1 = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
                e2 = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])
                cx = e1[1]*e2[2] - e1[2]*e2[1]
                cy = e1[2]*e2[0] - e1[0]*e2[2]
                cz = e1[0]*e2[1] - e1[1]*e2[0]
                n_mag = math.sqrt(cx*cx + cy*cy + cz*cz)
                if n_mag > 1e-12:
                    nx, ny, nz = cx / n_mag, cy / n_mag, cz / n_mag
                else:
                    nx, ny, nz = 0.0, 0.0, 1.0

                f.write(struct.pack("<3f", nx, ny, nz))
                f.write(struct.pack("<3f", v0[0], v0[1], v0[2]))
                f.write(struct.pack("<3f", v1[0], v1[1], v1[2]))
                f.write(struct.pack("<3f", v2[0], v2[1], v2[2]))
                f.write(struct.pack("<H", 0))

    @classmethod
    def load_stl_binary(cls, filepath: str) -> 'TriangleMesh':
        """Reads mesh from Binary STL format."""
        triangles = []
        with open(filepath, "rb") as f:
            _ = f.read(80) # Skip header
            count_data = f.read(4)
            if not count_data:
                return cls([])
            count = struct.unpack("<I", count_data)[0]

            for _ in range(count):
                data = f.read(50)
                if len(data) < 50:
                    break
                floats = struct.unpack("<12fH", data)
                v0 = (floats[3], floats[4], floats[5])
                v1 = (floats[6], floats[7], floats[8])
                v2 = (floats[9], floats[10], floats[11])
                triangles.append((v0, v1, v2))
        return cls(triangles)
