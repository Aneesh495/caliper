"""
tests/python/test_geometry.py
Unit tests for PointCloud, TriangleMesh, and spatial geometry algorithms.
"""

import math
import os
import tempfile
import unittest

from caliper.geometry import PointCloud, TriangleMesh


class TestPointCloud(unittest.TestCase):
    def test_bounding_box(self):
        pts = [(0.0, 1.0, 2.0), (10.0, -5.0, 3.0), (-2.0, 4.0, 15.0)]
        cloud = PointCloud(points=pts)
        min_b, max_b = cloud.compute_bounds()
        self.assertEqual(min_b, (-2.0, -5.0, 2.0))
        self.assertEqual(max_b, (10.0, 4.0, 15.0))

    def test_empty_cloud_bounds(self):
        cloud = PointCloud(points=[])
        min_b, max_b = cloud.compute_bounds()
        self.assertEqual(min_b, (0.0, 0.0, 0.0))
        self.assertEqual(max_b, (0.0, 0.0, 0.0))

    def test_normals_estimation(self):
        # Generate planar points z = 0 with small perturbations
        pts = [(float(x), float(y), 0.0) for x in range(5) for y in range(5)]
        cloud = PointCloud(points=pts)
        cloud.estimate_normals(k=6)
        self.assertTrue(cloud.has_normals)
        self.assertEqual(len(cloud.normals), 25)
        for n in cloud.normals:
            # Planar normal should point predominantly along +/- Z
            self.assertGreater(abs(n[2]), 0.90)

    def test_ply_round_trip(self):
        pts = [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)]
        normals = [(0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0)]
        devs = [0.012, -0.045, 0.108]
        colors = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]

        cloud = PointCloud(points=pts, normals=normals, deviations=devs, colors=colors)

        with tempfile.NamedTemporaryFile(suffix=".ply", delete=False) as f:
            temp_path = f.name

        try:
            cloud.save_ply(temp_path)
            loaded = PointCloud.load_ply(temp_path)
            self.assertEqual(len(loaded), 3)
            self.assertTrue(loaded.has_normals)
            self.assertTrue(loaded.has_deviations)
            self.assertTrue(loaded.has_colors)
            for i in range(3):
                for j in range(3):
                    self.assertAlmostEqual(loaded.points[i][j], pts[i][j], places=4)
                    self.assertAlmostEqual(loaded.normals[i][j], normals[i][j], places=4)
                self.assertAlmostEqual(loaded.deviations[i], devs[i], places=4)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)


class TestTriangleMesh(unittest.TestCase):
    def test_cube_surface_area_and_volume(self):
        # 10x10x10 cube centered at origin
        # Total surface area: 6 * (10 * 10) = 600.0
        # Total volume: 10 * 10 * 10 = 1000.0
        s = 5.0
        v = [
            (-s, -s, -s), (s, -s, -s), (s, s, -s), (-s, s, -s),
            (-s, -s, s), (s, -s, s), (s, s, s), (-s, s, s)
        ]
        # 12 triangles outward wound
        triangles = [
            (v[0], v[2], v[1]), (v[0], v[3], v[2]), # bottom (-Z)
            (v[4], v[5], v[6]), (v[4], v[6], v[7]), # top (+Z)
            (v[0], v[1], v[5]), (v[0], v[5], v[4]), # front (-Y)
            (v[2], v[3], v[7]), (v[2], v[7], v[6]), # back (+Y)
            (v[0], v[4], v[7]), (v[0], v[7], v[3]), # left (-X)
            (v[1], v[2], v[6]), (v[1], v[6], v[5]), # right (+X)
        ]
        mesh = TriangleMesh(triangles=triangles)
        area = mesh.compute_surface_area()
        self.assertAlmostEqual(area, 600.0, places=3)
        vol = mesh.compute_volume()
        self.assertAlmostEqual(vol, 1000.0, places=3)

    def test_stl_binary_round_trip(self):
        t1 = ((0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (0.0, 10.0, 0.0))
        t2 = ((10.0, 0.0, 0.0), (10.0, 10.0, 0.0), (0.0, 10.0, 0.0))
        mesh = TriangleMesh(triangles=[t1, t2])

        with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
            temp_path = f.name

        try:
            mesh.save_stl_binary(temp_path)
            loaded = TriangleMesh.load_stl_binary(temp_path)
            self.assertEqual(len(loaded), 2)
            for i in range(2):
                for j in range(3):
                    for k in range(3):
                        self.assertAlmostEqual(loaded.triangles[i][j][k], mesh.triangles[i][j][k], places=4)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    def test_signed_distance_evaluation(self):
        # Horizontal plane z = 0
        t1 = ((-10.0, -10.0, 0.0), (10.0, -10.0, 0.0), (10.0, 10.0, 0.0))
        t2 = ((-10.0, -10.0, 0.0), (10.0, 10.0, 0.0), (-10.0, 10.0, 0.0))
        mesh = TriangleMesh(triangles=[t1, t2])

        query_pts = [(0.0, 0.0, 2.5), (2.0, 3.0, -1.8)]
        s_dists = mesh.evaluate_signed_distances(query_pts)
        self.assertAlmostEqual(s_dists[0], 2.5, places=4)
        self.assertAlmostEqual(s_dists[1], -1.8, places=4)


if __name__ == "__main__":
    unittest.main()
