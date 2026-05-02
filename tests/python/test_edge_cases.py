"""
tests/python/test_edge_cases.py
Edge case and stress tests for degenerate meshes and singular configurations.
"""

import math
import unittest

from caliper.geometry import PointCloud, TriangleMesh
from caliper.gdt_rules import FlatnessCheck, CylindricityCheck, TruePositionCheck
from caliper.pipeline import InspectionPipeline, InspectionConfig


class TestEdgeCases(unittest.TestCase):
    def test_empty_point_cloud(self):
        cloud = PointCloud(points=[])
        cad = TriangleMesh(triangles=[((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))])
        pipeline = InspectionPipeline(tolerance_mm=0.05, run_sor=False)
        result = pipeline.run(cad, cloud)
        self.assertEqual(result.statistics.point_count, 0)
        self.assertEqual(result.statistics.rmse, 0.0)
        self.assertEqual(result.disposition, "PASS")

    def test_single_point(self):
        cloud = PointCloud(points=[(5.0, 5.0, 0.02)])
        v0 = (0.0, 0.0, 0.0)
        v1 = (10.0, 0.0, 0.0)
        v2 = (0.0, 10.0, 0.0)
        cad = TriangleMesh(triangles=[(v0, v1, v2)])
        pipeline = InspectionPipeline(tolerance_mm=0.05, run_sor=False)
        result = pipeline.run(cad, cloud)
        self.assertEqual(result.statistics.point_count, 1)
        self.assertAlmostEqual(result.statistics.rmse, 0.02, places=4)

    def test_query_point_on_vertex_and_edge(self):
        # Triangle in XY plane
        v0 = (0.0, 0.0, 0.0)
        v1 = (10.0, 0.0, 0.0)
        v2 = (0.0, 10.0, 0.0)
        cad = TriangleMesh(triangles=[(v0, v1, v2)])

        # Query exactly on vertex v0, and midpoint of edge (v0-v1)
        pts = [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0)]
        dists = cad.evaluate_signed_distances(pts)
        self.assertAlmostEqual(dists[0], 0.0, places=5)
        self.assertAlmostEqual(dists[1], 0.0, places=5)

    def test_degenerate_zero_area_triangle(self):
        # Collinear points: triangle has zero surface area
        v0 = (0.0, 0.0, 0.0)
        v1 = (5.0, 0.0, 0.0)
        v2 = (10.0, 0.0, 0.0)
        cad = TriangleMesh(triangles=[(v0, v1, v2)])
        self.assertAlmostEqual(cad.compute_surface_area(), 0.0, places=6)

        # Distance to line segment should still return orthogonal distance
        pts = [(5.0, 3.0, 0.0)]
        dists = cad.evaluate_signed_distances(pts)
        self.assertAlmostEqual(abs(dists[0]), 3.0, places=4)

    def test_needle_triangle_high_aspect_ratio(self):
        # Extremely thin needle triangle: base = 1000.0, height = 0.001
        v0 = (0.0, 0.0, 0.0)
        v1 = (1000.0, 0.0, 0.0)
        v2 = (500.0, 0.001, 0.0)
        cad = TriangleMesh(triangles=[(v0, v1, v2)])
        area = cad.compute_surface_area()
        self.assertGreater(area, 0.0)
        self.assertAlmostEqual(area, 0.5 * 1000.0 * 0.001, places=3)

    def test_extreme_coordinate_magnitudes(self):
        # Giant coordinate system: 1,000,000 mm (1 km)
        offset = 1_000_000.0
        v0 = (offset, offset, offset)
        v1 = (offset + 10.0, offset, offset)
        v2 = (offset, offset + 10.0, offset)
        cad = TriangleMesh(triangles=[(v0, v1, v2)])

        pts = [(offset + 2.0, offset + 2.0, offset + 0.05)]
        dists = cad.evaluate_signed_distances(pts)
        self.assertAlmostEqual(dists[0], 0.05, places=4)

    def test_flatness_insufficient_points(self):
        check = FlatnessCheck(feature_name="EmptyDatum", tolerance_mm=0.05)
        finding = check.evaluate([])
        self.assertTrue(finding.passed)
        self.assertEqual(finding.measured_value, 0.0)

    def test_cylindricity_insufficient_points(self):
        check = CylindricityCheck(feature_name="SmallCyl", tolerance_mm=0.05)
        finding = check.evaluate([(10.0, 0.0, 0.0), (0.0, 10.0, 0.0)])
        self.assertTrue(finding.passed)
        self.assertEqual(finding.measured_value, 0.0)

    def test_sor_isolated_outlier_removal(self):
        # 100 closely clustered points around origin
        pts = [(float(x) * 0.1, float(y) * 0.1, 0.0) for x in range(10) for y in range(10)]
        # Add 1 massive outlier far away
        pts.append((500.0, 500.0, 500.0))
        cloud = PointCloud(points=pts)
        filtered = cloud.filter_sor(k=8, alpha=1.5)
        self.assertEqual(len(filtered), 100)
        self.assertNotIn((500.0, 500.0, 500.0), filtered.points)


if __name__ == "__main__":
    unittest.main()
