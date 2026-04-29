"""
tests/python/test_pipeline_e2e.py
End-to-end inspection pipeline integration tests and report generator validation.
"""

import tempfile
import unittest
from pathlib import Path

from caliper.geometry import PointCloud, TriangleMesh
from caliper.pipeline import InspectionConfig, InspectionPipeline
from caliper.report import (
    export_point_deviations_csv,
    generate_json_manifest,
    generate_markdown_certificate,
    generate_svg_cross_section,
)
from caliper.viewer import generate_html_viewer


class TestPipelineE2E(unittest.TestCase):
    def setUp(self):
        # Create a horizontal reference square CAD plate: [0, 10] x [0, 10] at z = 0
        v0 = (0.0, 0.0, 0.0)
        v1 = (10.0, 0.0, 0.0)
        v2 = (10.0, 10.0, 0.0)
        v3 = (0.0, 10.0, 0.0)
        self.cad = TriangleMesh(triangles=[(v0, v1, v2), (v0, v2, v3)])

    def test_pipeline_nominal_pass(self):
        # Points directly on the plate with tiny 0.005 mm noise
        pts = [(float(x), float(y), 0.005 * ((x + y) % 2)) for x in range(1, 10) for y in range(1, 10)]
        scan = PointCloud(points=pts)

        config = InspectionConfig(tolerance_mm=0.05, sor_filter_active=False)
        pipeline = InspectionPipeline(config=config)
        result = pipeline.run(self.cad, scan, metadata={"part_id": "TEST-PLATE-01"})

        self.assertEqual(result.disposition, "PASS")
        self.assertLess(result.statistics.rmse, 0.02)
        self.assertEqual(len(result.defects), 0)

    def test_pipeline_with_burr_defect(self):
        # Points on plate plus a cluster of excess material (+0.12 mm)
        pts = [(float(x), float(y), 0.0) for x in range(1, 10) for y in range(1, 10)]
        for i in range(5):
            pts.append((5.0 + 0.1 * i, 5.0, 0.12))

        scan = PointCloud(points=pts)
        config = InspectionConfig(tolerance_mm=0.05, sor_filter_active=False)
        pipeline = InspectionPipeline(config=config)
        result = pipeline.run(self.cad, scan, metadata={"part_id": "TEST-PLATE-BURR"})

        self.assertIn(result.disposition, ["REWORK", "SCRAP"])
        self.assertGreater(len(result.defects), 0)
        self.assertIn("Flash/Burr", result.defects[0].defect_type)

    def test_reporters_and_exporters(self):
        pts = [(1.0, 1.0, 0.02), (2.0, 2.0, -0.01), (3.0, 3.0, 0.0)]
        scan = PointCloud(points=pts)
        pipeline = InspectionPipeline(tolerance_mm=0.05, run_sor=False)
        result = pipeline.run(self.cad, scan, metadata={"part_id": "REPORTER-TEST"})

        md_cert = generate_markdown_certificate(result)
        self.assertIn("Dimensional Inspection Certificate: REPORTER-TEST", md_cert)
        self.assertIn("Root Mean Square Error", md_cert)

        json_manifest = generate_json_manifest(result)
        self.assertIn('"caliper_version": "1.4.0"', json_manifest)
        self.assertIn('"standard": "ASME Y14.5-2018 / ISO 1101:2017"', json_manifest)

        svg = generate_svg_cross_section(result, plane_axis="z", plane_coord=0.0, tolerance_band=0.05)
        self.assertIn("<svg", svg)
        self.assertIn("Section Profile", svg)

        html = generate_html_viewer(result, tolerance_band=0.05)
        self.assertIn("<!DOCTYPE html>", html)
        self.assertIn("three.min.js", html)


        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_path = Path(tmp_dir) / "test_out.csv"
            export_point_deviations_csv(result, csv_path)
            self.assertTrue(csv_path.exists())
            content = csv_path.read_text(encoding="utf-8")
            self.assertIn("signed_deviation_mm", content)


if __name__ == "__main__":
    unittest.main()
