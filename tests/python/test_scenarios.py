"""
tests/python/test_scenarios.py
Integration tests executing procedural industrial metrology scenarios.
"""

import tempfile
import unittest
from pathlib import Path

from caliper.scenarios.engine_piston import run_engine_piston_scenario
from caliper.scenarios.satellite_bracket import run_satellite_bracket_scenario
from caliper.scenarios.turbine_blade import run_turbine_blade_scenario


class TestScenarios(unittest.TestCase):
    def test_turbine_blade_scenario_execution(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_p = Path(tmp_dir) / "turbine"
            res = run_turbine_blade_scenario(output_dir=out_p, num_points=250, inject_defects=True)

            self.assertIsNotNone(res)
            self.assertIn(res.disposition, ["PASS", "REWORK", "SCRAP"])
            self.assertTrue((out_p / "dashboard.html").exists())
            self.assertTrue((out_p / "cross_section_z50.svg").exists())
            self.assertTrue((out_p / "certificate.md").exists())
            self.assertTrue((out_p / "manifest.json").exists())
            self.assertTrue((out_p / "points_log.csv").exists())

    def test_engine_piston_scenario_execution(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_p = Path(tmp_dir) / "piston"
            res = run_engine_piston_scenario(output_dir=out_p, num_points=250, inject_defects=True)

            self.assertIsNotNone(res)
            self.assertIn(res.disposition, ["PASS", "REWORK", "SCRAP"])
            self.assertTrue((out_p / "dashboard.html").exists())
            self.assertTrue((out_p / "cross_section_z65.svg").exists())
            self.assertTrue((out_p / "certificate.md").exists())

    def test_satellite_bracket_scenario_execution(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_p = Path(tmp_dir) / "bracket"
            res = run_satellite_bracket_scenario(output_dir=out_p, num_points=250, inject_defects=True)

            self.assertIsNotNone(res)
            self.assertIn(res.disposition, ["PASS", "REWORK", "SCRAP"])
            self.assertTrue((out_p / "dashboard.html").exists())
            self.assertTrue((out_p / "cross_section_y30.svg").exists())
            self.assertTrue((out_p / "certificate.md").exists())


if __name__ == "__main__":
    unittest.main()
