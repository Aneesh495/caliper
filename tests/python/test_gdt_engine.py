"""
tests/python/test_gdt_engine.py
Unit tests for ASME Y14.5 and ISO 1101 GD&T rule evaluation engine.
"""

import math
import unittest

from caliper.gdt_rules import (
    CylindricityCheck,
    FlatnessCheck,
    GDTSpecificationBuilder,
    MaterialCondition,
    PerpendicularityCheck,
    TruePositionCheck,
)


class TestGdtRules(unittest.TestCase):
    def test_flatness_check_pass(self):
        # Planar surface with +/- 0.015 mm variation
        pts = [(float(x), float(y), 0.01 * (x % 2)) for x in range(10) for y in range(10)]
        check = FlatnessCheck(feature_name="DatumA", tolerance_mm=0.030)
        finding = check.evaluate(pts)
        self.assertTrue(finding.passed)
        self.assertAlmostEqual(finding.measured_value, 0.01, places=5)
        self.assertEqual(finding.characteristic, "Flatness")

    def test_flatness_check_fail(self):
        pts = [(float(x), float(y), 0.15 * (x % 2)) for x in range(10) for y in range(10)]
        check = FlatnessCheck(feature_name="DatumA", tolerance_mm=0.050)
        finding = check.evaluate(pts)
        self.assertFalse(finding.passed)
        self.assertAlmostEqual(finding.measured_value, 0.15, places=5)

    def test_cylindricity_check_pass(self):
        # Cylinder radius 25.0 mm with tiny noise
        pts = []
        for h in range(10):
            z = float(h) * 5.0
            for a in range(24):
                theta = 2.0 * math.pi * float(a) / 24.0
                r = 25.0 + 0.005 * math.sin(theta * 2.0)
                pts.append((r * math.cos(theta), r * math.sin(theta), z))

        check = CylindricityCheck(feature_name="BoreSurface", tolerance_mm=0.020)
        finding = check.evaluate(pts)
        self.assertTrue(finding.passed)
        self.assertLess(finding.measured_value, 0.020)

    def test_perpendicularity_check(self):
        check = PerpendicularityCheck(
            feature_name="SideWall",
            datum_reference="DatumA",
            feature_length_mm=100.0,
            tolerance_mm=0.080,
        )
        # Normal strictly perpendicular to Datum A (Z axis) -> along X
        finding_pass = check.evaluate(feature_normal=(1.0, 0.0, 0.0), datum_normal=(0.0, 0.0, 1.0))
        self.assertTrue(finding_pass.passed)

        # Angled normal
        finding_fail = check.evaluate(feature_normal=(0.9, 0.0, 0.4358), datum_normal=(0.0, 0.0, 1.0))
        self.assertFalse(finding_fail.passed)

    def test_true_position_rfs(self):
        check = TruePositionCheck(
            feature_name="Hole_1",
            nominal_x=50.0,
            nominal_y=50.0,
            base_tolerance_mm=0.10,
            nominal_diameter=10.0,
            material_condition=MaterialCondition.RFS,
        )
        # Actual center at (50.03, 50.02) -> offset sqrt(0.0009 + 0.0004) = 0.036 mm -> diametral 0.072 mm <= 0.10
        finding = check.evaluate(actual_x=50.03, actual_y=50.02, actual_diameter=10.05)
        self.assertTrue(finding.passed)
        self.assertEqual(finding.bonus_tolerance, 0.0)

    def test_true_position_mmc_bonus(self):
        # MMC for hole (internal feature) is smallest allowable diameter (10.00 mm)
        check = TruePositionCheck(
            feature_name="Hole_MMC",
            nominal_x=50.0,
            nominal_y=50.0,
            base_tolerance_mm=0.10,
            nominal_diameter=10.0,
            mmc_diameter=10.0,
            material_condition=MaterialCondition.MMC,
        )
        # Actual hole produced at 10.08 mm gives 0.08 mm bonus tolerance!
        # Total allowable position tolerance: 0.10 + 0.08 = 0.18 mm
        finding = check.evaluate(actual_x=50.06, actual_y=50.04, actual_diameter=10.08)
        self.assertTrue(finding.passed)
        self.assertAlmostEqual(finding.bonus_tolerance, 0.08, places=4)
        self.assertAlmostEqual(finding.total_tolerance, 0.18, places=4)

    def test_gdt_builder_and_report(self):
        spec = GDTSpecificationBuilder(part_name="TestGear", drawing_number="DWG-999")
        spec.add_flatness("DatumFace", tolerance_mm=0.05)
        spec.add_cylindricity("PinBore", tolerance_mm=0.04)

        self.assertEqual(spec.part_name, "TestGear")
        self.assertEqual(len(spec.flatness_checks), 1)
        self.assertEqual(len(spec.cylindricity_checks), 1)


if __name__ == "__main__":
    unittest.main()
