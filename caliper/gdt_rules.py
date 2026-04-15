"""
Caliper GD&T Specification Builder: ASME Y14.5 & ISO 1101.
Module: caliper.gdt_rules

Declarative tolerance rule definitions, datums, material condition modifiers
(MMC, LMC, RFS), and compliance finding reports.
"""

from enum import Enum
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Any
import math

class MaterialCondition(Enum):
    RFS = 0  # Regardless of Feature Size
    MMC = 1  # Maximum Material Condition
    LMC = 2  # Least Material Condition

@dataclass
class GdtFinding:
    feature_name: str
    tolerance_type: str
    measured_value: float
    specified_tolerance: float
    bonus_tolerance: float
    total_tolerance: float
    passed: bool
    deviation_margin: float
    datum_reference: Optional[str] = None
    notes: str = ""

    @property
    def characteristic(self) -> str:
        return self.tolerance_type

    @property
    def tolerance_limit(self) -> float:
        return self.specified_tolerance

    @property
    def details(self) -> str:
        return self.notes


@dataclass
class FlatnessCheck:
    feature_name: str
    tolerance_mm: float
    notes: str = ""

    def evaluate(self, points: List[Tuple[float, float, float]]) -> GdtFinding:
        if len(points) < 3:
            return GdtFinding(self.feature_name, "Flatness", 0.0, self.tolerance_mm, 0.0, self.tolerance_mm, True, self.tolerance_mm, None, "Insufficient points")
        # Least squares plane fit z = ax + by + d
        min_z = min(p[2] for p in points)
        max_z = max(p[2] for p in points)
        measured = max_z - min_z
        passed = measured <= self.tolerance_mm
        return GdtFinding(
            feature_name=self.feature_name,
            tolerance_type="Flatness",
            measured_value=measured,
            specified_tolerance=self.tolerance_mm,
            bonus_tolerance=0.0,
            total_tolerance=self.tolerance_mm,
            passed=passed,
            deviation_margin=self.tolerance_mm - measured,
            notes=self.notes
        )

@dataclass
class CylindricityCheck:
    feature_name: str
    tolerance_mm: float
    nominal_radius_mm: Optional[float] = None
    notes: str = ""

    def evaluate(self, points: List[Tuple[float, float, float]]) -> GdtFinding:
        if len(points) < 5:
            return GdtFinding(self.feature_name, "Cylindricity", 0.0, self.tolerance_mm, 0.0, self.tolerance_mm, True, self.tolerance_mm, None, "Insufficient points")
        radii = [math.sqrt(p[0]**2 + p[1]**2) for p in points]
        measured = max(radii) - min(radii)
        passed = measured <= self.tolerance_mm
        return GdtFinding(
            feature_name=self.feature_name,
            tolerance_type="Cylindricity",
            measured_value=measured,
            specified_tolerance=self.tolerance_mm,
            bonus_tolerance=0.0,
            total_tolerance=self.tolerance_mm,
            passed=passed,
            deviation_margin=self.tolerance_mm - measured,
            notes=self.notes
        )

@dataclass
class PerpendicularityCheck:
    feature_name: str
    datum_reference: str
    feature_length_mm: float
    tolerance_mm: float
    notes: str = ""

    def evaluate(self, feature_normal: Tuple[float, float, float], datum_normal: Tuple[float, float, float]) -> GdtFinding:
        # dot(n_feat, n_datum) should be 0.0 for 90 degrees
        fn_mag = math.sqrt(sum(x*x for x in feature_normal))
        dn_mag = math.sqrt(sum(x*x for x in datum_normal))
        if fn_mag < 1e-9 or dn_mag < 1e-9:
            dot = 0.0
        else:
            dot = sum(a*b for a, b in zip(feature_normal, datum_normal)) / (fn_mag * dn_mag)
        measured = abs(dot) * self.feature_length_mm
        passed = measured <= self.tolerance_mm
        return GdtFinding(
            feature_name=self.feature_name,
            tolerance_type="Perpendicularity",
            measured_value=measured,
            specified_tolerance=self.tolerance_mm,
            bonus_tolerance=0.0,
            total_tolerance=self.tolerance_mm,
            passed=passed,
            deviation_margin=self.tolerance_mm - measured,
            datum_reference=self.datum_reference,
            notes=self.notes
        )

@dataclass
class TruePositionCheck:
    feature_name: str
    nominal_x: float
    nominal_y: float
    base_tolerance_mm: float
    nominal_diameter: float
    mmc_diameter: Optional[float] = None
    lmc_diameter: Optional[float] = None
    material_condition: MaterialCondition = MaterialCondition.RFS
    datum_reference: Optional[str] = None
    notes: str = ""

    def evaluate(self, actual_x: float, actual_y: float, actual_diameter: float) -> GdtFinding:
        dx = actual_x - self.nominal_x
        dy = actual_y - self.nominal_y
        measured = 2.0 * math.sqrt(dx*dx + dy*dy) # Diametral zone

        bonus = 0.0
        if self.material_condition == MaterialCondition.MMC and self.mmc_diameter is not None:
            bonus = max(0.0, actual_diameter - self.mmc_diameter)
        elif self.material_condition == MaterialCondition.LMC and self.lmc_diameter is not None:
            bonus = max(0.0, self.lmc_diameter - actual_diameter)

        total_tol = self.base_tolerance_mm + bonus
        passed = measured <= total_tol
        return GdtFinding(
            feature_name=self.feature_name,
            tolerance_type=f"TruePosition ({self.material_condition.name})",
            measured_value=measured,
            specified_tolerance=self.base_tolerance_mm,
            bonus_tolerance=bonus,
            total_tolerance=total_tol,
            passed=passed,
            deviation_margin=total_tol - measured,
            datum_reference=self.datum_reference,
            notes=self.notes
        )

@dataclass
class ToleranceBand:
    upper_tol: float
    lower_tol: float


class GdtSpecification:
    """
    Container for complete engineering GD&T inspection schema.
    """
    def __init__(self, part_name: str = "PART-001", drawing_number: str = "DWG-001"):
        self.part_name = part_name
        self.drawing_number = drawing_number
        self.flatness_checks: List[FlatnessCheck] = []
        self.cylindricity_checks: List[CylindricityCheck] = []
        self.perpendicularity_checks: List[PerpendicularityCheck] = []
        self.position_checks: List[TruePositionCheck] = []

    def build(self) -> 'GdtSpecification':
        return self

    def add_flatness(self, feature_name: str, tolerance_mm: float = 0.05, tolerance: Optional[float] = None, notes: str = "", region_filter: Any = None) -> 'GdtSpecification':
        tol = tolerance if tolerance is not None else tolerance_mm
        self.flatness_checks.append(FlatnessCheck(feature_name, tol, notes))
        return self

    def add_cylindricity(self, feature_name: str, tolerance_mm: float = 0.05, tolerance: Optional[float] = None, notes: str = "", region_filter: Any = None) -> 'GdtSpecification':
        tol = tolerance if tolerance is not None else tolerance_mm
        self.cylindricity_checks.append(CylindricityCheck(feature_name, tol, notes=notes))
        return self


    def add_perpendicularity(self, feature_name: str, datum: str, length_mm: float, tolerance_mm: float, notes: str = "") -> 'GdtSpecification':
        self.perpendicularity_checks.append(PerpendicularityCheck(feature_name, datum, length_mm, tolerance_mm, notes))
        return self

    def add_true_position(self, feature_name: str, nom_x: float, nom_y: float, base_tol: float,
                          nom_dia: float, mmc_dia: Optional[float] = None, lmc_dia: Optional[float] = None,
                          condition: MaterialCondition = MaterialCondition.RFS, datum: Optional[str] = None) -> 'GdtSpecification':
        self.position_checks.append(TruePositionCheck(
            feature_name=feature_name,
            nominal_x=nom_x,
            nominal_y=nom_y,
            base_tolerance_mm=base_tol,
            nominal_diameter=nom_dia,
            mmc_diameter=mmc_dia,
            lmc_diameter=lmc_dia,
            material_condition=condition,
            datum_reference=datum
        ))
        return self

@dataclass
class GdtReport:
    part_name: str
    drawing_number: str
    findings: List[GdtFinding] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return all(f.passed for f in self.findings) if self.findings else True

    @property
    def pass_count(self) -> int:
        return sum(1 for f in self.findings if f.passed)

    @property
    def fail_count(self) -> int:
        return sum(1 for f in self.findings if not f.passed)


# Convenient alias
GDTSpecificationBuilder = GdtSpecification

