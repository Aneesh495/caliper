"""
Caliper: High-Throughput Computational Metrology & 3D Optical Inspection Platform.

Core modules:
  caliper.core       - Native Futhark/C acceleration engine and CFFI bridge
  caliper.geometry   - PointCloud, TriangleMesh, and spatial data structures
  caliper.gdt_rules  - Declarative ASME Y14.5 and ISO 1101 GD&T rule specifications
  caliper.pipeline   - End-to-end optical inspection orchestrator
  caliper.viewer     - Self-contained interactive 3D WebGL dashboard generator
  caliper.report     - Formal QC certificates, CSV logs, and SVG cross-sections
  caliper.synthetic  - Procedural CAD models with calibrated manufacturing defects
"""

__version__ = "1.4.0"
__author__ = "Aneesh Krishna Parthasarathy"
__license__ = "Apache-2.0"

from caliper.geometry import PointCloud, TriangleMesh
from caliper.gdt_rules import (
    GdtSpecification,
    FlatnessCheck,
    CylindricityCheck,
    PerpendicularityCheck,
    TruePositionCheck,
    MaterialCondition,
)
from caliper.pipeline import InspectionPipeline, InspectionResult
from caliper.viewer import generate_inspection_dashboard
from caliper.report import generate_qc_certificate, export_svg_cross_section

__all__ = [
    "PointCloud",
    "TriangleMesh",
    "GdtSpecification",
    "FlatnessCheck",
    "CylindricityCheck",
    "PerpendicularityCheck",
    "TruePositionCheck",
    "MaterialCondition",
    "InspectionPipeline",
    "InspectionResult",
    "generate_inspection_dashboard",
    "generate_qc_certificate",
    "export_svg_cross_section",
]
