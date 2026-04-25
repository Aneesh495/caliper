"""
caliper/scenarios/satellite_bracket.py
Aerospace Additive Manufacturing Ti-6Al-4V Satellite Bracket Inspection Scenario.

Simulates dimensional verification of an electron beam melted (EBM) additive
titanium satellite structural bracket. Analyzes residual stress thermal warpage
on the primary mounting interface pad and verifies perpendicularity to datum A.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

from caliper.gdt_rules import GDTSpecificationBuilder, ToleranceBand
from caliper.pipeline import InspectionConfig, InspectionPipeline, InspectionResult
from caliper.report import (
    export_point_deviations_csv,
    generate_json_manifest,
    generate_markdown_certificate,
    save_svg_cross_section,
)
from caliper.synthetic import (
    generate_satellite_bracket_mesh,
    inject_sensor_noise_and_defects,
    sample_mesh_surface,
)
from caliper.viewer import save_html_viewer


def run_satellite_bracket_scenario(
    output_dir: Optional[str | Path] = None,
    num_points: int = 15000,
    inject_defects: bool = True,
    noise_sigma: float = 0.015,
) -> InspectionResult:
    """Execute complete satellite bracket inspection workflow.

    Args:
        output_dir: Directory to save generated certificates, SVGs, and 3D viewer.
        num_points: Number of surface points to scan.
        inject_defects: Whether to inject additive thermal warpage defect.
        noise_sigma: Optical sensor Gaussian noise standard deviation in mm.

    Returns:
        The evaluated InspectionResult object.
    """
    out_path = Path(output_dir) if output_dir else Path("artifacts/satellite_bracket")
    out_path.mkdir(parents=True, exist_ok=True)

    cad_mesh = generate_satellite_bracket_mesh(
        base_length=100.0,
        base_width=60.0,
        base_thickness=12.0,
        flange_height=80.0,
        flange_thickness=12.0,
    )

    clean_points = sample_mesh_surface(cad_mesh, target_count=num_points, seed=505)

    defect_specs = []
    if inject_defects:
        defect_specs.append({
            "type": "warpage",
            "center": (50.0, 30.0, 0.0),
            "radius": 40.0,
            "amplitude": 0.28,
        })

    scan_cloud = inject_sensor_noise_and_defects(
        clean_points,
        noise_sigma=noise_sigma,
        defects=defect_specs,
        seed=606,
    )

    gdt_builder = GDTSpecificationBuilder()
    gdt_builder.add_flatness(
        feature_name="BaseMountingPad",
        tolerance=0.15,
        region_filter=lambda p: p[2] < 2.0,
    )
    gdt_builder.add_perpendicularity(
        feature_name="PayloadInterfaceFlange",
        datum="DatumA",
        length_mm=80.0,
        tolerance_mm=0.18,
    )

    config = InspectionConfig(
        tolerance_band=ToleranceBand(upper_tol=0.20, lower_tol=-0.20),
        gdt_specification=gdt_builder.build(),
        sor_filter_active=True,
        sor_k_neighbors=16,
        sor_std_ratio=2.2,
        estimate_normals=True,
        defect_threshold_mm=0.20,
        defect_min_cluster_points=12,
        cluster_eps_mm=5.0,
        max_allowable_rmse=0.15,
    )

    pipeline = InspectionPipeline(config=config)
    result = pipeline.run(
        nominal_mesh=cad_mesh,
        scanned_points=scan_cloud,
        metadata={
            "part_id": "SAT-TI64-BRK-088",
            "part_description": "Additively Manufactured Satellite Bus Star-Tracker Bracket",
            "material": "Ti-6Al-4V Grade 5 Titanium",
            "operator": "Automated Non-Contact Optical CMM",
            "drawing_ref": "DWG-SAT-STR-088-REV4",
        },
    )

    save_html_viewer(
        result,
        output_path=out_path / "dashboard.html",
        title="Satellite Bracket Metrology Dashboard",
        tolerance_band=0.30,
    )
    save_svg_cross_section(
        result,
        output_path=out_path / "cross_section_y30.svg",
        plane_axis="y",
        plane_coord=30.0,
        slab_thickness=3.0,
        tolerance_band=0.20,
    )
    cert_md = generate_markdown_certificate(result)
    (out_path / "certificate.md").write_text(cert_md, encoding="utf-8")

    cert_json = generate_json_manifest(result)
    (out_path / "manifest.json").write_text(cert_json, encoding="utf-8")

    export_point_deviations_csv(result, out_path / "points_log.csv")

    cad_mesh.to_stl(out_path / "nominal_cad.stl")
    scan_cloud.to_ply(out_path / "scanned_bracket.ply")

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Satellite Bracket Inspection Scenario Runner")
    parser.add_argument("--output-dir", type=str, default="artifacts/satellite_bracket", help="Output artifact directory")
    parser.add_argument("--points", type=int, default=12000, help="Surface sampling point count")
    parser.add_argument("--no-defects", action="store_true", help="Omit synthetic defect injection")
    args = parser.parse_args()

    res = run_satellite_bracket_scenario(
        output_dir=args.output_dir,
        num_points=args.points,
        inject_defects=not args.no_defects,
    )
    print(f"Inspection Complete. Disposition: {res.disposition}")
    print(f"Points: {res.statistics.point_count:,} | RMSE: {res.statistics.rmse:.4f} mm | Defects: {len(res.defects)}")
