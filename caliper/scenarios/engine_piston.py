"""
caliper/scenarios/engine_piston.py
Automotive Hypereutectic Aluminum Piston Inspection Scenario.

Simulates dimensional verification of an automotive internal combustion engine
piston. Evaluates combustion crown top flatness, cylindrical skirt form
tolerances (cylindricity), and detects localized casting shrinkage porosity voids.
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
    generate_engine_piston_mesh,
    inject_sensor_noise_and_defects,
    sample_mesh_surface,
)
from caliper.viewer import save_html_viewer


def run_engine_piston_scenario(
    output_dir: Optional[str | Path] = None,
    num_points: int = 15000,
    inject_defects: bool = True,
    noise_sigma: float = 0.010,
) -> InspectionResult:
    """Execute complete engine piston inspection workflow.

    Args:
        output_dir: Directory to save generated certificates, SVGs, and 3D viewer.
        num_points: Number of surface points to scan.
        inject_defects: Whether to inject manufacturing shrinkage void.
        noise_sigma: Optical sensor Gaussian noise standard deviation in mm.

    Returns:
        The evaluated InspectionResult object.
    """
    out_path = Path(output_dir) if output_dir else Path("artifacts/engine_piston")
    out_path.mkdir(parents=True, exist_ok=True)

    cad_mesh = generate_engine_piston_mesh(
        radius=40.0,
        height=65.0,
        pin_bore_radius=11.0,
        pin_bore_height=28.0,
        radial_slices=48,
        height_slices=24,
    )

    clean_points = sample_mesh_surface(cad_mesh, target_count=num_points, seed=303)

    defect_specs = []
    if inject_defects:
        defect_specs.append({
            "type": "void",
            "center": (15.0, 10.0, 65.0),
            "radius": 6.0,
            "amplitude": -0.32,
        })
        defect_specs.append({
            "type": "chatter",
            "center": (0.0, 0.0, 30.0),
            "radius": 40.0,
            "amplitude": 0.04,
            "frequency": 8.0,
        })

    scan_cloud = inject_sensor_noise_and_defects(
        clean_points,
        noise_sigma=noise_sigma,
        defects=defect_specs,
        seed=404,
    )

    gdt_builder = GDTSpecificationBuilder()
    gdt_builder.add_flatness(
        feature_name="CombustionCrown",
        tolerance=0.035,
        region_filter=lambda p: p[2] > 64.0,
    )
    gdt_builder.add_cylindricity(
        feature_name="PistonSkirt",
        tolerance=0.045,
        region_filter=lambda p: p[2] < 55.0,
    )

    config = InspectionConfig(
        tolerance_band=ToleranceBand(upper_tol=0.050, lower_tol=-0.050),
        gdt_specification=gdt_builder.build(),
        sor_filter_active=True,
        sor_k_neighbors=16,
        sor_std_ratio=2.0,
        estimate_normals=True,
        defect_threshold_mm=0.060,
        defect_min_cluster_points=10,
        cluster_eps_mm=4.0,
        max_allowable_rmse=0.045,
    )

    pipeline = InspectionPipeline(config=config)
    result = pipeline.run(
        nominal_mesh=cad_mesh,
        scanned_points=scan_cloud,
        metadata={
            "part_id": "PISTON-AL-4032-V6",
            "part_description": "Hypereutectic 4032 Alloy 80mm Piston",
            "material": "Al-Si-Cu-Mg Hypereutectic Alloy",
            "operator": "Inline Production QC Cell #4",
            "drawing_ref": "DWG-ENG-P80-REV2",
        },
    )

    save_html_viewer(
        result,
        output_path=out_path / "dashboard.html",
        title="Engine Piston Metrology Dashboard",
        tolerance_band=0.10,
    )
    save_svg_cross_section(
        result,
        output_path=out_path / "cross_section_z65.svg",
        plane_axis="z",
        plane_coord=65.0,
        slab_thickness=2.0,
        tolerance_band=0.050,
    )
    cert_md = generate_markdown_certificate(result)
    (out_path / "certificate.md").write_text(cert_md, encoding="utf-8")

    cert_json = generate_json_manifest(result)
    (out_path / "manifest.json").write_text(cert_json, encoding="utf-8")

    export_point_deviations_csv(result, out_path / "points_log.csv")

    cad_mesh.to_stl(out_path / "nominal_cad.stl")
    scan_cloud.to_ply(out_path / "scanned_piston.ply")

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Engine Piston Inspection Scenario Runner")
    parser.add_argument("--output-dir", type=str, default="artifacts/engine_piston", help="Output artifact directory")
    parser.add_argument("--points", type=int, default=12000, help="Surface sampling point count")
    parser.add_argument("--no-defects", action="store_true", help="Omit synthetic defect injection")
    args = parser.parse_args()

    res = run_engine_piston_scenario(
        output_dir=args.output_dir,
        num_points=args.points,
        inject_defects=not args.no_defects,
    )
    print(f"Inspection Complete. Disposition: {res.disposition}")
    print(f"Points: {res.statistics.point_count:,} | RMSE: {res.statistics.rmse:.4f} mm | Defects: {len(res.defects)}")
