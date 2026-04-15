"""
caliper/scenarios/turbine_blade.py
Aerospace Inconel 718 High-Pressure Gas Turbine Blade Inspection Scenario.

Simulates automated dimensional inspection of a twisted, tapered transonic
turbine blade. Validates airfoil surface profile, root platform flatness,
and isolates localized leading-edge casting flash/burrs using ASME Y14.5 rules.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

from caliper.gdt_rules import FlatnessCheck, GDTSpecificationBuilder, ToleranceBand
from caliper.pipeline import InspectionConfig, InspectionPipeline, InspectionResult
from caliper.report import (
    export_point_deviations_csv,
    generate_json_manifest,
    generate_markdown_certificate,
    save_svg_cross_section,
)
from caliper.synthetic import (
    generate_turbine_blade_mesh,
    inject_sensor_noise_and_defects,
    sample_mesh_surface,
)
from caliper.viewer import save_html_viewer


def run_turbine_blade_scenario(
    output_dir: Optional[str | Path] = None,
    num_points: int = 15000,
    inject_defects: bool = True,
    noise_sigma: float = 0.012,
) -> InspectionResult:
    """Execute complete turbine blade inspection workflow.

    Args:
        output_dir: Directory to save generated certificates, SVGs, and 3D viewer.
        num_points: Number of surface points to scan.
        inject_defects: Whether to inject manufacturing flash defect.
        noise_sigma: Optical sensor Gaussian noise standard deviation in mm.

    Returns:
        The evaluated InspectionResult object.
    """
    out_path = Path(output_dir) if output_dir else Path("artifacts/turbine_blade")
    out_path.mkdir(parents=True, exist_ok=True)

    cad_mesh = generate_turbine_blade_mesh(
        span=100.0,
        root_chord=60.0,
        tip_chord=35.0,
        twist_deg=25.0,
        num_span_stations=20,
        points_per_station=36,
    )

    clean_points = sample_mesh_surface(cad_mesh, target_count=num_points, seed=101)

    defect_specs = []
    if inject_defects:
        defect_specs.append({
            "type": "burr",
            "center": (15.0, 5.0, 50.0),
            "radius": 8.0,
            "amplitude": 0.38,
        })

    scan_cloud = inject_sensor_noise_and_defects(
        clean_points,
        noise_sigma=noise_sigma,
        defects=defect_specs,
        seed=202,
    )

    gdt_builder = GDTSpecificationBuilder()
    gdt_builder.add_flatness(
        feature_name="RootPlatform",
        tolerance=0.08,
        region_filter=lambda p: p[2] < 3.0,
    )
    gdt_builder.add_flatness(
        feature_name="TipShroud",
        tolerance=0.10,
        region_filter=lambda p: p[2] > 97.0,
    )

    config = InspectionConfig(
        tolerance_band=ToleranceBand(upper_tol=0.25, lower_tol=-0.25),
        gdt_specification=gdt_builder.build(),
        sor_filter_active=True,
        sor_k_neighbors=16,
        sor_std_ratio=2.2,
        estimate_normals=True,
        defect_threshold_mm=0.25,
        defect_min_cluster_points=10,
        cluster_eps_mm=4.0,
        max_allowable_rmse=0.18,
    )

    pipeline = InspectionPipeline(config=config)
    result = pipeline.run(
        nominal_mesh=cad_mesh,
        scanned_points=scan_cloud,
        metadata={
            "part_id": "INCONEL-718-TB-402",
            "part_description": "HP Turbine Stage 1 Airfoil Rotor Blade",
            "material": "Inconel 718 Nickel Superalloy",
            "operator": "Coordinate Metrology Automation",
            "drawing_ref": "DWG-TB-718-A3",
        },
    )

    save_html_viewer(
        result,
        output_path=out_path / "dashboard.html",
        title="Inconel 718 Turbine Blade Metrology Dashboard",
        tolerance_band=0.30,
    )
    save_svg_cross_section(
        result,
        output_path=out_path / "cross_section_z50.svg",
        plane_axis="z",
        plane_coord=50.0,
        slab_thickness=2.5,
        tolerance_band=0.25,
    )
    cert_md = generate_markdown_certificate(result)
    (out_path / "certificate.md").write_text(cert_md, encoding="utf-8")

    cert_json = generate_json_manifest(result)
    (out_path / "manifest.json").write_text(cert_json, encoding="utf-8")

    export_point_deviations_csv(result, out_path / "points_log.csv")

    cad_mesh.to_stl(out_path / "nominal_cad.stl")
    scan_cloud.to_ply(out_path / "scanned_blade.ply")

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Turbine Blade Inspection Metrology Runner")
    parser.add_argument("--output-dir", type=str, default="artifacts/turbine_blade", help="Output artifact directory")
    parser.add_argument("--points", type=int, default=12000, help="Surface sampling point count")
    parser.add_argument("--no-defects", action="store_true", help="Omit synthetic defect injection")
    args = parser.parse_args()

    res = run_turbine_blade_scenario(
        output_dir=args.output_dir,
        num_points=args.points,
        inject_defects=not args.no_defects,
    )
    print(f"Inspection Complete. Disposition: {res.disposition}")
    print(f"Points: {res.statistics.point_count:,} | RMSE: {res.statistics.rmse:.4f} mm | Defects: {len(res.defects)}")
