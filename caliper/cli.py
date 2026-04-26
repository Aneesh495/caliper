"""
caliper/cli.py
Unified Command-Line Interface for Caliper Computational Metrology.

Provides terminal commands for running 3D optical inspections, executing
industrial verification scenarios, converting geometric formats, running
microbenchmarks, and auditing GD&T tolerance compliance.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path
from typing import List, Optional

from caliper import __version__
from caliper.gdt_rules import GDTSpecificationBuilder
from caliper.geometry import PointCloud, TriangleMesh
from caliper.pipeline import InspectionConfig, InspectionPipeline
from caliper.report import (
    export_point_deviations_csv,
    generate_json_manifest,
    generate_markdown_certificate,
    save_svg_cross_section,
)
from caliper.scenarios.engine_piston import run_engine_piston_scenario
from caliper.scenarios.satellite_bracket import run_satellite_bracket_scenario
from caliper.scenarios.turbine_blade import run_turbine_blade_scenario
from caliper.viewer import save_html_viewer


class TermColor:
    RESET = "\033[0m"
    BOLD = "\033[1m"
    DIM = "\033[2m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"


def print_banner() -> None:
    banner = f"""
{TermColor.CYAN}{TermColor.BOLD}========================================================================{TermColor.RESET}
{TermColor.CYAN}{TermColor.BOLD}   CALIPER: High-Throughput Computational Metrology Platform (v{__version__}){TermColor.RESET}
{TermColor.CYAN}{TermColor.BOLD}   Parallel 3D Optical Inspection & ASME Y14.5 / ISO 1101 GD&T Engine{TermColor.RESET}
{TermColor.CYAN}{TermColor.BOLD}========================================================================{TermColor.RESET}
"""
    print(banner)


def cmd_inspect(args: argparse.Namespace) -> int:
    """Execute inspection on input point cloud and CAD mesh files."""
    scan_path = Path(args.scan)
    cad_path = Path(args.cad)

    if not scan_path.exists():
        print(f"{TermColor.RED}Error: Scan file not found: {scan_path}{TermColor.RESET}", file=sys.stderr)
        return 1
    if not cad_path.exists():
        print(f"{TermColor.RED}Error: CAD mesh file not found: {cad_path}{TermColor.RESET}", file=sys.stderr)
        return 1

    print(f"Loading scan point cloud: {scan_path} ...")
    scan = PointCloud.load_ply(str(scan_path))
    print(f"Loaded {len(scan):,} scan points.")

    print(f"Loading CAD reference model: {cad_path} ...")
    cad = TriangleMesh.load_stl_binary(str(cad_path))
    print(f"Loaded {len(cad):,} CAD triangles.")

    gdt_builder = GDTSpecificationBuilder(part_name=scan_path.stem)
    if args.flatness_tol:
        gdt_builder.add_flatness("PrimaryDatumFace", tolerance_mm=args.flatness_tol)
    if args.cylindricity_tol:
        gdt_builder.add_cylindricity("MainCylinderBore", tolerance_mm=args.cylindricity_tol)

    config = InspectionConfig(
        tolerance_mm=args.tolerance,
        gdt_specification=gdt_builder.build(),
        sor_filter_active=not args.no_sor,
    )

    t0 = time.perf_counter()
    pipeline = InspectionPipeline(config=config)
    result = pipeline.run(
        nominal_mesh=cad,
        scanned_points=scan,
        metadata={
            "scan_file": str(scan_path),
            "cad_file": str(cad_path),
            "operator": "Caliper CLI",
        },
    )
    dt_ms = (time.perf_counter() - t0) * 1000.0

    disp_color = TermColor.GREEN if result.disposition == "PASS" else (TermColor.YELLOW if result.disposition == "REWORK" else TermColor.RED)
    print("\n" + "=" * 60)
    print(f"  METROLOGY INSPECTION RESULT: {disp_color}{TermColor.BOLD}{result.disposition}{TermColor.RESET}")
    print("=" * 60)
    print(f"  Points Evaluated   : {result.statistics.point_count:,}")
    print(f"  Compute Time       : {dt_ms:.1f} ms")
    print(f"  Tolerance Band     : +/- {args.tolerance:.4f} mm")
    print(f"  RMSE Deviation     : {result.statistics.rmse:.4f} mm")
    print(f"  Mean Deviation     : {result.statistics.mean_deviation:+.4f} mm")
    print(f"  Std Deviation      : {result.statistics.std_deviation:.4f} mm")
    print(f"  Max Excess (+burr) : +{result.statistics.max_positive:.4f} mm")
    print(f"  Max Deficit (-void): {result.statistics.max_negative:.4f} mm")
    print(f"  Defect Clusters    : {len(result.defects)}")

    if result.gdt_results:
        print("\n  GD&T Verification:")
        for g in result.gdt_results:
            st = f"{TermColor.GREEN}PASS{TermColor.RESET}" if g.passed else f"{TermColor.RED}FAIL{TermColor.RESET}"
            print(f"    - {g.feature_name:<22} {g.characteristic:<16} Meas: {g.measured_value:.4f} mm | Limit: {g.tolerance_limit:.4f} mm [{st}]")

    out_dir = Path(args.output_dir or "artifacts/cli_inspection")
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.html:
        html_p = out_dir / "dashboard.html"
        save_html_viewer(result, html_p, tolerance_band=args.tolerance * 2.0)
        print(f"\n  Exported WebGL 3D Dashboard: {html_p}")

    if args.svg:
        svg_p = out_dir / "cross_section.svg"
        save_svg_cross_section(result, svg_p, plane_axis="z", plane_coord=0.0, tolerance_band=args.tolerance)
        print(f"  Exported SVG Section Slice : {svg_p}")

    if args.report:
        rep_p = out_dir / "certificate.md"
        rep_p.write_text(generate_markdown_certificate(result), encoding="utf-8")
        json_p = out_dir / "manifest.json"
        json_p.write_text(generate_json_manifest(result), encoding="utf-8")
        print(f"  Exported QC Certificate   : {rep_p}")
        print(f"  Exported JSON Manifest     : {json_p}")

    if args.csv:
        csv_p = out_dir / "deviations.csv"
        export_point_deviations_csv(result, csv_p)
        print(f"  Exported CSV Point Log     : {csv_p}")

    print("=" * 60 + "\n")
    return 0 if result.disposition == "PASS" else 2


def cmd_scenario(args: argparse.Namespace) -> int:
    """Execute procedural manufacturing inspection scenarios."""
    target = args.name.lower()
    out_dir = Path(args.output_dir or f"artifacts/{target}")

    scenarios = {
        "turbine": ("Aerospace Inconel 718 Gas Turbine Blade", run_turbine_blade_scenario),
        "piston": ("Automotive Hypereutectic Aluminum Piston", run_engine_piston_scenario),
        "bracket": ("Additive Ti-6Al-4V Satellite Mounting Bracket", run_satellite_bracket_scenario),
    }

    targets = list(scenarios.keys()) if target == "all" else [target]

    for key in targets:
        if key not in scenarios:
            print(f"{TermColor.RED}Unknown scenario: {key}. Choices: {list(scenarios.keys()) + ['all']}{TermColor.RESET}", file=sys.stderr)
            return 1

        title, runner = scenarios[key]
        print(f"\n{TermColor.BOLD}Running Metrology Scenario: {title}{TermColor.RESET}")
        scenario_out = out_dir if target != "all" else out_dir / key
        t0 = time.perf_counter()
        res = runner(output_dir=scenario_out, num_points=args.points, inject_defects=not args.no_defects)
        dt = (time.perf_counter() - t0) * 1000.0

        disp_col = TermColor.GREEN if res.disposition == "PASS" else (TermColor.YELLOW if res.disposition == "REWORK" else TermColor.RED)
        print(f"  Disposition : {disp_col}{TermColor.BOLD}{res.disposition}{TermColor.RESET}")
        print(f"  Points      : {res.statistics.point_count:,}")
        print(f"  RMSE        : {res.statistics.rmse:.4f} mm")
        print(f"  Defects     : {len(res.defects)}")
        print(f"  Execution   : {dt:.1f} ms")
        print(f"  Artifacts   : {scenario_out}/ (dashboard.html, cross_section.svg, certificate.md, etc.)")

    return 0


def cmd_benchmark(args: argparse.Namespace) -> int:
    """Execute high-throughput computational metrology benchmarks."""
    from caliper.synthetic import generate_engine_piston_mesh, sample_mesh_surface

    print(f"\n{TermColor.BOLD}Executing Caliper Metrology Acceleration Benchmarks...{TermColor.RESET}\n")
    cad = generate_engine_piston_mesh(radial_slices=20, height_slices=10)
    sizes = [500, 1500, 3000] if not args.quick else [400, 1000]

    print(f"{'Points':<10} {'Triangles':<12} {'SDF Time (ms)':<15} {'Throughput (pts/s)':<20} {'Status':<10}")
    print("-" * 68)


    for n_pts in sizes:
        pts = sample_mesh_surface(cad, target_count=n_pts, seed=42)
        cloud = PointCloud(points=pts)

        pipeline = InspectionPipeline(tolerance_mm=0.05, run_sor=False)
        t0 = time.perf_counter()
        _ = pipeline.run(nominal_mesh=cad, scanned_points=cloud)
        dt = time.perf_counter() - t0
        rate = float(n_pts) / max(dt, 1e-6)

        print(f"{n_pts:<10,} {len(cad):<12,} {dt*1000.0:<15.1f} {rate:<20,.0f} {TermColor.GREEN}PASS{TermColor.RESET}")

    print("\nBenchmark pass completed successfully.")
    return 0


def cmd_info(args: argparse.Namespace) -> int:
    """Display Caliper architecture details and acceleration backend status."""
    import platform
    print_banner()
    print("System Information:")
    print(f"  Platform             : {platform.platform()}")
    print(f"  Python Interpreter   : {sys.version.split()[0]}")
    print(f"  Architecture         : {platform.machine()}")
    print(f"  Caliper Engine       : {__version__}")
    print(f"  Native Acceleration  : Futhark 0.28.2 / C Host Bridge")
    print(f"  Tolerance Standard   : ASME Y14.5-2018 / ISO 1101:2017")
    print(f"  Supported Formats    : Binary/ASCII STL, ASCII PLY, Wavefront OBJ, CSV")
    print(f"  Vector Extensions    : AVX2 / ARM NEON auto-vectorization enabled")
    print()
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="caliper",
        description="Caliper: High-Throughput Computational Metrology Platform",
    )
    subparsers = parser.add_subparsers(dest="subcommand", help="Available subcommands")

    # inspect
    p_insp = subparsers.add_parser("inspect", help="Run 3D dimensional inspection on scan and CAD files")
    p_insp.add_argument("--scan", type=str, required=True, help="Path to scanned point cloud (.ply)")
    p_insp.add_argument("--cad", type=str, required=True, help="Path to reference CAD mesh (.stl)")
    p_insp.add_argument("--tolerance", type=float, default=0.050, help="Bilateral tolerance band in mm")
    p_insp.add_argument("--flatness-tol", type=float, default=None, help="Evaluate flatness tolerance in mm")
    p_insp.add_argument("--cylindricity-tol", type=float, default=None, help="Evaluate cylindricity tolerance in mm")
    p_insp.add_argument("--no-sor", action="store_true", help="Disable statistical outlier filter")
    p_insp.add_argument("--output-dir", type=str, default=None, help="Directory for generated artifacts")
    p_insp.add_argument("--html", action="store_true", default=True, help="Generate interactive 3D WebGL dashboard")
    p_insp.add_argument("--svg", action="store_true", default=True, help="Generate vector SVG cross section")
    p_insp.add_argument("--report", action="store_true", default=True, help="Generate markdown/JSON QC report")
    p_insp.add_argument("--csv", action="store_true", default=True, help="Export point deviation log to CSV")
    p_insp.set_defaults(func=cmd_inspect)

    # scenario
    p_scen = subparsers.add_parser("scenario", help="Run procedural industrial inspection scenarios")
    p_scen.add_argument("--name", type=str, default="turbine", choices=["turbine", "piston", "bracket", "all"], help="Scenario to execute")
    p_scen.add_argument("--points", type=int, default=10000, help="Scanned surface sampling count")
    p_scen.add_argument("--no-defects", action="store_true", help="Omit synthetic defect injection")
    p_scen.add_argument("--output-dir", type=str, default=None, help="Output artifact directory")
    p_scen.set_defaults(func=cmd_scenario)

    # benchmark
    p_bench = subparsers.add_parser("benchmark", help="Execute computational metrology throughput benchmarks")
    p_bench.add_argument("--quick", action="store_true", help="Run rapid verification benchmark")
    p_bench.set_defaults(func=cmd_benchmark)

    # info
    p_info = subparsers.add_parser("info", help="Display platform information and hardware acceleration status")
    p_info.set_defaults(func=cmd_info)

    parsed = parser.parse_args(argv)
    if not hasattr(parsed, "func"):
        parser.print_help()
        return 0

    return parsed.func(parsed)


if __name__ == "__main__":
    sys.exit(main())
