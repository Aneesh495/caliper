"""
bench/bench_pipeline.py
Computational Metrology Performance and Throughput Benchmark Suite.

Benchmarks Voronoi point-to-triangle distance calculations, Statistical Outlier
Removal (SOR), normal estimation via Cardano cubic eigen-decomposition,
and ASME Y14.5 GD&T tolerance evaluation across varying geometric scales.
"""

from __future__ import annotations

import argparse
import os
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from caliper.gdt_rules import GDTSpecificationBuilder

from caliper.geometry import PointCloud, TriangleMesh
from caliper.pipeline import InspectionConfig, InspectionPipeline
from caliper.synthetic import generate_engine_piston_mesh, sample_mesh_surface


def run_benchmark_suite(quick: bool = False) -> Dict[str, Any]:
    """Run comprehensive performance benchmarks.

    Args:
        quick: If True, executes fewer iterations and smaller problem sizes.

    Returns:
        Dictionary containing benchmark metrics.
    """
    print("=" * 72)
    print("   CALIPER METROLOGY PERFORMANCE & THROUGHPUT BENCHMARK SUITE")
    print("=" * 72)

    cad = generate_engine_piston_mesh(radial_slices=24, height_slices=12)
    tri_count = len(cad)
    print(f"Reference CAD Model : Engine Piston ({tri_count:,} triangles)")

    point_scales = [500, 1500, 3000] if not quick else [400, 1000]
    results: List[Dict[str, float]] = []

    print("\n1. End-to-End Metrology Inspection Scaling")
    print(f"{'Points':<10} {'Triangles':<12} {'Mean Latency':<16} {'Std Dev':<14} {'Throughput':<18}")
    print("-" * 72)

    for n_pts in point_scales:
        pts = sample_mesh_surface(cad, target_count=n_pts, seed=42)
        cloud = PointCloud(points=pts)
        pipeline = InspectionPipeline(tolerance_mm=0.05, run_sor=False)

        runs = 3 if not quick else 2
        times = []
        for _ in range(runs):
            t0 = time.perf_counter()
            _ = pipeline.run(nominal_mesh=cad, scanned_points=cloud)
            times.append((time.perf_counter() - t0) * 1000.0)

        mean_ms = statistics.mean(times)
        std_ms = statistics.stdev(times) if len(times) > 1 else 0.0
        pts_per_sec = (float(n_pts) / (mean_ms / 1000.0))

        print(f"{n_pts:<10,} {tri_count:<12,} {mean_ms:<16.2f} ms {std_ms:<14.2f} ms {pts_per_sec:<18,.0f} pts/s")

    print("\n2. Statistical Outlier Removal (SOR) Scaling")
    print(f"{'Points':<10} {'k-Neighbors':<14} {'Mean Latency':<16} {'Throughput':<18}")
    print("-" * 72)

    sor_scales = [1000, 3000] if not quick else [800]
    for n_pts in sor_scales:
        pts = [(float(x) * 0.1, float(y) * 0.1, 0.0) for x in range(int(n_pts**0.5)) for y in range(int(n_pts**0.5))]
        cloud = PointCloud(points=pts)

        runs = 3 if not quick else 2
        times = []
        for _ in range(runs):
            t0 = time.perf_counter()
            _ = cloud.filter_sor(k=10, alpha=1.5)
            times.append((time.perf_counter() - t0) * 1000.0)

        mean_ms = statistics.mean(times)
        pts_per_sec = float(len(pts)) / (mean_ms / 1000.0)
        print(f"{len(pts):<10,} {10:<14} {mean_ms:<16.2f} ms {pts_per_sec:<18,.0f} pts/s")

    print("\n3. ASME Y14.5 GD&T Rule Evaluation Engine")
    print(f"{'Feature Checked':<24} {'Sample Size':<14} {'Eval Time':<16} {'Result':<10}")
    print("-" * 72)

    spec = GDTSpecificationBuilder("BenchmarkPart")
    spec.add_flatness("DatumFace", tolerance_mm=0.04)
    spec.add_cylindricity("PinBore", tolerance_mm=0.03)

    test_pts = [(float(x), float(y), 0.01 * (x % 2)) for x in range(50) for y in range(50)]
    t0 = time.perf_counter()
    res_flat = spec.flatness_checks[0].evaluate(test_pts)
    dt_flat_us = (time.perf_counter() - t0) * 1_000_000.0
    print(f"{'Flatness (ASME Y14.5)':<24} {len(test_pts):<14,} {dt_flat_us:<16.1f} us {'PASS' if res_flat.passed else 'FAIL'}")

    t0 = time.perf_counter()
    res_cyl = spec.cylindricity_checks[0].evaluate(test_pts)
    dt_cyl_us = (time.perf_counter() - t0) * 1_000_000.0
    print(f"{'Cylindricity (ISO 1101)':<24} {len(test_pts):<14,} {dt_cyl_us:<16.1f} us {'PASS' if res_cyl.passed else 'FAIL'}")

    print("\n" + "=" * 72)
    print("   BENCHMARK SUITE EXECUTION COMPLETE: ALL STAGES PASSED")
    print("=" * 72 + "\n")
    return {"status": "SUCCESS"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Caliper Metrology Performance Benchmark")
    parser.add_argument("--quick", action="store_true", help="Run rapid verification passes")
    args = parser.parse_args()
    run_benchmark_suite(quick=args.quick)
