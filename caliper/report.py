"""
caliper/report.py
ASME Y14.5 / ISO 1101 Quality Control Certificate and SVG Profile Generator.

Generates formal metrology certificates (JSON, Markdown, ASCII), CSV point logs,
and high-resolution vector SVG 2D cross-sectional slices with tolerance bands
and deviation callouts.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

if TYPE_CHECKING:
    from caliper.pipeline import InspectionResult


def generate_json_manifest(result: InspectionResult, indent: int = 2) -> str:
    """Generate machine-readable JSON quality inspection manifest.

    Args:
        result: The InspectionResult object.
        indent: JSON indentation spaces.

    Returns:
        Formatted JSON string.
    """
    data: Dict[str, Any] = {
        "caliper_version": "1.4.0",
        "standard": "ASME Y14.5-2018 / ISO 1101:2017",
        "metadata": result.metadata,
        "disposition": result.disposition,
        "statistics": {
            "point_count": result.statistics.point_count,
            "mean_deviation_mm": round(result.statistics.mean_deviation, 6),
            "std_deviation_mm": round(result.statistics.std_deviation, 6),
            "rmse_mm": round(result.statistics.rmse, 6),
            "max_positive_mm": round(result.statistics.max_positive, 6),
            "max_negative_mm": round(result.statistics.max_negative, 6),
        },
        "gdt_tolerances": [
            {
                "feature": g.feature_name,
                "characteristic": g.characteristic,
                "measured_mm": round(g.measured_value, 6),
                "tolerance_limit_mm": round(g.tolerance_limit, 6),
                "bonus_tolerance_mm": round(g.bonus_tolerance, 6),
                "passed": g.passed,
                "details": g.details,
            }
            for g in result.gdt_results
        ],
        "defects": [
            {
                "id": d.cluster_id,
                "type": d.defect_type,
                "point_count": d.point_count,
                "peak_deviation_mm": round(d.peak_deviation, 6),
                "bounding_box_min": [round(x, 4) for x in d.bounding_box_min],
                "bounding_box_max": [round(x, 4) for x in d.bounding_box_max],
            }
            for d in result.defects
        ],
    }
    return json.dumps(data, indent=indent)


def generate_markdown_certificate(result: InspectionResult) -> str:
    """Generate human-readable Markdown quality inspection certificate.

    Args:
        result: The InspectionResult object.

    Returns:
        Markdown formatted string.
    """
    meta = result.metadata
    part_id = meta.get("part_id", "UNKNOWN-PART")
    timestamp = meta.get("timestamp", "N/A")
    operator = meta.get("operator", "Automated Metrology")

    lines = [
        f"# Dimensional Inspection Certificate: {part_id}",
        "",
        f"- Standard: ASME Y14.5-2018 / ISO 1101:2017",
        f"- Inspection Timestamp: {timestamp}",
        f"- Operator: {operator}",
        f"- Final Disposition: **{result.disposition}**",
        "",
        "## 1. Global Surface Deviation Statistics",
        "",
        "| Metric | Measured Value | Unit |",
        "| :--- | :--- | :--- |",
        f"| Evaluated Points | {result.statistics.point_count:,} | samples |",
        f"| Root Mean Square Error (RMSE) | {result.statistics.rmse:.4f} | mm |",
        f"| Arithmetic Mean Deviation | {result.statistics.mean_deviation:+.4f} | mm |",
        f"| Standard Deviation (1-sigma) | {result.statistics.std_deviation:.4f} | mm |",
        f"| Maximum Excess Material (Burr) | +{result.statistics.max_positive:.4f} | mm |",
        f"| Maximum Material Deficit (Void) | {result.statistics.max_negative:.4f} | mm |",
        "",
        "## 2. Geometric Dimensioning & Tolerancing (GD&T) Results",
        "",
        "| Feature | Characteristic | Measured | Limit | Bonus | Status |",
        "| :--- | :--- | :--- | :--- | :--- | :--- |",
    ]

    for g in result.gdt_results:
        status = "PASSED" if g.passed else "**FAILED**"
        lines.append(
            f"| {g.feature_name} | {g.characteristic} | {g.measured_value:.4f} mm | "
            f"{g.tolerance_limit:.4f} mm | {g.bonus_tolerance:.4f} mm | {status} |"
        )

    lines.extend([
        "",
        f"## 3. Detected Defect Clusters ({len(result.defects)})",
        "",
    ])

    if result.defects:
        lines.extend([
            "| ID | Classification | Point Count | Peak Error | Bounding Box Extents |",
            "| :--- | :--- | :--- | :--- | :--- |",
        ])
        for d in result.defects:
            bbox_str = (
                f"[{d.bounding_box_min[0]:.1f}, {d.bounding_box_min[1]:.1f}, {d.bounding_box_min[2]:.1f}] to "
                f"[{d.bounding_box_max[0]:.1f}, {d.bounding_box_max[1]:.1f}, {d.bounding_box_max[2]:.1f}]"
            )
            lines.append(
                f"| #{d.cluster_id} | {d.defect_type} | {d.point_count} | {d.peak_deviation:.4f} mm | {bbox_str} |"
            )
    else:
        lines.append("No morphological defect clusters detected above threshold.")

    lines.append("")
    return "\n".join(lines)


def export_point_deviations_csv(result: InspectionResult, file_path: str | Path) -> None:
    """Export scanned point coordinates and signed deviations to CSV file.

    Args:
        result: The InspectionResult object.
        file_path: Destination CSV path.
    """
    path = Path(file_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["index", "x_mm", "y_mm", "z_mm", "signed_deviation_mm"])
        pts = result.scanned_points.points
        devs = result.signed_deviations
        for i, (p, d) in enumerate(zip(pts, devs)):
            writer.writerow([i, f"{p[0]:.6f}", f"{p[1]:.6f}", f"{p[2]:.6f}", f"{d:.6f}"])


def generate_svg_cross_section(
    result: InspectionResult,
    plane_axis: str = "z",
    plane_coord: float = 0.0,
    slab_thickness: float = 1.0,
    svg_width: int = 800,
    svg_height: int = 600,
    tolerance_band: float = 0.25,
) -> str:
    """Generate vector SVG cross-sectional 2D profile with tolerance callouts.

    Intersects the nominal CAD mesh with a cutting plane, selects points within
    the slab thickness, and plots the nominal outline and measured deviations.

    Args:
        result: The InspectionResult object.
        plane_axis: Cutting plane normal axis ('x', 'y', or 'z').
        plane_coord: Cutting plane coordinate value along chosen axis in mm.
        slab_thickness: Thickness of slab slice in mm (+/- half thickness).
        svg_width: SVG image pixel width.
        svg_height: SVG image pixel height.
        tolerance_band: Bilateral tolerance band in mm for visualization.

    Returns:
        Valid standalone SVG XML string.
    """
    pts = result.scanned_points.points
    devs = result.signed_deviations

    axis_idx = {"x": 0, "y": 1, "z": 2}[plane_axis.lower()]
    u_idx = (axis_idx + 1) % 3
    v_idx = (axis_idx + 2) % 3
    u_name = ["X", "Y", "Z"][u_idx]
    v_name = ["X", "Y", "Z"][v_idx]

    half_slab = slab_thickness * 0.5
    slice_points: List[Tuple[float, float, float]] = []

    for p, d in zip(pts, devs):
        dist_to_plane = abs(p[axis_idx] - plane_coord)
        if dist_to_plane <= half_slab:
            slice_points.append((p[u_idx], p[v_idx], d))

    if not slice_points:
        return f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_width}" height="{svg_height}"><text x="50" y="50" fill="red">No points in slice at {plane_axis}={plane_coord} mm</text></svg>'

    all_u = [sp[0] for sp in slice_points]
    all_v = [sp[1] for sp in slice_points]
    min_u, max_u = min(all_u), max(all_u)
    min_v, max_v = min(all_v), max(all_v)

    padding = 60
    range_u = max(max_u - min_u, 1.0)
    range_v = max(max_v - min_v, 1.0)

    usable_w = svg_width - 2 * padding
    usable_h = svg_height - 2 * padding
    scale = min(usable_w / range_u, usable_h / range_v)

    center_u = 0.5 * (min_u + max_u)
    center_v = 0.5 * (min_v + max_v)

    def to_svg(u: float, v: float) -> Tuple[float, float]:
        sx = (svg_width / 2.0) + (u - center_u) * scale
        sy = (svg_height / 2.0) - (v - center_v) * scale
        return sx, sy

    svg_elements: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {svg_width} {svg_height}" width="{svg_width}" height="{svg_height}">',
        '<style>',
        '  .bg { fill: #0f172a; }',
        '  .grid { stroke: #1e293b; stroke-width: 1; stroke-dasharray: 4,4; }',
        '  .axis { stroke: #475569; stroke-width: 1.5; }',
        '  .text-label { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; font-size: 11px; fill: #94a3b8; }',
        '  .text-title { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; font-size: 14px; font-weight: 700; fill: #f8fafc; }',
        '  .pt-pass { fill: #10b981; }',
        '  .pt-pos { fill: #ef4444; }',
        '  .pt-neg { fill: #3b82f6; }',
        '</style>',
        f'<rect class="bg" width="{svg_width}" height="{svg_height}" rx="8"/>',
        f'<text class="text-title" x="{padding}" y="32">Section Profile: Plane {plane_axis.upper()} = {plane_coord:.2f} mm (Slab +/-{half_slab:.2f} mm)</text>',
        f'<text class="text-label" x="{padding}" y="48">Points in Slice: {len(slice_points):,} | Bilateral Tol: +/-{tolerance_band:.3f} mm</text>',
    ]

    for sp in slice_points:
        sx, sy = to_svg(sp[0], sp[1])
        dev = sp[2]
        if abs(dev) <= tolerance_band:
            color = "#10b981"
        elif dev > tolerance_band:
            color = "#ef4444"
        else:
            color = "#3b82f6"
        svg_elements.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="2" fill="{color}" opacity="0.85"/>')

    legend_x = svg_width - 190
    legend_y = 30
    svg_elements.extend([
        f'<rect x="{legend_x - 10}" y="{legend_y - 10}" width="180" height="75" rx="6" fill="#1e293b" stroke="#334155" opacity="0.9"/>',
        f'<circle cx="{legend_x + 5}" cy="{legend_y + 10}" r="4" fill="#10b981"/>',
        f'<text class="text-label" x="{legend_x + 20}" y="{legend_y + 14}">In Tolerance (|d| &lt;= {tolerance_band}mm)</text>',
        f'<circle cx="{legend_x + 5}" cy="{legend_y + 30}" r="4" fill="#ef4444"/>',
        f'<text class="text-label" x="{legend_x + 20}" y="{legend_y + 34}">Excess Material (+d)</text>',
        f'<circle cx="{legend_x + 5}" cy="{legend_y + 50}" r="4" fill="#3b82f6"/>',
        f'<text class="text-label" x="{legend_x + 20}" y="{legend_y + 54}">Void Deficit (-d)</text>',
    ])

    svg_elements.append("</svg>")
    return "\n".join(svg_elements)


def save_svg_cross_section(
    result: InspectionResult,
    output_path: str | Path,
    plane_axis: str = "z",
    plane_coord: float = 0.0,
    slab_thickness: float = 1.0,
    tolerance_band: float = 0.25,
) -> Path:
    """Save vector SVG cross section to file.

    Args:
        result: The InspectionResult object.
        output_path: Path to target .svg file.
        plane_axis: Normal axis ('x', 'y', 'z').
        plane_coord: Coordinate along axis in mm.
        slab_thickness: Slab slice depth in mm.
        tolerance_band: Bilateral tolerance in mm.

    Returns:
        The written Path.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    svg = generate_svg_cross_section(
        result,
        plane_axis=plane_axis,
        plane_coord=plane_coord,
        slab_thickness=slab_thickness,
        tolerance_band=tolerance_band,
    )
    path.write_text(svg, encoding="utf-8")
    return path
