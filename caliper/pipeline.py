"""
Caliper End-to-End Inspection Pipeline Orchestrator.
Module: caliper.pipeline

Coordinates point cloud pre-processing, rigid registration, signed distance
field analysis, false-color scalar field generation, and ASME Y14.5 GD&T evaluation.
"""

import math
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Any
from caliper.geometry import PointCloud, TriangleMesh
from caliper.gdt_rules import GdtSpecification, GdtReport
from caliper.core import deviation_to_color_rgb

@dataclass
class InspectionConfig:
    tolerance_band: Any = None
    tolerance_mm: float = 0.050
    gdt_specification: Optional[GdtSpecification] = None
    sor_filter_active: bool = True
    sor_k_neighbors: int = 16
    sor_std_ratio: float = 2.0
    estimate_normals: bool = True
    defect_threshold_mm: float = 0.20
    defect_min_cluster_points: int = 10
    cluster_eps_mm: float = 5.0
    max_allowable_rmse: float = 0.15

    def __post_init__(self):
        if self.tolerance_band is not None and hasattr(self.tolerance_band, "upper_tol"):
            self.tolerance_mm = abs(self.tolerance_band.upper_tol)


@dataclass
class DefectInfo:
    cluster_id: int
    defect_type: str
    point_count: int
    peak_deviation: float
    bounding_box_min: Tuple[float, float, float]
    bounding_box_max: Tuple[float, float, float]


@dataclass
class StatisticsInfo:
    point_count: int
    rmse: float
    mean_deviation: float
    std_deviation: float
    max_positive: float
    max_negative: float


@dataclass
class InspectionResult:
    evaluated_points: int
    tolerance_mm: float
    mean_deviation_mm: float
    std_dev_mm: float
    rmse_mm: float
    max_positive_mm: float
    max_negative_mm: float
    in_tolerance_pct: float
    passed: bool
    scan_pointcloud: PointCloud
    cad_mesh: TriangleMesh
    gdt_report: Optional[GdtReport] = None
    defect_clusters: List[Dict[str, Any]] = field(default_factory=list)
    quality_disposition: str = "PASS"  # PASS, REWORK, SCRAP
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def scanned_points(self) -> PointCloud:
        return self.scan_pointcloud

    @property
    def nominal_mesh(self) -> TriangleMesh:
        return self.cad_mesh

    @property
    def signed_deviations(self) -> List[float]:
        return self.scan_pointcloud.deviations or []

    @property
    def disposition(self) -> str:
        return self.quality_disposition

    @property
    def statistics(self) -> StatisticsInfo:
        return StatisticsInfo(
            point_count=self.evaluated_points,
            rmse=self.rmse_mm,
            mean_deviation=self.mean_deviation_mm,
            std_deviation=self.std_dev_mm,
            max_positive=self.max_positive_mm,
            max_negative=self.max_negative_mm,
        )

    @property
    def defects(self) -> List[DefectInfo]:
        res = []
        for c in self.defect_clusters:
            centroid = c.get("centroid", (0.0, 0.0, 0.0))
            rad = 5.0
            res.append(DefectInfo(
                cluster_id=c.get("cluster_id", 1),
                defect_type=c.get("classification", "Defect"),
                point_count=c.get("point_count", 0),
                peak_deviation=c.get("peak_dev", 0.0),
                bounding_box_min=(centroid[0] - rad, centroid[1] - rad, centroid[2] - rad),
                bounding_box_max=(centroid[0] + rad, centroid[1] + rad, centroid[2] + rad),
            ))
        return res

    @property
    def gdt_results(self) -> List[Any]:
        if self.gdt_report is not None:
            return self.gdt_report.findings
        return []


class InspectionPipeline:
    def __init__(self, tolerance_mm: float = 0.050, run_sor: bool = True, config: Optional[InspectionConfig] = None):
        if config is not None:
            self.config = config
            self.tolerance_mm = config.tolerance_mm
            self.run_sor = config.sor_filter_active
        else:
            self.config = InspectionConfig(tolerance_mm=tolerance_mm, sor_filter_active=run_sor)
            self.tolerance_mm = tolerance_mm
            self.run_sor = run_sor

    def run(self, nominal_mesh: TriangleMesh, scanned_points: PointCloud, metadata: Optional[Dict[str, Any]] = None) -> InspectionResult:
        res = self.inspect(scanned_points, nominal_mesh, gdt_spec=self.config.gdt_specification)
        if metadata:
            res.metadata.update(metadata)
        return res

    def inspect(self,
                scan: PointCloud,
                cad: TriangleMesh,
                gdt_spec: Optional[GdtSpecification] = None) -> InspectionResult:

        """
        Executes end-to-end metrology inspection:
        1. Pre-process scan with SOR outlier filter
        2. Estimate surface normals
        3. Evaluate exact Signed Distance Field against CAD triangles
        4. Generate false-color RGB tolerance mapping
        5. Evaluate GD&T specification rules
        6. Cluster out-of-tolerance defect regions
        7. Compute quality disposition decision
        """
        processed_scan = scan
        if self.run_sor and len(scan) > 20:
            processed_scan = scan.filter_sor(k=10, alpha=1.5)

        if not processed_scan.has_normals:
            processed_scan.estimate_normals(k=8)

        # Compute signed distances to CAD mesh
        s_dists = cad.evaluate_signed_distances(processed_scan.points)
        processed_scan.deviations = s_dists

        # Map to false-color RGB
        colors = [deviation_to_color_rgb(d, self.tolerance_mm) for d in s_dists]
        processed_scan.colors = colors

        # Statistical analysis
        n = len(s_dists)
        if n > 0:
            mean_dev = sum(s_dists) / n
            rmse = math.sqrt(sum(d*d for d in s_dists) / n)
            var = sum((d - mean_dev)**2 for d in s_dists) / n
            std_dev = math.sqrt(var)
            max_pos = max(s_dists)
            max_neg = min(s_dists)
            in_tol_count = sum(1 for d in s_dists if abs(d) <= self.tolerance_mm)
            in_tol_pct = (in_tol_count / n) * 100.0
        else:
            mean_dev = std_dev = rmse = max_pos = max_neg = 0.0
            in_tol_pct = 100.0

        # Cluster defect regions
        defect_clusters = self._cluster_defects(processed_scan.points, s_dists, self.tolerance_mm)

        # Evaluate GD&T rules if specified
        gdt_report = None
        if gdt_spec:
            findings = []
            for fc in gdt_spec.flatness_checks:
                findings.append(fc.evaluate(processed_scan.points))
            for cc in gdt_spec.cylindricity_checks:
                findings.append(cc.evaluate(processed_scan.points))
            gdt_report = GdtReport(
                part_name=gdt_spec.part_name,
                drawing_number=gdt_spec.drawing_number,
                findings=findings
            )

        passed = (in_tol_pct >= 95.0 and
                  max_pos <= self.tolerance_mm * 1.5 and
                  max_neg >= -self.tolerance_mm * 1.5 and
                  (gdt_report.passed if gdt_report else True))

        # Disposition
        if passed:
            disposition = "PASS"
        elif all(c.get("peak_dev", 0.0) > 0.0 for c in defect_clusters) and max_pos < self.tolerance_mm * 3.0:
            disposition = "REWORK"
        else:
            disposition = "SCRAP"

        return InspectionResult(
            evaluated_points=n,
            tolerance_mm=self.tolerance_mm,
            mean_deviation_mm=mean_dev,
            std_dev_mm=std_dev,
            rmse_mm=rmse,
            max_positive_mm=max_pos,
            max_negative_mm=max_neg,
            in_tolerance_pct=in_tol_pct,
            passed=passed,
            scan_pointcloud=processed_scan,
            cad_mesh=cad,
            gdt_report=gdt_report,
            defect_clusters=defect_clusters,
            quality_disposition=disposition
        )

    def _cluster_defects(self,
                         points: List[Tuple[float, float, float]],
                         deviations: List[float],
                         tol: float) -> List[Dict[str, Any]]:
        """
        Groups out-of-tolerance points into spatial defect clusters.
        """
        defect_indices = [i for i, d in enumerate(deviations) if abs(d) > tol]
        if not defect_indices:
            return []

        # Simple spatial connected component clustering
        clusters = []
        visited = set()
        cluster_radius = 5.0 # mm

        for idx in defect_indices:
            if idx in visited:
                continue
            cluster_members = [idx]
            visited.add(idx)

            p = points[idx]
            for other_idx in defect_indices:
                if other_idx in visited:
                    continue
                q = points[other_idx]
                dist = math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2)
                if dist <= cluster_radius:
                    cluster_members.append(other_idx)
                    visited.add(other_idx)

            c_pts = [points[i] for i in cluster_members]
            c_devs = [deviations[i] for i in cluster_members]
            cx = sum(pt[0] for pt in c_pts) / len(c_pts)
            cy = sum(pt[1] for pt in c_pts) / len(c_pts)
            cz = sum(pt[2] for pt in c_pts) / len(c_pts)
            peak = max(c_devs, key=abs)

            # Morphology classification
            if peak > 0.0:
                classification = "Excess Material (Flash/Burr)"
            else:
                classification = "Missing Material (Sink Mark/Porosity)"

            clusters.append({
                "cluster_id": len(clusters) + 1,
                "point_count": len(cluster_members),
                "centroid": (round(cx, 3), round(cy, 3), round(cz, 3)),
                "peak_dev": round(peak, 4),
                "classification": classification
            })

        return clusters
