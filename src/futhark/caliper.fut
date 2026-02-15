-- ============================================================================
-- Caliper Metrology Engine: Master Orchestration Module
-- Module: caliper.fut
--
-- Top-level data-parallel entry points orchestrating end-to-end optical
-- inspection, stereo 3D reconstruction, point-to-plane ICP alignment,
-- signed distance field deviation mapping, and ASME Y14.5 GD&T evaluation.
-- ============================================================================

import "math"
import "calib"
import "census"
import "sgm"
import "filter"
import "normals"
import "mesh"
import "bvh"
import "distance"
import "icp"
import "fit"
import "gdt"
import "defects"

-- ----------------------------------------------------------------------------
-- Top-Level Pipeline: Stereo Images -> Dense Calibrated 3D Point Cloud
-- ----------------------------------------------------------------------------

entry caliper_stereo_pipeline [h] [w]
    (fx: f64, fy: f64, cx: f64, cy: f64, baseline_mm: f64)
    (max_disp: i64, p1: i32, p2: i32, min_disp: f64)
    (left_img: [h][w]f64)
    (right_img: [h][w]f64): ([h][w]f64, [h][w]f64, [h][w]f64, [h][w]bool) =
  let c_left = compute_census_map left_img
  let c_right = compute_census_map right_img
  let raw_cost = compute_cost_volume c_left c_right max_disp
  let agg_cost = aggregate_sgm_4path p1 p2 raw_cost
  let disp_map = compute_refined_disparity agg_cost
  let cam: camera_intrinsics = { fx = fx, fy = fy, cx = cx, cy = cy, skew = 0.0 }
  let pts_3d = backproject_to_3d cam baseline_mm min_disp disp_map
  let xs = map (map (\(p, _) -> p.0)) pts_3d
  let ys = map (map (\(p, _) -> p.1)) pts_3d
  let zs = map (map (\(p, _) -> p.2)) pts_3d
  let valids = map (map (\(_, v) -> v)) pts_3d
  in (xs, ys, zs, valids)

-- ----------------------------------------------------------------------------
-- Top-Level Pipeline: Point Cloud Filtering & Normal Estimation
-- ----------------------------------------------------------------------------

entry caliper_filter_and_normals [n]
    (xs: [n]f64, ys: [n]f64, zs: [n]f64)
    (k_sor: i64, alpha_sor: f64)
    (k_normals: i64)
    (vp_x: f64, vp_y: f64, vp_z: f64): ([n]bool, [n]f64, [n]f64, [n]f64, [n]f64) =
  let pts = tabulate n (\i -> (xs[i], ys[i], zs[i]))
  let inlier_flags = sor_filter pts k_sor alpha_sor
  let vp: vec3 = (vp_x, vp_y, vp_z)
  let (norms, curvs) = estimate_normals pts k_normals vp
  let nxs = map (\n_vec -> n_vec.0) norms
  let nys = map (\n_vec -> n_vec.1) norms
  let nzs = map (\n_vec -> n_vec.2) norms
  in (inlier_flags, nxs, nys, nzs, curvs)

-- ----------------------------------------------------------------------------
-- Top-Level Pipeline: Scanned Points vs Nominal CAD Mesh Inspection
-- ----------------------------------------------------------------------------

entry caliper_inspect_mesh [n] [m]
    (px: [n]f64, py: [n]f64, pz: [n]f64)
    (v0x: [m]f64, v0y: [m]f64, v0z: [m]f64)
    (v1x: [m]f64, v1y: [m]f64, v1z: [m]f64)
    (v2x: [m]f64, v2y: [m]f64, v2z: [m]f64)
    (tolerance_mm: f64): ([n]f64, [n]f64, [n]f64, [n]f64, [n]f64, f64, f64, f64, f64, f64, f64) =
  let pts = tabulate n (\i -> (px[i], py[i], pz[i]))
  let tris = tabulate m (\i ->
    ((v0x[i], v0y[i], v0z[i]), (v1x[i], v1y[i], v1z[i]), (v2x[i], v2y[i], v2z[i]))
  )
  let devs = compute_deviations pts tris
  let dists = map (\d -> d.distance) devs
  let s_dists = map (\d -> d.signed_distance) devs
  let colors = map (\d -> deviation_to_rgb d.signed_distance tolerance_mm) devs
  let red_channel = map (\(r, _, _) -> r) colors
  let green_channel = map (\(_, g, _) -> g) colors
  let blue_channel = map (\(_, _, b) -> b) colors
  let stats = compute_deviation_stats devs tolerance_mm
  in ( dists
     , s_dists
     , red_channel
     , green_channel
     , blue_channel
     , stats.mean_deviation
     , stats.std_dev
     , stats.rmse
     , stats.max_positive
     , stats.max_negative
     , stats.in_tolerance_pct
     )

-- ----------------------------------------------------------------------------
-- Top-Level Pipeline: Automated GD&T Verification Suite
-- ----------------------------------------------------------------------------

entry caliper_gdt_suite [n]
    (xs: [n]f64, ys: [n]f64, zs: [n]f64)
    (flatness_tol: f64)
    (cylindricity_tol: f64): (f64, bool, f64, bool) =
  let pts = tabulate n (\i -> (xs[i], ys[i], zs[i]))
  let f_res = check_flatness pts flatness_tol
  let c_res = check_cylindricity pts cylindricity_tol
  in (f_res.measured_value, f_res.passed, c_res.measured_value, c_res.passed)
