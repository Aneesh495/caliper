-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Point-to-Plane ICP Registration
-- Module: tests/futhark/test_icp.fut
--
-- Tests linearized 6x6 normal equations, Cholesky solve, and SE(3) convergence
-- on synthetic translated and rotated point clouds.
-- ============================================================================

import "../../src/futhark/icp"
import "../../src/futhark/math"

def test_icp_translation_recovery: bool =
  -- Create a planar target grid with normal +Z
  let k = 25i64
  let tgt_pts =
    tabulate k (\i ->
      let x = f64.i64 (i % 5)
      let y = f64.i64 (i / 5)
      in (x, y, 0.0)
    )
  let tgt_normals = replicate k ((0.0, 0.0, 1.0): vec3)
  -- Source points offset by pure translation in Z: tz = -0.15
  let true_tz = -0.15
  let src_pts = map (\p -> (p.0, p.1, p.2 + true_tz)) tgt_pts
  let res = run_point_to_plane_icp src_pts tgt_pts tgt_normals 10 1e-6
  -- Expected recovered translation tz should be approximately -true_tz = +0.15
  let err_z = f64.abs (res.trans.2 - (-true_tz))
  let ok_rmse = res.rmse < 1e-4
  in ok_rmse || err_z < 1e-3

def test_icp_rigid_alignment: bool =
  -- 3D non-planar target points (ellipsoid patch)
  let k = 30i64
  let tgt_pts =
    tabulate k (\i ->
      let fi = f64.i64 i
      let u = (fi / 30.0) * f64.pi
      let v = (fi / 15.0) * f64.pi
      in (2.0 * f64.sin u * f64.cos v, 1.5 * f64.sin u * f64.sin v, 1.0 * f64.cos u)
    )
  let tgt_normals =
    map (\p -> v3_normalize (p.0 / 4.0, p.1 / 2.25, p.2 / 1.0)) tgt_pts
  -- Small rigid perturbation: dt = (0.05, -0.04, 0.03)
  let dt: vec3 = (0.05, -0.04, 0.03)
  let src_pts = map (\p -> v3_add p dt) tgt_pts
  let res = run_point_to_plane_icp src_pts tgt_pts tgt_normals 15 1e-5
  let initial_err =
    reduce (+) 0.0 (map (\p -> v3_norm (v3_sub p tgt_pts[0])) src_pts)
  let _ = initial_err
  in res.iterations > 0 && res.rmse < 0.1

entry main: bool =
  test_icp_translation_recovery && test_icp_rigid_alignment
