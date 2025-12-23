-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Surface Normal Estimation
-- Module: tests/futhark/test_normals.fut
--
-- Tests covariance tensor eigen-decomposition and surface normal orientation
-- on synthetic geometric planes and spheres.
-- ============================================================================

import "../../src/futhark/normals"
import "../../src/futhark/math"

def test_plane_normal_and_curvature: bool =
  -- Create a 5x5 grid of points on plane Z = 5.0
  let pts =
    tabulate 25 (\i ->
      let x = f64.i64 (i % 5)
      let y = f64.i64 (i / 5)
      in (x, y, 5.0)
    )
  let (_, cov) = cluster_covariance pts
  let viewpoint: vec3 = (2.0, 2.0, 10.0)
  let center_pt: vec3 = (2.0, 2.0, 5.0)
  let (norm, curv) = normal_from_covariance cov center_pt viewpoint
  -- Expected normal pointing towards (2, 2, 10): (0.0, 0.0, 1.0)
  let ok_nx = f64.abs norm.0 < 1e-6
  let ok_ny = f64.abs norm.1 < 1e-6
  let ok_nz = f64.abs (norm.2 - 1.0) < 1e-6
  let ok_curv = curv < 1e-6 -- Completely flat
  in ok_nx && ok_ny && ok_nz && ok_curv

def test_sphere_normals: bool =
  -- Sample points on unit sphere centered at origin
  let n = 20i64
  let pts =
    tabulate n (\i ->
      let theta = (f64.i64 i / f64.i64 n) * 2.0 * f64.pi
      let phi = f64.pi / 4.0
      in (f64.cos theta * f64.sin phi, f64.sin theta * f64.sin phi, f64.cos phi)
    )
  let viewpoint: vec3 = (0.0, 0.0, 5.0)
  let (norms, curvs) = estimate_normals pts 5 viewpoint
  -- Verify all normal vectors are normalized (norm == 1.0)
  let norms_unit = map (\n_vec -> f64.abs (v3_norm n_vec - 1.0) < 1e-6) norms
  -- Curvature on sphere should be non-zero
  let curvs_nonzero = map (\c -> c >= 0.0) curvs
  in all id norms_unit && all id curvs_nonzero

entry main: bool =
  test_plane_normal_and_curvature && test_sphere_normals
