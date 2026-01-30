-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Primitive Fitting
-- Module: tests/futhark/test_fit.fut
--
-- Tests Orthogonal Distance Regression (ODR) plane fitting, sphere least-squares,
-- and cylinder axis extraction on synthetic geometric primitives.
-- ============================================================================

import "../../src/futhark/fit"
import "../../src/futhark/math"

def test_plane_fit: bool =
  -- Points on plane z = 0.5 * x - 0.2 * y + 3.0
  let n = 25i64
  let pts =
    tabulate n (\i ->
      let x = f64.i64 (i % 5) * 2.0
      let y = f64.i64 (i / 5) * 2.0
      let z = 0.5 * x - 0.2 * y + 3.0
      in (x, y, z)
    )
  let res = fit_plane pts
  -- Normal is (-0.5, 0.2, 1.0) / norm
  let target_n = v3_normalize (-0.5, 0.2, 1.0)
  let err_nx = f64.abs (res.normal.0 - target_n.0)
  let err_ny = f64.abs (res.normal.1 - target_n.1)
  let err_nz = f64.abs (res.normal.2 - target_n.2)
  let ok_flatness = res.flatness_error < 1e-6
  let ok_rmse = res.rmse < 1e-6
  in err_nx < 1e-6 && err_ny < 1e-6 && err_nz < 1e-6 && ok_flatness && ok_rmse

def test_sphere_fit: bool =
  -- Sphere with center (1.0, 2.0, 3.0) and radius 5.0
  let center: vec3 = (1.0, 2.0, 3.0)
  let radius = 5.0
  let n_phi = 6i64
  let n_theta = 8i64
  let pts =
    tabulate (n_phi * n_theta) (\i ->
      let iphi = i / n_theta
      let itheta = i % n_theta
      let phi = (f64.i64 (iphi + 1) / f64.i64 (n_phi + 1)) * f64.pi
      let theta = (f64.i64 itheta / f64.i64 n_theta) * 2.0 * f64.pi
      let p_sphere = (radius * f64.sin phi * f64.cos theta, radius * f64.sin phi * f64.sin theta, radius * f64.cos phi)
      in v3_add center p_sphere
    )
  let res = fit_sphere pts
  let err_c = v3_dist res.center center
  let err_r = f64.abs (res.radius - radius)
  in err_c < 1e-3 && err_r < 1e-3 && res.roundness_error < 1e-2

def test_cylinder_fit: bool =
  -- Cylinder aligned along Z axis with radius 3.0 centered at (2.0, -1.0, 0.0)
  let n_z = 6i64
  let n_th = 8i64
  let radius = 3.0
  let pts =
    tabulate (n_z * n_th) (\i ->
      let iz = i / n_th
      let ith = i % n_th
      let theta = (f64.i64 ith / f64.i64 n_th) * 2.0 * f64.pi
      let z = f64.i64 iz * 4.0
      in (2.0 + radius * f64.cos theta, -1.0 + radius * f64.sin theta, z)
    )
  let res = fit_cylinder pts
  let err_r = f64.abs (res.radius - radius)
  let ok_dir = f64.abs (f64.abs res.axis_dir.2 - 1.0) < 1e-3
  in err_r < 1e-3 && ok_dir && res.cylindricity_error < 1e-3

entry main: bool =
  test_plane_fit && test_sphere_fit && test_cylinder_fit
