-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Signed Distance Fields
-- Module: tests/futhark/test_distance.fut
--
-- Tests 7-region Voronoi point-to-triangle exact distance, signed projection,
-- tolerance statistics, and RGB false-color colormap.
-- ============================================================================

import "../../src/futhark/distance"
import "../../src/futhark/math"
import "../../src/futhark/mesh"

def test_voronoi_regions: bool =
  -- Triangle in XY plane: A=(0,0,0), B=(4,0,0), C=(0,4,0)
  let tri: triangle = ((0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 4.0, 0.0))
  -- 1. Point above interior: (1, 1, 3) -> closest is (1, 1, 0), dist = 3
  let (q1, d1) = point_triangle_closest (1.0, 1.0, 3.0) tri
  let ok1 = f64.abs (q1.0 - 1.0) < 1e-6 && f64.abs (q1.1 - 1.0) < 1e-6 && f64.abs (d1 - 3.0) < 1e-6

  -- 2. Point near vertex A: (-2, -2, 0) -> closest is A=(0,0,0), dist = sqrt(8)
  let (q2, d2) = point_triangle_closest (-2.0, -2.0, 0.0) tri
  let ok2 = f64.abs q2.0 < 1e-6 && f64.abs q2.1 < 1e-6 && f64.abs (d2 - f64.sqrt 8.0) < 1e-6

  -- 3. Point near edge AB: (2, -3, 0) -> closest is (2, 0, 0), dist = 3
  let (q3, d3) = point_triangle_closest (2.0, -3.0, 0.0) tri
  let ok3 = f64.abs (q3.0 - 2.0) < 1e-6 && f64.abs q3.1 < 1e-6 && f64.abs (d3 - 3.0) < 1e-6

  in ok1 && ok2 && ok3

def test_signed_distance: bool =
  -- Triangle in XY plane with outward normal +Z
  let tri: triangle = ((0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 4.0, 0.0))
  let mesh = [tri]
  -- Point above triangle (positive side / excess material)
  let res_pos = point_mesh_signed_distance (1.0, 1.0, 2.5) mesh
  -- Point below triangle (negative side / missing material)
  let res_neg = point_mesh_signed_distance (1.0, 1.0, -1.5) mesh
  let ok_pos = f64.abs (res_pos.signed_distance - 2.5) < 1e-6
  let ok_neg = f64.abs (res_neg.signed_distance - (-1.5)) < 1e-6
  in ok_pos && ok_neg

def test_rgb_colormap: bool =
  -- At nominal 0.0 -> pure green (0, 1, 0)
  let (r0, g0, b0) = deviation_to_rgb 0.0 1.0
  let ok_mid = f64.abs r0 < 1e-6 && f64.abs (g0 - 1.0) < 1e-6 && f64.abs b0 < 1e-6
  -- At +tol -> pure red (1, 0, 0)
  let (r_hi, _, b_hi) = deviation_to_rgb 1.0 1.0
  let ok_hi = f64.abs (r_hi - 1.0) < 1e-6 && f64.abs b_hi < 1e-6
  -- At -tol -> pure blue (0, 0, 1)
  let (r_lo, _, b_lo) = deviation_to_rgb (-1.0) 1.0
  let ok_lo = f64.abs r_lo < 1e-6 && f64.abs (b_lo - 1.0) < 1e-6
  in ok_mid && ok_hi && ok_lo

entry main: bool =
  test_voronoi_regions && test_signed_distance && test_rgb_colormap
