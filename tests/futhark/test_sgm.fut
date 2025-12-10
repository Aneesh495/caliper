-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Semi-Global Matching (SGM)
-- Module: tests/futhark/test_sgm.fut
--
-- Tests SGM dynamic programming transitions, smoothness penalty propagation,
-- sub-pixel parabolic refinement, and 3D back-projection geometry.
-- ============================================================================

import "../../src/futhark/sgm"
import "../../src/futhark/math"
import "../../src/futhark/calib"

def test_sgm_dp_transitions: bool =
  let p1 = 10i32
  let p2 = 40i32
  -- Previous path cost vector with minimum at disparity 4
  let prev_L: [8]i32 = [100, 80, 50, 20, 0, 20, 50, 80]
  -- Current raw matching costs with minimum at disparity 4
  let cost_p: [8]i32 = [60, 50, 40, 20, 5, 20, 40, 60]
  let next_L = sgm_dp_step p1 p2 cost_p prev_L
  -- Minimum should stay at disparity 4
  let (min_d, min_c) =
    loop (best_d, min_v) = (0, next_L[0]) for i in 1..<8 do
      if next_L[i] < min_v then (i, next_L[i]) else (best_d, min_v)
  in min_d == 4 && min_c == 5

def test_subpixel_refinement: bool =
  -- Symmetric V-curve costs centered at d = 5 -> sub-pixel delta should be exactly 0.0
  let costs_sym: [10]i32 = [50, 40, 30, 20, 10, 0, 10, 20, 30, 40]
  let d_sym = refine_subpixel costs_sym
  let ok_sym = f64.abs (d_sym - 5.0) < 1e-6

  -- Asymmetric costs: c_minus = 12, c0 = 2, c_plus = 6 -> minimum shifted towards right (delta > 0)
  let costs_asym: [10]i32 = [50, 40, 30, 20, 12, 2, 6, 20, 30, 40]
  let d_asym = refine_subpixel costs_asym
  let ok_asym = d_asym > 5.0 && d_asym < 5.5
  in ok_sym && ok_asym

def test_backprojection_geometry: bool =
  let cam: camera_intrinsics =
    { fx = 1000.0
    , fy = 1000.0
    , cx = 500.0
    , cy = 500.0
    , skew = 0.0
    }
  let baseline_mm = 100.0
  let min_disp = 1.0
  let h = 4i64
  let w = 4i64
  -- Disparity map with constant 10.0 pixels
  let disp_map = replicate h (replicate w 10.0)
  let pts = backproject_to_3d cam baseline_mm min_disp disp_map
  -- Expected Z = (1000 * 100) / 10 = 10000 mm
  let (p00, valid00) = pts[0, 0]
  let expected_z = 10000.0
  let expected_x = (0.0 - 500.0) * expected_z / 1000.0 -- -5000.0
  let expected_y = (0.0 - 500.0) * expected_z / 1000.0 -- -5000.0
  let ok_z = f64.abs (p00.2 - expected_z) < 1e-6
  let ok_x = f64.abs (p00.0 - expected_x) < 1e-6
  let ok_y = f64.abs (p00.1 - expected_y) < 1e-6
  in valid00 && ok_z && ok_x && ok_y

entry main: bool =
  test_sgm_dp_transitions && test_subpixel_refinement && test_backprojection_geometry
