-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Camera Calibration & Rectification
-- Module: tests/futhark/test_calib.fut
--
-- Tests Brown-Conrady forward/inverse contraction mapping, bicubic sampling,
-- and epipolar ray projection round trips.
-- ============================================================================

import "../../src/futhark/calib"
import "../../src/futhark/math"

def test_distortion_roundtrip: bool =
  let dist: lens_distortion =
    { k1 = -0.15
    , k2 = 0.05
    , k3 = -0.01
    , p1 = 0.001
    , p2 = -0.002
    , s1 = 0.0
    , s2 = 0.0
    }
  -- Test a grid of normalized coordinates within [-0.5, 0.5]
  let test_points =
    [ (-0.4, -0.3)
    , (-0.2, 0.1)
    , (0.0, 0.0)
    , (0.3, -0.2)
    , (0.4, 0.4)
    ]
  let check_pt (x: f64, y: f64): bool =
    let (xd, yd) = distort_point dist x y
    let (xr, yr) = undistort_point dist xd yd
    let err = f64.sqrt ((xr - x)**2.0 + (yr - y)**2.0)
    in err < 1e-6
  in all check_pt test_points

def test_ray_pixel_roundtrip: bool =
  let cam: camera_intrinsics =
    { fx = 1200.0
    , fy = 1200.0
    , cx = 640.0
    , cy = 480.0
    , skew = 0.0
    }
  let dist: lens_distortion =
    { k1 = -0.1
    , k2 = 0.02
    , k3 = 0.0
    , p1 = 0.0
    , p2 = 0.0
    , s1 = 0.0
    , s2 = 0.0
    }
  let pixels =
    [ (100.0, 100.0)
    , (640.0, 480.0)
    , (1200.0, 900.0)
    ]
  let check_px (u: f64, v: f64): bool =
    let ray = pixel_to_ray cam dist u v
    let (ur, vr) = ray_to_pixel cam dist ray
    let err = f64.sqrt ((ur - u)**2.0 + (vr - v)**2.0)
    in err < 1e-5
  in all check_px pixels

def test_bicubic_interpolation: bool =
  -- Create a known smooth 2D quadratic field f(x, y) = x^2 + y^2
  let h = 8i64
  let w = 8i64
  let img = tabulate_2d h w (\y x -> (f64.i64 x)**2.0 + (f64.i64 y)**2.0)
  -- Sample at fractional coordinate (3.5, 4.5)
  let val_sampled = sample_bicubic img 3.5 4.5
  -- Analytical value for quadratic is exact under Catmull-Rom up to small spline bias
  let expected = 3.5**2.0 + 4.5**2.0 -- 12.25 + 20.25 = 32.5
  let err = f64.abs (val_sampled - expected)
  in err < 0.25

entry main: bool =
  test_distortion_roundtrip && test_ray_pixel_roundtrip && test_bicubic_interpolation
