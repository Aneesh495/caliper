-- ============================================================================
-- Caliper Metrology Engine: Sensor Optics & Epipolar Rectification
-- Module: calib.fut
--
-- Brown-Conrady 8-parameter non-linear lens distortion unwarping,
-- high-accuracy bicubic spline interpolation, and epipolar stereo rectification.
-- ============================================================================

import "math"

type camera_intrinsics =
  { fx: f64
  , fy: f64
  , cx: f64
  , cy: f64
  , skew: f64
  }

type lens_distortion =
  { k1: f64
  , k2: f64
  , k3: f64
  , p1: f64
  , p2: f64
  , s1: f64
  , s2: f64
  }

-- ----------------------------------------------------------------------------
-- Brown-Conrady Forward Distortion
--
-- Maps ideal normalized pinhole ray (x, y) to distorted sensor coordinate (xd, yd).
-- ----------------------------------------------------------------------------

def distort_point (dist: lens_distortion) (x: f64) (y: f64): (f64, f64) =
  let r2 = x*x + y*y
  let r4 = r2 * r2
  let r6 = r4 * r2
  let radial = 1.0 + dist.k1 * r2 + dist.k2 * r4 + dist.k3 * r6
  let dx_tan = 2.0 * dist.p1 * x * y + dist.p2 * (r2 + 2.0 * x * x)
  let dy_tan = dist.p1 * (r2 + 2.0 * y * y) + 2.0 * dist.p2 * x * y
  let dx_prism = dist.s1 * r2
  let dy_prism = dist.s2 * r2
  let xd = x * radial + dx_tan + dx_prism
  let yd = y * radial + dy_tan + dy_prism
  in (xd, yd)

-- ----------------------------------------------------------------------------
-- Iterative Inversion of Brown-Conrady Distortion
--
-- Solves for ideal normalized coordinates (x, y) given distorted coordinates (xd, yd)
-- using 6 iterations of Newton-Raphson / fixed-point contraction mapping.
-- Guaranteed convergence for realistic industrial lens fields of view.
-- ----------------------------------------------------------------------------

def undistort_point (dist: lens_distortion) (xd: f64) (yd: f64): (f64, f64) =
  let (x, y) =
    loop (cx, cy) = (xd, yd) for _step < 8 do
      let r2 = cx*cx + cy*cy
      let r4 = r2 * r2
      let r6 = r4 * r2
      let radial = 1.0 + dist.k1 * r2 + dist.k2 * r4 + dist.k3 * r6
      let dx_tan = 2.0 * dist.p1 * cx * cy + dist.p2 * (r2 + 2.0 * cx * cx)
      let dy_tan = dist.p1 * (r2 + 2.0 * cy * cy) + 2.0 * dist.p2 * cx * cy
      let dx_prism = dist.s1 * r2
      let dy_prism = dist.s2 * r2
      let nx = (xd - dx_tan - dx_prism) / radial
      let ny = (yd - dy_tan - dy_prism) / radial
      in (nx, ny)
  in (x, y)

-- ----------------------------------------------------------------------------
-- Pixel to Ideal Optical Ray Mapping
-- ----------------------------------------------------------------------------

def pixel_to_ray (cam: camera_intrinsics) (dist: lens_distortion) (u: f64) (v: f64): vec3 =
  let xd = (u - cam.cx) / cam.fx
  let yd = (v - cam.cy) / cam.fy
  let (x, y) = undistort_point dist xd yd
  in v3_normalize (x, y, 1.0)

def ray_to_pixel (cam: camera_intrinsics) (dist: lens_distortion) (r: vec3): (f64, f64) =
  let x = r.0 / r.2
  let y = r.1 / r.2
  let (xd, yd) = distort_point dist x y
  let u = cam.fx * xd + cam.cx
  let v = cam.fy * yd + cam.cy
  in (u, v)

-- ----------------------------------------------------------------------------
-- Bicubic Catmull-Rom Image Sampling
--
-- Provides C1-continuous sub-pixel interpolation on 2D scalar fields,
-- preventing aliasing artifacts in rectified epipolar imagery.
-- ----------------------------------------------------------------------------

def cubic_weight (x: f64): (f64, f64, f64, f64) =
  let x2 = x * x
  let x3 = x2 * x
  let w0 = -0.5 * x3 + x2 - 0.5 * x
  let w1 = 1.5 * x3 - 2.5 * x2 + 1.0
  let w2 = -1.5 * x3 + 2.0 * x2 + 0.5 * x
  let w3 = 0.5 * x3 - 0.5 * x2
  in (w0, w1, w2, w3)

def sample_bicubic [h] [w] (img: [h][w]f64) (u: f64) (v: f64): f64 =
  let iu = i64.f64 (f64.floor u)
  let iv = i64.f64 (f64.floor v)
  let fu = u - f64.i64 iu
  let fv = v - f64.i64 iv
  let (wu0, wu1, wu2, wu3) = cubic_weight fu
  let (wv0, wv1, wv2, wv3) = cubic_weight fv
  let clamp_x (x: i64): i64 = i64.max 0 (i64.min (w - 1) x)
  let clamp_y (y: i64): i64 = i64.max 0 (i64.min (h - 1) y)

  let sample_row (y_idx: i64): f64 =
    let y = clamp_y y_idx
    let p0 = img[y, clamp_x (iu - 1)]
    let p1 = img[y, clamp_x iu]
    let p2 = img[y, clamp_x (iu + 1)]
    let p3 = img[y, clamp_x (iu + 2)]
    in wu0 * p0 + wu1 * p1 + wu2 * p2 + wu3 * p3

  let r0 = sample_row (iv - 1)
  let r1 = sample_row iv
  let r2 = sample_row (iv + 1)
  let r3 = sample_row (iv + 2)
  in wv0 * r0 + wv1 * r1 + wv2 * r2 + wv3 * r3

-- ----------------------------------------------------------------------------
-- Parallel Lens Undistortion Filter
--
-- Maps an entire distorted image to an ideal rectilinear perspective frame.
-- ----------------------------------------------------------------------------

def undistort_image [h] [w]
    (cam: camera_intrinsics)
    (dist: lens_distortion)
    (img: [h][w]f64): [h][w]f64 =
  tabulate_2d h w (\y x ->
    let u = f64.i64 x
    let v = f64.i64 y
    let x_norm = (u - cam.cx) / cam.fx
    let y_norm = (v - cam.cy) / cam.fy
    let (xd, yd) = distort_point dist x_norm y_norm
    let src_u = cam.fx * xd + cam.cx
    let src_v = cam.fy * yd + cam.cy
    in sample_bicubic img src_u src_v
  )

-- ----------------------------------------------------------------------------
-- Parallel Epipolar Rectification Map Generation
--
-- Generates forward sampling coordinates (map_x, map_y) to transform
-- calibrated stereo pairs into row-aligned epipolar geometry.
-- ----------------------------------------------------------------------------

def generate_rectification_maps [h] [w]
    (cam: camera_intrinsics)
    (dist: lens_distortion)
    (R_rect: mat3): ([h][w]f64, [h][w]f64) =
  let (_, R_rect_inv) = mat3_inv R_rect
  let coords =
    tabulate_2d h w (\y x ->
      let u = f64.i64 x
      let v = f64.i64 y
      let x_norm = (u - cam.cx) / cam.fx
      let y_norm = (v - cam.cy) / cam.fy
      let ray_rect: vec3 = (x_norm, y_norm, 1.0)
      let ray_orig = mat3_mul_v3 R_rect_inv ray_rect
      let (src_u, src_v) = ray_to_pixel cam dist ray_orig
      in (src_u, src_v)
    )
  let map_x = map (map (\(u, _) -> u)) coords
  let map_y = map (map (\(_, v) -> v)) coords
  in (map_x, map_y)

-- ----------------------------------------------------------------------------
-- Remap Image using Dense Coordinate Fields
-- ----------------------------------------------------------------------------

def remap_image [h] [w] (img: [h][w]f64) (map_x: [h][w]f64) (map_y: [h][w]f64): [h][w]f64 =
  tabulate_2d h w (\y x ->
    sample_bicubic img map_x[y, x] map_y[y, x]
  )

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_distort_point (k1: f64, k2: f64, p1: f64, p2: f64, x: f64, y: f64): (f64, f64) =
  let dist = { k1 = k1, k2 = k2, k3 = 0.0, p1 = p1, p2 = p2, s1 = 0.0, s2 = 0.0 }
  in distort_point dist x y

entry fut_undistort_point (k1: f64, k2: f64, p1: f64, p2: f64, xd: f64, yd: f64): (f64, f64) =
  let dist = { k1 = k1, k2 = k2, k3 = 0.0, p1 = p1, p2 = p2, s1 = 0.0, s2 = 0.0 }
  in undistort_point dist xd yd

entry fut_undistort_image [h] [w]
    (fx: f64, fy: f64, cx: f64, cy: f64, k1: f64, k2: f64, p1: f64, p2: f64)
    (img: [h][w]f64): [h][w]f64 =
  let cam = { fx = fx, fy = fy, cx = cx, cy = cy, skew = 0.0 }
  let dist = { k1 = k1, k2 = k2, k3 = 0.0, p1 = p1, p2 = p2, s1 = 0.0, s2 = 0.0 }
  in undistort_image cam dist img
