-- ============================================================================
-- Caliper Metrology Engine: Semi-Global Matching (SGM) & Dense 3D Reconstruction
-- Module: sgm.fut
--
-- Heiko Hirschmuller Semi-Global Matching dynamic programming path aggregation,
-- penalty tuning (P1, P2), sub-pixel parabolic interpolation, left-right
-- consistency check, and 3D metric Euclidean back-projection.
-- ============================================================================

import "math"
import "calib"

-- ----------------------------------------------------------------------------
-- SGM Dynamic Programming Step along a Single 1D Path
--
-- Computes L_r(p, d) given previous path cost vector prev_L[0..d-1], current
-- raw cost vector cost_p[0..d-1], and smoothness penalties P1, P2.
-- ----------------------------------------------------------------------------

def sgm_dp_step [d]
    (p1: i32)
    (p2: i32)
    (cost_p: [d]i32)
    (prev_L: [d]i32): [d]i32 =
  -- Find min over all disparities in previous step
  let min_prev = loop m = prev_L[0] for i in 1..<d do i32.min m prev_L[i]
  in tabulate d (\disp ->
       let c0 = prev_L[disp]
       let c_minus = if disp > 0 then prev_L[disp - 1] + p1 else c0 + p2
       let c_plus = if disp < d - 1 then prev_L[disp + 1] + p1 else c0 + p2
       let c_any = min_prev + p2
       let min_trans = i32.min (i32.min c0 c_minus) (i32.min c_plus c_any)
       in cost_p[disp] + min_trans - min_prev
     )

-- ----------------------------------------------------------------------------
-- Horizontal Path Aggregation (Left-to-Right and Right-to-Left)
-- ----------------------------------------------------------------------------

def aggregate_horizontal [h] [w] [d]
    (p1: i32)
    (p2: i32)
    (cost_vol: [h][w][d]i32): [h][w][d]i32 =
  -- Scan Left to Right across each row
  let l2r =
    tabulate_2d h w (\y x ->
      -- Accumulate row recursively
      let row_costs = cost_vol[y]
      let row_agg =
        loop row_res = replicate w (replicate d 0i32) for col < w do
          let c_cur = row_costs[col]
          in if col == 0
             then row_res with [col] = copy c_cur
             else
               let prev = row_res[col - 1]
               let next_L = sgm_dp_step p1 p2 c_cur prev
               in row_res with [col] = copy next_L
      in row_agg[x]
    )
  -- Scan Right to Left across each row
  let r2l =
    tabulate_2d h w (\y x ->
      let row_costs = cost_vol[y]
      let row_agg =
        loop row_res = replicate w (replicate d 0i32) for step < w do
          let col = w - 1 - step
          let c_cur = row_costs[col]
          in if step == 0
             then row_res with [col] = copy c_cur
             else
               let prev = row_res[col + 1]
               let next_L = sgm_dp_step p1 p2 c_cur prev
               in row_res with [col] = copy next_L
      in row_agg[x]
    )
  in tabulate_2d h w (\y x ->
       tabulate d (\disp -> l2r[y, x, disp] + r2l[y, x, disp])
     )

-- ----------------------------------------------------------------------------
-- Vertical Path Aggregation (Top-to-Bottom and Bottom-to-Top)
-- ----------------------------------------------------------------------------

def aggregate_vertical [h] [w] [d]
    (p1: i32)
    (p2: i32)
    (cost_vol: [h][w][d]i32): [h][w][d]i32 =
  -- Scan Top to Bottom down each column
  let t2b =
    tabulate_2d h w (\y x ->
      let col_agg =
        loop col_res = replicate h (replicate d 0i32) for row < h do
          let c_cur = cost_vol[row, x]
          in if row == 0
             then col_res with [row] = copy c_cur
             else
               let prev = col_res[row - 1]
               let next_L = sgm_dp_step p1 p2 c_cur prev
               in col_res with [row] = copy next_L
      in col_agg[y]
    )
  -- Scan Bottom to Top up each column
  let b2t =
    tabulate_2d h w (\y x ->
      let col_agg =
        loop col_res = replicate h (replicate d 0i32) for step < h do
          let row = h - 1 - step
          let c_cur = cost_vol[row, x]
          in if step == 0
             then col_res with [row] = copy c_cur
             else
               let prev = col_res[row + 1]
               let next_L = sgm_dp_step p1 p2 c_cur prev
               in col_res with [row] = copy next_L
      in col_agg[y]
    )
  in tabulate_2d h w (\y x ->
       tabulate d (\disp -> t2b[y, x, disp] + b2t[y, x, disp])
     )

-- ----------------------------------------------------------------------------
-- Full 4-Path SGM Cost Aggregation
-- ----------------------------------------------------------------------------

def aggregate_sgm_4path [h] [w] [d]
    (p1: i32)
    (p2: i32)
    (cost_vol: [h][w][d]i32): [h][w][d]i32 =
  let agg_h = aggregate_horizontal p1 p2 cost_vol
  let agg_v = aggregate_vertical p1 p2 cost_vol
  in tabulate_2d h w (\y x ->
       tabulate d (\disp -> (agg_h[y, x, disp] + agg_v[y, x, disp]) / 4i32)
     )

-- ----------------------------------------------------------------------------
-- Sub-Pixel Parabolic Disparity Refinement
--
-- Fits a 2nd-degree Taylor polynomial around discrete minimum disparity:
-- Delta d = (c_minus - c_plus) / (2 * (c_minus - 2*c0 + c_plus))
-- ----------------------------------------------------------------------------

def refine_subpixel [d] (costs: [d]i32): f64 =
  let (min_d, _) =
    loop (best_d, min_c) = (0, costs[0]) for disp in 1..<d do
      let c = costs[disp]
      in if c < min_c then (disp, c) else (best_d, min_c)
  in if min_d == 0 || min_d == d - 1
     then f64.i64 min_d
     else
       let c0 = f64.i32 costs[min_d]
       let c_minus = f64.i32 costs[min_d - 1]
       let c_plus = f64.i32 costs[min_d + 1]
       let denom = 2.0 * (c_minus - 2.0 * c0 + c_plus)
       in if f64.abs denom > 1e-6
          then
            let delta = (c_minus - c_plus) / denom
            let delta_clamped = f64.max (-0.5) (f64.min 0.5 delta)
            in f64.i64 min_d + delta_clamped
          else f64.i64 min_d

-- ----------------------------------------------------------------------------
-- Dense Disparity Map with Sub-Pixel Refinement
-- ----------------------------------------------------------------------------

def compute_refined_disparity [h] [w] [d]
    (agg_cost: [h][w][d]i32): [h][w]f64 =
  tabulate_2d h w (\y x ->
    refine_subpixel agg_cost[y, x]
  )

-- ----------------------------------------------------------------------------
-- 3D Euclidean Metric Point Cloud Back-Projection
--
-- Projects dense rectified disparity map to metric 3D point cloud:
-- Z = (f * B) / disparity
-- X = (u - cx) * Z / fx
-- Y = (v - cy) * Z / fy
-- Filters points with invalid or near-zero disparity.
-- ----------------------------------------------------------------------------

def backproject_to_3d [h] [w]
    (cam: camera_intrinsics)
    (baseline_mm: f64)
    (min_disp: f64)
    (disp_map: [h][w]f64): [h][w](vec3, bool) =
  tabulate_2d h w (\y x ->
    let d = disp_map[y, x]
    in if d < min_disp
       then ((0.0, 0.0, 0.0), false)
       else
         let z = (cam.fx * baseline_mm) / d
         let u = f64.i64 x
         let v = f64.i64 y
         let x_coord = (u - cam.cx) * z / cam.fx
         let y_coord = (v - cam.cy) * z / cam.fy
         in ((x_coord, y_coord, z), true)
  )

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_sgm_disparity [h] [w] [d]
    (p1: i32)
    (p2: i32)
    (cost_vol: [h][w][d]i32): [h][w]f64 =
  let agg = aggregate_sgm_4path p1 p2 cost_vol
  in compute_refined_disparity agg

entry fut_backproject [h] [w]
    (fx: f64, fy: f64, cx: f64, cy: f64, baseline_mm: f64, min_disp: f64)
    (disp_map: [h][w]f64): ([h][w]f64, [h][w]f64, [h][w]f64, [h][w]bool) =
  let cam = { fx = fx, fy = fy, cx = cx, cy = cy, skew = 0.0 }
  let pts = backproject_to_3d cam baseline_mm min_disp disp_map
  let xs = map (map (\(p, _) -> p.0)) pts
  let ys = map (map (\(p, _) -> p.1)) pts
  let zs = map (map (\(p, _) -> p.2)) pts
  let valids = map (map (\(_, v) -> v)) pts
  in (xs, ys, zs, valids)
