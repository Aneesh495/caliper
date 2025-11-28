-- ============================================================================
-- Caliper Metrology Engine: Census Transform & Stereo Cost Volume
-- Module: census.fut
--
-- High-throughput bitwise Census Transform (9x7 window, 63 bits) and
-- parallel Hamming distance 3D cost volume construction for stereo matching.
-- ============================================================================

-- ----------------------------------------------------------------------------
-- 9x7 Window Census Transform Kernel
--
-- Evaluates 63 neighbor offsets around center pixel (y, x).
-- Produces a 64-bit integer where bit k is 1 if neighbor < center.
-- Invariant to monotonic radiometric gain and bias transformations.
-- ----------------------------------------------------------------------------

def census_pixel_9x7 [h] [w] (img: [h][w]f64) (y: i64) (x: i64): u64 =
  let center_val = img[y, x]
  -- Window dimensions: radius_y = 4, radius_x = 3 (height 9, width 7)
  let code =
    loop acc = 0u64 for dy_idx < 9 do
      let dy = dy_idx - 4
      let py = i64.max 0 (i64.min (h - 1) (y + dy))
      in loop inner_acc = acc for dx_idx < 7 do
        let dx = dx_idx - 3
        let px = i64.max 0 (i64.min (w - 1) (x + dx))
        in if dy == 0 && dx == 0 then
             inner_acc
           else
             let neighbor_val = img[py, px]
             let bit = if neighbor_val < center_val then 1u64 else 0u64
             in (inner_acc << 1u64) | bit
  in code

-- ----------------------------------------------------------------------------
-- Parallel 2D Census Map Computation
-- ----------------------------------------------------------------------------

def compute_census_map [h] [w] (img: [h][w]f64): [h][w]u64 =
  tabulate_2d h w (\y x -> census_pixel_9x7 img y x)

-- ----------------------------------------------------------------------------
-- Parallel 3D Stereo Cost Volume Construction
--
-- Computes the Hamming distance cost volume across all pixels and
-- candidate disparity hypotheses d in [0, max_disp - 1].
-- ----------------------------------------------------------------------------

def compute_cost_volume [h] [w]
    (left_census: [h][w]u64)
    (right_census: [h][w]u64)
    (max_disp: i64): [h][w][]i32 =
  tabulate_2d h w (\y x ->
    tabulate max_disp (\d ->
      let right_x = x - d
      in if right_x >= 0 && right_x < w
         then
           let c_left = left_census[y, x]
           let c_right = right_census[y, right_x]
           in u64.popc (c_left ^ c_right)
         else
           64i32 -- Maximum penalty for occluded / out of frame hypotheses
    )
  )

-- ----------------------------------------------------------------------------
-- Box Filter Cost Volume Aggregation
--
-- Applies a 2D separable box filter of radius r across the cost volume
-- to integrate spatial support before dynamic programming path aggregation.
-- ----------------------------------------------------------------------------

def box_filter_cost_volume [h] [w] [d]
    (cost_vol: [h][w][d]i32)
    (radius: i64): [h][w][d]i32 =
  -- Horizontal pass
  let horiz =
    tabulate_2d h w (\y x ->
      tabulate d (\disp ->
        let x_start = i64.max 0 (x - radius)
        let x_end = i64.min (w - 1) (x + radius)
        let count = x_end - x_start + 1
        let sum = loop s = 0 for k in x_start...x_end do s + cost_vol[y, k, disp]
        in sum / i32.i64 count
      )
    )
  -- Vertical pass
  in tabulate_2d h w (\y x ->
       tabulate d (\disp ->
         let y_start = i64.max 0 (y - radius)
         let y_end = i64.min (h - 1) (y + radius)
         let count = y_end - y_start + 1
         let sum = loop s = 0 for k in y_start...y_end do s + horiz[k, x, disp]
         in sum / i32.i64 count
       )
     )

-- ----------------------------------------------------------------------------
-- Winner-Takes-All (WTA) Baseline Disparity Selection
--
-- Extracts the argmin disparity hypothesis for each pixel from the cost volume.
-- ----------------------------------------------------------------------------

def wta_disparity [h] [w] [d] (cost_vol: [h][w][d]i32): [h][w]i32 =
  tabulate_2d h w (\y x ->
    let costs = cost_vol[y, x]
    let (best_d, _) =
      loop (min_d, min_c) = (0, costs[0]) for disp in 1..<d do
        let c = costs[disp]
        in if c < min_c then (i32.i64 disp, c) else (min_d, min_c)
    in best_d
  )

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_census_map [h] [w] (img: [h][w]f64): [h][w]u64 =
  compute_census_map img

entry fut_cost_volume [h] [w]
    (left_img: [h][w]f64)
    (right_img: [h][w]f64)
    (max_disp: i64): [h][w][]i32 =
  let c_left = compute_census_map left_img
  let c_right = compute_census_map right_img
  in compute_cost_volume c_left c_right max_disp
