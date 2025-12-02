-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Census Transform & Cost Volume
-- Module: tests/futhark/test_census.fut
--
-- Tests radiometric invariance of Census Transform and WTA disparity recovery
-- on synthetic stereo pairs with known shift.
-- ============================================================================

import "../../src/futhark/census"

def test_census_radiometric_invariance: bool =
  let h = 16i64
  let w = 16i64
  -- Base pattern with varied gradients
  let img1 = tabulate_2d h w (\y x -> f64.i64 (y * 10 + x * 3))
  -- Linearly scaled and offset pattern (simulating camera gain and bias)
  let img2 = tabulate_2d h w (\y x -> img1[y, x] * 2.5 + 45.0)
  let c1 = compute_census_map img1
  let c2 = compute_census_map img2
  -- The interior pixels should match exactly
  let matches =
    tabulate_2d h w (\y x ->
      if y >= 4 && y < h - 4 && x >= 3 && x < w - 3
      then c1[y, x] == c2[y, x]
      else true
    )
  in all (all id) matches

def test_cost_volume_and_wta: bool =
  let h = 16i64
  let w = 32i64
  let true_disp = 3i64
  -- Construct synthetic left image with high-frequency texture
  let left_img =
    tabulate_2d h w (\y x ->
      let fx = f64.i64 x
      let fy = f64.i64 y
      in f64.sin (fx * 1.2) * 50.0 + f64.cos (fy * 0.8) * 40.0 + 128.0
    )
  -- Construct right image shifted by true_disp
  let right_img =
    tabulate_2d h w (\y x ->
      let src_x = x + true_disp
      in if src_x < w then left_img[y, src_x] else 128.0
    )
  let c_left = compute_census_map left_img
  let c_right = compute_census_map right_img
  let max_disp = 8i64
  let cost_vol = compute_cost_volume c_left c_right max_disp
  let filtered = box_filter_cost_volume cost_vol 1
  let disp_map = wta_disparity filtered
  -- Interior pixels where texture exists should recover true_disp = 3
  let ok_count =
    loop count = 0 for y in 5...10 do
      loop inner = count for x in 8...20 do
        if disp_map[y, x] == i32.i64 true_disp then inner + 1 else inner
  -- Expect high recovery rate in interior
  in ok_count > 30

entry main: bool =
  test_census_radiometric_invariance && test_cost_volume_and_wta
