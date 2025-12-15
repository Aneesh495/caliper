-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Point Cloud Filtering
-- Module: tests/futhark/test_filter.fut
--
-- Tests Statistical Outlier Removal (SOR) on known outliers and bilateral
-- smoothing variance reduction along surface normals.
-- ============================================================================

import "../../src/futhark/filter"
import "../../src/futhark/math"

def test_sor_outlier_detection: bool =
  -- Create a cluster of 50 points centered around (0, 0, 0) within radius 1.0
  let cluster =
    tabulate 50 (\i ->
      let fi = f64.i64 i
      in (f64.sin fi * 0.8, f64.cos fi * 0.8, (fi / 50.0) * 0.5)
    )
  -- Add 2 distant outliers at (50.0, 50.0, 50.0) and (-40.0, -40.0, -40.0)
  let outliers: [2]vec3 = [(50.0, 50.0, 50.0), (-40.0, -40.0, -40.0)]
  let all_pts = concat cluster outliers
  let inlier_flags = sor_filter all_pts 10 1.5
  -- Cluster points should be true, outliers should be false
  let cluster_ok = all id inlier_flags[0:50]
  let outlier1_bad = !inlier_flags[50]
  let outlier2_bad = !inlier_flags[51]
  in cluster_ok && outlier1_bad && outlier2_bad

def test_bilateral_smoothing: bool =
  -- Create points on plane Z = 0 with random-like noise in Z
  let n = 20i64
  let pts =
    tabulate n (\i ->
      let fi = f64.i64 i
      let noise_z = if i % 2 == 0 then 0.5 else -0.5
      in (fi * 0.5, 0.0, noise_z)
    )
  let normals = replicate n ((0.0, 0.0, 1.0): vec3)
  let smoothed = bilateral_filter pts normals 2.0 1.0
  -- Check that Z coordinates moved closer to 0.0 (variance reduction)
  let initial_z_err = reduce (+) 0.0 (map (\p -> f64.abs p.2) pts)
  let smoothed_z_err = reduce (+) 0.0 (map (\p -> f64.abs p.2) smoothed)
  in smoothed_z_err < initial_z_err

entry main: bool =
  test_sor_outlier_detection && test_bilateral_smoothing
