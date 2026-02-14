-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Defect Clustering & Classification
-- Module: tests/futhark/test_defects.fut
--
-- Tests spatial connected-component label propagation, morphological defect
-- classification, and automated manufacturing quality disposition gates.
-- ============================================================================

import "../../src/futhark/defects"
import "../../src/futhark/math"

def test_spatial_clustering: bool =
  -- Cluster 1: 5 points around (0, 0, 0) within distance 1.0
  let c1: [5]vec3 =
    [ (0.0, 0.0, 0.0)
    , (0.5, 0.2, 0.1)
    , (-0.4, 0.3, -0.2)
    , (0.1, -0.5, 0.3)
    , (0.2, 0.4, -0.1)
    ]
  -- Cluster 2: 5 points around (50, 50, 50) within distance 1.0
  let c2: [5]vec3 =
    [ (50.0, 50.0, 50.0)
    , (50.3, 50.1, 49.8)
    , (49.7, 50.4, 50.2)
    , (50.2, 49.6, 50.1)
    , (49.9, 50.2, 49.9)
    ]
  let all_pts = concat c1 c2
  let adj = compute_cluster_adjacency all_pts 2.0
  let labels = label_propagation adj 10
  -- All points in cluster 1 should have identical label
  let l1 = labels[0]
  let c1_ok = all (\idx -> labels[idx] == l1) (0...4)
  -- All points in cluster 2 should have identical label
  let l2 = labels[5]
  let c2_ok = all (\idx -> labels[idx] == l2) (5...9)
  -- Two clusters must have distinct labels
  let distinct_ok = l1 != l2
  in c1_ok && c2_ok && distinct_ok

def test_morphology_classification: bool =
  -- Elongated thin ridge of excess material (Burr / Flash): X spans 20 mm, Y spans 0.2 mm
  let pts: [6]vec3 =
    [ (0.0, 0.0, 0.0)
    , (4.0, 0.05, 0.0)
    , (8.0, -0.05, 0.0)
    , (12.0, 0.02, 0.0)
    , (16.0, -0.02, 0.0)
    , (20.0, 0.0, 0.0)
    ]
  let devs: [6]f64 = [0.12, 0.15, 0.18, 0.14, 0.13, 0.11] -- Positive excess material
  let cluster = classify_cluster_morphology pts devs 0
  -- Should be classified as 1 (FLASH_BURR)
  let ok_class = cluster.classification == 1i32
  let ok_extent = cluster.bounding_extent >= 19.0
  let ok_peak = f64.abs (cluster.peak_deviation - 0.18) < 1e-6
  in ok_class && ok_extent && ok_peak

def test_quality_disposition: bool =
  let dummy_c =
    { cluster_id = 0
    , point_count = 5
    , centroid = ((0.0, 0.0, 0.0): vec3)
    , peak_deviation = 0.080
    , mean_deviation = 0.060
    , bounding_extent = 5.0
    , classification = 1i32 -- Excess material
    }
  -- 1 excess material defect -> REWORKABLE (1)
  let disp_rework = evaluate_quality_disposition [dummy_c] 5 (-0.100)
  -- Fatal negative void defect (-0.250 mm deep) -> SCRAP (2)
  let void_c =
    { cluster_id = 1
    , point_count = 5
    , centroid = ((0.0, 0.0, 0.0): vec3)
    , peak_deviation = -0.250
    , mean_deviation = -0.200
    , bounding_extent = 5.0
    , classification = 2i32
    }
  let disp_scrap = evaluate_quality_disposition [void_c] 5 (-0.100)
  in disp_rework == 1i32 && disp_scrap == 2i32

entry main: bool =
  test_spatial_clustering && test_morphology_classification && test_quality_disposition
