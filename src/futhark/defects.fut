-- ============================================================================
-- Caliper Metrology Engine: Defect Clustering & Morphological Classification
-- Module: defects.fut
--
-- Data-parallel spatial clustering of out-of-tolerance scan points,
-- geometric bounding, morphological feature extraction, and automated
-- manufacturing defect classification (Burr/Flash, Sink/Void, Warpage, Chatter).
-- ============================================================================

import "math"
import "normals"

-- Defect class constants
-- 0: IN_TOLERANCE
-- 1: EXCESS_MATERIAL_FLASH_BURR
-- 2: MISSING_MATERIAL_SINK_VOID
-- 3: SURFACE_WARPAGE_TWIST
-- 4: MACHINING_CHATTER_RIPPLE

type defect_cluster =
  { cluster_id: i64
  , point_count: i64
  , centroid: vec3
  , peak_deviation: f64
  , mean_deviation: f64
  , bounding_extent: f64
  , classification: i32
  }

-- ----------------------------------------------------------------------------
-- Identification of Out-of-Tolerance Scan Points
-- ----------------------------------------------------------------------------

def flag_out_of_tolerance [n]
    (signed_deviations: [n]f64)
    (upper_tol_mm: f64)
    (lower_tol_mm: f64): [n]bool =
  map (\d -> d > upper_tol_mm || d < lower_tol_mm) signed_deviations

-- ----------------------------------------------------------------------------
-- Spatial Distance Matrix for Defect Point Grouping
--
-- Computes connected component connectivity graph for points where dist < eps.
-- ----------------------------------------------------------------------------

def compute_cluster_adjacency [k]
    (defect_pts: [k]vec3)
    (cluster_radius_mm: f64): [k][k]bool =
  let r_sq = cluster_radius_mm * cluster_radius_mm
  in tabulate_2d k k (\i j ->
       v3_dist_sq defect_pts[i] defect_pts[j] <= r_sq
     )

-- ----------------------------------------------------------------------------
-- Connected Component Cluster Assignment via Parallel Label Propagation
-- ----------------------------------------------------------------------------

def label_propagation [k]
    (adj: [k][k]bool)
    (max_iters: i64): [k]i64 =
  let init_labels = tabulate k (\i -> i)
  in loop labels = init_labels for _step < max_iters do
       tabulate k (\i ->
         loop min_label = labels[i] for j < k do
           if adj[i, j] then i64.min min_label labels[j] else min_label
       )

-- ----------------------------------------------------------------------------
-- Morphological Classification of a Defect Cluster
-- ----------------------------------------------------------------------------

def classify_cluster_morphology [k]
    (cluster_pts: [k]vec3)
    (cluster_devs: [k]f64)
    (cluster_id: i64): defect_cluster =
  if k == 0 then
    { cluster_id = cluster_id
    , point_count = 0
    , centroid = (0.0, 0.0, 0.0)
    , peak_deviation = 0.0
    , mean_deviation = 0.0
    , bounding_extent = 0.0
    , classification = 0
    }
  else
    let sum_p = reduce v3_add (0.0, 0.0, 0.0) cluster_pts
    let centroid = v3_scale (1.0 / f64.i64 k) sum_p
    let sum_d = reduce (+) 0.0 cluster_devs
    let mean_d = sum_d / f64.i64 k
    let peak_d =
      loop best = cluster_devs[0] for i in 1..<k do
        if f64.abs cluster_devs[i] > f64.abs best then cluster_devs[i] else best

    -- Bounding diameter
    let dists_to_c = map (\p -> v3_dist p centroid) cluster_pts
    let max_radius = reduce f64.max 0.0 dists_to_c
    let extent = 2.0 * max_radius

    -- Covariance PCA for elongation
    let (_, cov) = cluster_covariance cluster_pts
    let (_, l2, l3) = sym_mat3_eigenvalues cov
    let elongation = if l2 > 1e-9 then l3 / l2 else 1.0

    -- Classification heuristics:
    -- If peak > 0 and high elongation -> Flash / Burr along edge
    -- If peak < 0 -> Sink / Porosity void
    -- If extent is broad with low variance -> Warpage
    let class_id =
      if peak_d > 0.0 && elongation > 3.0 then 1i32 -- FLASH_BURR
      else if peak_d < 0.0 then 2i32                 -- SINK_VOID
      else if extent > 20.0 then 3i32                -- WARPAGE
      else 4i32                                      -- CHATTER

    in { cluster_id = cluster_id
       , point_count = k
       , centroid = centroid
       , peak_deviation = peak_d
       , mean_deviation = mean_d
       , bounding_extent = extent
       , classification = class_id
       }

-- ----------------------------------------------------------------------------
-- Automated Factory Quality Disposition Gate
--
-- Evaluates whether inspected part is PASS, REWORKABLE (excess material),
-- or SCRAP (missing material below critical tolerance).
-- Return: 0 = PASS, 1 = REWORK, 2 = SCRAP
-- ----------------------------------------------------------------------------

def evaluate_quality_disposition [k]
    (clusters: [k]defect_cluster)
    (max_allowable_clusters: i64)
    (critical_scrap_neg_tol_mm: f64): i32 =
  if k == 0 then 0i32 -- Zero defects -> PASS
  else if k > max_allowable_clusters then 2i32 -- Excessive defects -> SCRAP
  else
    let has_fatal_sink =
      reduce (||) false
        (map (\c -> c.classification == 2i32 && c.peak_deviation < critical_scrap_neg_tol_mm) clusters)
    in if has_fatal_sink then 2i32 -- Unrecoverable void -> SCRAP
       else
         let all_excess =
           reduce (&&) true (map (\c -> c.peak_deviation > 0.0) clusters)
         in if all_excess then 1i32 -- Excess material can be deburred/machined -> REWORK
            else 2i32               -- Non-repairable -> SCRAP

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_flag_defects [n]
    (signed_devs: [n]f64)
    (upper_tol: f64)
    (lower_tol: f64): [n]bool =
  flag_out_of_tolerance signed_devs upper_tol lower_tol

entry fut_cluster_labels [k]
    (px: [k]f64, py: [k]f64, pz: [k]f64)
    (radius: f64)
    (max_iters: i64): [k]i64 =
  let pts = tabulate k (\i -> (px[i], py[i], pz[i]))
  let adj = compute_cluster_adjacency pts radius
  in label_propagation adj max_iters
