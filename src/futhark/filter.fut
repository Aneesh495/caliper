-- ============================================================================
-- Caliper Metrology Engine: Point Cloud Filtering & Outlier Removal
-- Module: filter.fut
--
-- Data-parallel spatial voxel grid filtering, Statistical Outlier Removal (SOR),
-- and normal-aligned bilateral point cloud smoothing.
-- ============================================================================

import "math"

-- ----------------------------------------------------------------------------
-- Statistical Outlier Removal (SOR)
--
-- Computes the mean Euclidean distance of each point to its k nearest neighbors
-- (or local spatial sample), derives global Gaussian distribution (mu, sigma),
-- and flags points beyond (mu + alpha * sigma) as noise artifacts.
-- ----------------------------------------------------------------------------

def sor_filter [n] (pts: [n]vec3) (k_sample: i64) (alpha: f64): [n]bool =
  if n <= 1 then replicate n true
  else
    let effective_k = i64.min (n - 1) (i64.max 1 k_sample)
    -- Compute average distance of each point to sampled neighboring points
    let mean_dists =
      map (\p ->
        -- Compute sorted distances to a window of points
        let stride = i64.max 1 (n / effective_k)
        let sampled =
          tabulate effective_k (\idx ->
            let target_idx = (idx * stride) % n
            in v3_dist p pts[target_idx]
          )
        let sum_d = reduce (+) 0.0 sampled
        in sum_d / f64.i64 effective_k
      ) pts

    let global_mu = (reduce (+) 0.0 mean_dists) / f64.i64 n
    let var_sum =
      reduce (+) 0.0 (map (\d -> (d - global_mu)**2.0) mean_dists)
    let global_sigma = f64.sqrt (var_sum / f64.i64 n)
    let threshold = global_mu + alpha * global_sigma

    in map (\d -> d <= threshold) mean_dists

-- ----------------------------------------------------------------------------
-- Spatial Voxel Grid Downsampling via Uniform Bounding Box Discretization
--
-- Maps points into regular 3D grid bins and calculates per-voxel centroids.
-- ----------------------------------------------------------------------------

def compute_bounds [n] (pts: [n]vec3): (vec3, vec3) =
  if n == 0 then ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
  else
    let init_min = pts[0]
    let init_max = pts[0]
    let (bmin, bmax) =
      reduce (\(min_a, max_a) (min_b, max_b) ->
        (v3_min min_a min_b, v3_max max_a max_b)
      ) (init_min, init_max) (zip pts pts)
    in (bmin, bmax)

def voxel_hash (voxel_size: f64) (bmin: vec3) (p: vec3): (i64, i64, i64) =
  let vx = i64.f64 (f64.floor ((p.0 - bmin.0) / voxel_size))
  let vy = i64.f64 (f64.floor ((p.1 - bmin.1) / voxel_size))
  let vz = i64.f64 (f64.floor ((p.2 - bmin.2) / voxel_size))
  in (vx, vy, vz)

-- ----------------------------------------------------------------------------
-- Bilateral Point Cloud Smoothing
--
-- Preserves high-curvature geometric edges while suppressing sensor noise:
-- delta_p = sum( Ws(|p - q|) * Wr((q - p) . n) * ((q - p) . n) ) / sum( Ws * Wr )
-- ----------------------------------------------------------------------------

def bilateral_smooth_point [n]
    (p: vec3)
    (n_vec: vec3)
    (pts: [n]vec3)
    (sigma_s: f64)
    (sigma_r: f64): vec3 =
  let inv_2_sig_s2 = 1.0 / (2.0 * sigma_s * sigma_s)
  let inv_2_sig_r2 = 1.0 / (2.0 * sigma_r * sigma_r)
  let (weighted_disp, total_weight) =
    loop (w_disp, w_sum) = (0.0, 0.0) for j < n do
      let q = pts[j]
      let diff = v3_sub q p
      let dist_sq = v3_norm_sq diff
      in if dist_sq > 9.0 * sigma_s * sigma_s then
           (w_disp, w_sum)
         else
           let proj = v3_dot diff n_vec
           let ws = f64.exp (-dist_sq * inv_2_sig_s2)
           let wr = f64.exp (-(proj * proj) * inv_2_sig_r2)
           let w = ws * wr
           in (w_disp + w * proj, w_sum + w)

  in if total_weight > 1e-9
     then v3_add p (v3_scale (weighted_disp / total_weight) n_vec)
     else p

def bilateral_filter [n]
    (pts: [n]vec3)
    (normals: [n]vec3)
    (sigma_s: f64)
    (sigma_r: f64): [n]vec3 =
  tabulate n (\i ->
    bilateral_smooth_point pts[i] normals[i] pts sigma_s sigma_r
  )

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_sor_filter [n] (xs: [n]f64, ys: [n]f64, zs: [n]f64, k: i64, alpha: f64): [n]bool =
  let pts = tabulate n (\i -> (xs[i], ys[i], zs[i]))
  in sor_filter pts k alpha

entry fut_bilateral_filter [n]
    (xs: [n]f64, ys: [n]f64, zs: [n]f64)
    (nxs: [n]f64, nys: [n]f64, nzs: [n]f64)
    (sigma_s: f64, sigma_r: f64): ([n]f64, [n]f64, [n]f64) =
  let pts = tabulate n (\i -> (xs[i], ys[i], zs[i]))
  let norms = tabulate n (\i -> (nxs[i], nys[i], nzs[i]))
  let smoothed = bilateral_filter pts norms sigma_s sigma_r
  let out_x = map (\p -> p.0) smoothed
  let out_y = map (\p -> p.1) smoothed
  let out_z = map (\p -> p.2) smoothed
  in (out_x, out_y, out_z)
