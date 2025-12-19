-- ============================================================================
-- Caliper Metrology Engine: Surface Normal & Curvature Estimation
-- Module: normals.fut
--
-- Data-parallel 3x3 covariance tensor estimation, analytic Cardano cubic
-- eigen-decomposition, surface curvature metrics, and viewpoint orientation.
-- ============================================================================

import "math"

-- ----------------------------------------------------------------------------
-- Covariance Tensor Estimation for a Local Point Cluster
--
-- Computes the symmetric 3x3 covariance matrix:
-- C = (1/K) * sum( (p_j - centroid) * (p_j - centroid)^T )
-- ----------------------------------------------------------------------------

def cluster_covariance [k] (pts: [k]vec3): (vec3, mat3) =
  if k == 0 then
    ((0.0, 0.0, 0.0), mat3_zeros)
  else
    -- Compute centroid
    let sum_p =
      reduce v3_add (0.0, 0.0, 0.0) pts
    let centroid = v3_scale (1.0 / f64.i64 k) sum_p
    -- Accumulate outer products
    let cov =
      loop acc = mat3_zeros for i < k do
        let diff = v3_sub pts[i] centroid
        let dx = diff.0
        let dy = diff.1
        let dz = diff.2
        let outer =
          [ [dx*dx, dx*dy, dx*dz]
          , [dy*dx, dy*dy, dy*dz]
          , [dz*dx, dz*dy, dz*dz]
          ]
        in mat3_add acc outer
    let inv_k = 1.0 / f64.i64 k
    let norm_cov = mat3_scale inv_k cov
    in (centroid, norm_cov)

-- ----------------------------------------------------------------------------
-- Surface Normal & Curvature from Local Covariance
--
-- Evaluates eigenvalues lambda1 <= lambda2 <= lambda3 via Cardano closed form.
-- Normal is eigenvector of lambda1. Curvature sigma = lambda1 / (lambda1 + lambda2 + lambda3).
-- Orients normal towards specified viewpoint.
-- ----------------------------------------------------------------------------

def normal_from_covariance
    (cov: mat3)
    (point: vec3)
    (viewpoint: vec3): (vec3, f64) =
  let (l1, l2, l3) = sym_mat3_eigenvalues cov
  let raw_normal = sym_mat3_min_eigenvector cov l1
  -- Orient normal towards sensor viewpoint: dot(normal, viewpoint - point) > 0
  let view_dir = v3_sub viewpoint point
  let normal =
    if v3_dot raw_normal view_dir < 0.0
    then v3_scale (-1.0) raw_normal
    else raw_normal
  let total_var = l1 + l2 + l3
  let curvature =
    if total_var > 1e-12
    then f64.max 0.0 (f64.min 1.0 (l1 / total_var))
    else 0.0
  in (normal, curvature)

-- ----------------------------------------------------------------------------
-- Point Cloud Normal Estimation across Stride-Sampled Neighborhoods
--
-- Computes normal and curvature for every point in a point cloud using
-- k-nearest or stride-sampled neighborhood window.
-- ----------------------------------------------------------------------------

def estimate_normals [n]
    (pts: [n]vec3)
    (k_neighbors: i64)
    (viewpoint: vec3): ([n]vec3, [n]f64) =
  let effective_k = i64.min n (i64.max 3 k_neighbors)
  let stride = i64.max 1 (n / effective_k)
  let results =
    tabulate n (\i ->
      let p_cur = pts[i]
      -- Gather neighboring points around index i
      let neighbors =
        tabulate effective_k (\k_idx ->
          let offset = k_idx - (effective_k / 2)
          let n_idx = (i + offset * stride + n * 100) % n
          in pts[n_idx]
        )
      let (_, cov) = cluster_covariance neighbors
      in normal_from_covariance cov p_cur viewpoint
    )
  let normals = map (\(norm, _) -> norm) results
  let curvatures = map (\(_, curv) -> curv) results
  in (normals, curvatures)

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_estimate_normals [n]
    (xs: [n]f64, ys: [n]f64, zs: [n]f64)
    (k: i64)
    (vx: f64, vy: f64, vz: f64): ([n]f64, [n]f64, [n]f64, [n]f64) =
  let pts = tabulate n (\i -> (xs[i], ys[i], zs[i]))
  let vp = (vx, vy, vz)
  let (norms, curvs) = estimate_normals pts k vp
  let nxs = map (\n_vec -> n_vec.0) norms
  let nys = map (\n_vec -> n_vec.1) norms
  let nzs = map (\n_vec -> n_vec.2) norms
  in (nxs, nys, nzs, curvs)
