-- ============================================================================
-- Caliper Metrology Engine: Point-to-Plane Iterative Closest Point (ICP)
-- Module: icp.fut
--
-- Data-parallel linear point-to-plane registration, linearized 6x6 normal
-- equation accumulation, Levenberg-Marquardt damping, Cholesky solving,
-- and robust Tukey biweight outlier rejection.
-- ============================================================================

import "math"

type icp_transform =
  { rot: mat3
  , trans: vec3
  , rmse: f64
  , converged: bool
  , iterations: i64
  }

-- ----------------------------------------------------------------------------
-- Point-to-Plane Linearized Correspondence Matrix & Vector Accumulation
--
-- Computes outer product a_i * a_i^T and scalar residual a_i * b_i for point pair:
-- a_i = [ (p_i x n_i), n_i ]^T in R^6
-- b_i = (q_i - p_i) . n_i in R
-- ----------------------------------------------------------------------------

def linear_correspondence
    (p: vec3)
    (q: vec3)
    (n: vec3): ([6][6]f64, [6]f64, f64) =
  let cross_pn = v3_cross p n
  let a: [6]f64 = [cross_pn.0, cross_pn.1, cross_pn.2, n.0, n.1, n.2]
  let diff = v3_sub q p
  let b = v3_dot diff n
  let H =
    tabulate 6 (\r ->
      tabulate 6 (\c -> a[r] * a[c])
    )
  let g = map (\ar -> ar * b) a
  in (H, g, b * b)

-- ----------------------------------------------------------------------------
-- Parallel Normal Equations Accumulation with Tukey Robust Weighting
-- ----------------------------------------------------------------------------

def accumulate_icp_system [k]
    (src_pts: [k]vec3)
    (tgt_pts: [k]vec3)
    (tgt_normals: [k]vec3)
    (tukey_c: f64)
    (damping: f64): ([6][6]f64, [6]f64, f64) =
  let init_H = replicate 6 (replicate 6 0.0)
  let init_g = replicate 6 0.0
  let (sum_H, sum_g, sum_err) =
    loop (acc_H, acc_g, acc_err) = (init_H, init_g, 0.0) for i < k do
      let p = src_pts[i]
      let q = tgt_pts[i]
      let n = tgt_normals[i]
      let (H_i, g_i, r2) = linear_correspondence p q n
      let r = f64.sqrt r2
      -- Tukey biweight: w = (1 - (r/c)^2)^2 if |r| < c else 0
      let w =
        if r < tukey_c then
          let u = r / tukey_c
          in (1.0 - u * u)**2.0
        else 0.0
      let w_H = map (map (\elem -> elem * w)) H_i
      let w_g = map (\elem -> elem * w) g_i
      let next_H =
        tabulate 6 (\r_idx ->
          tabulate 6 (\c_idx -> acc_H[r_idx, c_idx] + w_H[r_idx, c_idx])
        )
      let next_g =
        tabulate 6 (\idx -> acc_g[idx] + w_g[idx])
      in (next_H, next_g, acc_err + r2)

  -- Add Levenberg-Marquardt diagonal damping lambda * I
  let damped_H =
    tabulate 6 (\r ->
      tabulate 6 (\c ->
        if r == c then sum_H[r, c] + damping else sum_H[r, c]
      )
    )
  let rmse = f64.sqrt (sum_err / f64.i64 (i64.max 1 k))
  in (damped_H, sum_g, rmse)

-- ----------------------------------------------------------------------------
-- Single Step ICP Solve & SE(3) Transform Update
--
-- Solves H * delta = g -> delta = [alpha, beta, gamma, tx, ty, tz]
-- ----------------------------------------------------------------------------

def icp_step [k]
    (src_pts: [k]vec3)
    (tgt_pts: [k]vec3)
    (tgt_normals: [k]vec3)
    (damping: f64): (bool, mat3, vec3, f64) =
  let (H, g, rmse) = accumulate_icp_system src_pts tgt_pts tgt_normals 10.0 damping
  let (ok, delta) = cholesky_solve_6x6 H g
  in if !ok then
       (false, mat3_identity, (0.0, 0.0, 0.0), rmse)
     else
       let alpha = delta[0]
       let beta = delta[1]
       let gamma = delta[2]
       let tx = delta[3]
       let ty = delta[4]
       let tz = delta[5]
       -- Small angle rotation matrix approximation
       let dR: mat3 =
         [ [1.0, -gamma, beta]
         , [gamma, 1.0, -alpha]
         , [-beta, alpha, 1.0]
         ]
       -- Orthonormalize dR via Gram-Schmidt
       let c0: vec3 = v3_normalize (dR[0,0], dR[1,0], dR[2,0])
       let c1_raw: vec3 = (dR[0,1], dR[1,1], dR[2,1])
       let c1 = v3_normalize (v3_sub c1_raw (v3_scale (v3_dot c1_raw c0) c0))
       let c2 = v3_cross c0 c1
       let R_ortho: mat3 =
         [ [c0.0, c1.0, c2.0]
         , [c0.1, c1.1, c2.1]
         , [c0.2, c1.2, c2.2]
         ]
       let t: vec3 = (tx, ty, tz)
       in (true, R_ortho, t, rmse)

-- ----------------------------------------------------------------------------
-- Multi-Iteration Point-to-Plane ICP Registration Loop
-- ----------------------------------------------------------------------------

def run_point_to_plane_icp [k]
    (src_pts: [k]vec3)
    (tgt_pts: [k]vec3)
    (tgt_normals: [k]vec3)
    (max_iters: i64)
    (tol: f64): icp_transform =
  let (cur_R, cur_t, _, final_rmse, conv, iters) =
    loop (R, t, pts, prev_rmse, converged, iter) =
      (mat3_identity, ((0.0, 0.0, 0.0): vec3), src_pts, 1e30, false, 0i64)
    for _step < max_iters do
      if converged then
        (R, t, pts, prev_rmse, true, iter)
      else
        let (ok, dR, dt, step_rmse) = icp_step pts tgt_pts tgt_normals 1e-4
        in if !ok then
             (R, t, pts, prev_rmse, false, iter)
           else
             let next_R = mat3_mul dR R
             let next_t = v3_add (mat3_mul_v3 dR t) dt
             let next_pts = map (\p -> transform_point dR dt p) pts
             let delta_rmse = f64.abs (prev_rmse - step_rmse)
             let is_conv = delta_rmse < tol || v3_norm dt < 1e-6
             in (next_R, next_t, next_pts, step_rmse, is_conv, iter + 1)

  in { rot = cur_R
     , trans = cur_t
     , rmse = final_rmse
     , converged = conv
     , iterations = iters
     }

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_icp_step [k]
    (px: [k]f64, py: [k]f64, pz: [k]f64)
    (qx: [k]f64, qy: [k]f64, qz: [k]f64)
    (nx: [k]f64, ny: [k]f64, nz: [k]f64): (bool, [3][3]f64, f64, f64, f64, f64) =
  let p = tabulate k (\i -> (px[i], py[i], pz[i]))
  let q = tabulate k (\i -> (qx[i], qy[i], qz[i]))
  let n = tabulate k (\i -> (nx[i], ny[i], nz[i]))
  let (ok, dR, dt, rmse) = icp_step p q n 1e-4
  in (ok, dR, dt.0, dt.1, dt.2, rmse)
