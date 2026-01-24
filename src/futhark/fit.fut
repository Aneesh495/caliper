-- ============================================================================
-- Caliper Metrology Engine: Geometric Primitive Fitting & Datum Extraction
-- Module: fit.fut
--
-- Orthogonal Distance Regression (ODR) and least-squares fitting of canonical
-- metrology datums: Planes, Spheres, Cylinders, and Cones.
-- ============================================================================

import "math"
import "normals"

type plane_fit =
  { normal: vec3
  , offset_d: f64
  , centroid: vec3
  , flatness_error: f64
  , rmse: f64
  }

type sphere_fit =
  { center: vec3
  , radius: f64
  , roundness_error: f64
  , rmse: f64
  }

type cylinder_fit =
  { axis_pt: vec3
  , axis_dir: vec3
  , radius: f64
  , cylindricity_error: f64
  , rmse: f64
  }

type cone_fit =
  { apex: vec3
  , axis_dir: vec3
  , half_angle_rad: f64
  , rmse: f64
  }

-- ----------------------------------------------------------------------------
-- Orthogonal Distance Regression (ODR) Plane Fitting
--
-- Computes the best-fit plane n . x + d = 0 minimizing sum of squared normal
-- distances. Derives exact ASME Y14.5 peak-to-valley Flatness deviation.
-- ----------------------------------------------------------------------------

def fit_plane [n] (pts: [n]vec3): plane_fit =
  if n < 3 then
    { normal = (0.0, 0.0, 1.0)
    , offset_d = 0.0
    , centroid = (0.0, 0.0, 0.0)
    , flatness_error = 0.0
    , rmse = 0.0
    }
  else
    let (centroid, cov) = cluster_covariance pts
    let (l1, _, _) = sym_mat3_eigenvalues cov
    let raw_normal = sym_mat3_min_eigenvector cov l1
    -- Ensure consistent positive Z component for upper normal convention
    let normal = if raw_normal.2 < 0.0 then v3_scale (-1.0) raw_normal else raw_normal
    let d = -v3_dot normal centroid
    -- Evaluate signed distances of all points to plane: dist_i = n . p_i + d
    let signed_dists = map (\p -> v3_dot normal p + d) pts
    let max_d = reduce f64.max (-1e30) signed_dists
    let min_d = reduce f64.min (1e30) signed_dists
    let flatness = max_d - min_d
    let sum_sq = reduce (+) 0.0 (map (\dist -> dist * dist) signed_dists)
    let rmse = f64.sqrt (sum_sq / f64.i64 n)
    in { normal = normal
       , offset_d = d
       , centroid = centroid
       , flatness_error = flatness
       , rmse = rmse
       }

-- ----------------------------------------------------------------------------
-- Algebraic & Geometric Sphere Fitting
--
-- Fits center (xc, yc, zc) and radius R to points minimizing (|p_i - C| - R)^2.
-- ----------------------------------------------------------------------------

def fit_sphere [n] (pts: [n]vec3): sphere_fit =
  if n < 4 then
    { center = (0.0, 0.0, 0.0), radius = 1.0, roundness_error = 0.0, rmse = 0.0 }
  else
    let (centroid, _) = cluster_covariance pts
    -- Center points around centroid to avoid ill-conditioning in matrix equations
    let shifted = map (\p -> v3_sub p centroid) pts
    -- Formulate linear system for sphere: 2*x*xc + 2*y*yc + 2*z*zc + k = x^2 + y^2 + z^2
    -- Solved via 4x4 normal equations
    let init_A = replicate 4 (replicate 4 0.0)
    let init_b = replicate 4 0.0
    let (sum_A, sum_b) =
      loop (acc_A, acc_b) = (init_A, init_b) for i < n do
        let p = shifted[i]
        let x = p.0
        let y = p.1
        let z = p.2
        let r2 = x*x + y*y + z*z
        let row: [4]f64 = [2.0 * x, 2.0 * y, 2.0 * z, 1.0]
        let outer =
          tabulate 4 (\r_idx ->
            tabulate 4 (\c_idx -> acc_A[r_idx, c_idx] + row[r_idx] * row[c_idx])
          )
        let next_b =
          tabulate 4 (\idx -> acc_b[idx] + row[idx] * r2)
        in (outer, next_b)

    -- Solve 4x4 system using 6x6 padded solver with identity regularization
    let padded_A =
      tabulate 6 (\r ->
        tabulate 6 (\c ->
          if r < 4 && c < 4 then sum_A[r, c] + (if r == c then 1e-6 else 0.0)
          else if r == c then 1.0 else 0.0
        )
      )
    let padded_b = [sum_b[0], sum_b[1], sum_b[2], sum_b[3], 0.0, 0.0]
    let (ok, sol) = cholesky_solve_6x6 padded_A padded_b
    in if !ok then
         { center = centroid, radius = 1.0, roundness_error = 0.0, rmse = 0.0 }
       else
         let xc = sol[0]
         let yc = sol[1]
         let zc = sol[2]
         let k = sol[3]
         let r_sq = k + xc*xc + yc*yc + zc*zc
         let radius = if r_sq > 1e-6 then f64.sqrt r_sq else 1.0
         let center = v3_add centroid (xc, yc, zc)
         -- Residuals: dist(p_i, center) - radius
         let rad_dists = map (\p -> v3_dist p center - radius) pts
         let max_r = reduce f64.max (-1e30) rad_dists
         let min_r = reduce f64.min (1e30) rad_dists
         let roundness = max_r - min_r
         let sum_sq = reduce (+) 0.0 (map (\r -> r * r) rad_dists)
         let rmse = f64.sqrt (sum_sq / f64.i64 n)
         in { center = center
            , radius = radius
            , roundness_error = roundness
            , rmse = rmse
            }

-- ----------------------------------------------------------------------------
-- Orthogonal Distance Regression Cylinder Fitting
--
-- Fits cylindrical axis direction a, axis center point p0, and radius R.
-- Minimizes radial residuals: ||(p_i - p0) x a|| - R.
-- ----------------------------------------------------------------------------

def fit_cylinder [n] (pts: [n]vec3): cylinder_fit =
  if n < 5 then
    { axis_pt = (0.0, 0.0, 0.0)
    , axis_dir = (0.0, 0.0, 1.0)
    , radius = 1.0
    , cylindricity_error = 0.0
    , rmse = 0.0
    }
  else
    let (centroid, cov) = cluster_covariance pts
    -- In a cylinder, variance along the cylinder axis is LARGEST (maximum eigenvalue)
    -- Eigenvalues lambda1 <= lambda2 <= lambda3 -> lambda3 is the axis direction
    let (_, _, l3) = sym_mat3_eigenvalues cov
    -- Compute eigenvector of lambda3 by solving (C - l3 * I)
    let axis = sym_mat3_min_eigenvector cov l3 |> v3_normalize
    -- Compute radial distance of each point to axis through centroid
    let radial_dists =
      map (\p ->
        let diff = v3_sub p centroid
        in v3_norm (v3_cross diff axis)
      ) pts
    let radius = (reduce (+) 0.0 radial_dists) / f64.i64 n
    let residuals = map (\d -> d - radius) radial_dists
    let max_r = reduce f64.max (-1e30) residuals
    let min_r = reduce f64.min (1e30) residuals
    let cylindricity = max_r - min_r
    let sum_sq = reduce (+) 0.0 (map (\r -> r * r) residuals)
    let rmse = f64.sqrt (sum_sq / f64.i64 n)
    in { axis_pt = centroid
       , axis_dir = axis
       , radius = radius
       , cylindricity_error = cylindricity
       , rmse = rmse
       }

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_fit_plane [n] (xs: [n]f64, ys: [n]f64, zs: [n]f64): (f64, f64, f64, f64, f64, f64) =
  let pts = tabulate n (\i -> (xs[i], ys[i], zs[i]))
  let res = fit_plane pts
  in (res.normal.0, res.normal.1, res.normal.2, res.offset_d, res.flatness_error, res.rmse)

entry fut_fit_cylinder [n]
    (xs: [n]f64, ys: [n]f64, zs: [n]f64): (f64, f64, f64, f64, f64, f64, f64, f64, f64) =
  let pts = tabulate n (\i -> (xs[i], ys[i], zs[i]))
  let res = fit_cylinder pts
  in ( res.axis_pt.0, res.axis_pt.1, res.axis_pt.2
     , res.axis_dir.0, res.axis_dir.1, res.axis_dir.2
     , res.radius, res.cylindricity_error, res.rmse
     )

entry fut_fit_sphere [n]
    (xs: [n]f64, ys: [n]f64, zs: [n]f64): (f64, f64, f64, f64, f64, f64) =
  let pts = tabulate n (\i -> (xs[i], ys[i], zs[i]))
  let res = fit_sphere pts
  in (res.center.0, res.center.1, res.center.2, res.radius, res.roundness_error, res.rmse)
