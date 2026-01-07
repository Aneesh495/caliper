-- ============================================================================
-- Caliper Metrology Engine: Signed Distance Fields (SDF) & Deviation Analysis
-- Module: distance.fut
--
-- Exact point-to-triangle Euclidean distance via 7-region Voronoi feature
-- classification, signed surface normal projection, micron-level deviation
-- statistics, and false-color RGB scalar field mapping.
-- ============================================================================

import "math"
import "mesh"
import "bvh"

type closest_result =
  { closest_point: vec3
  , distance: f64
  , signed_distance: f64
  , triangle_index: i64
  }

type deviation_stats =
  { mean_deviation: f64
  , std_dev: f64
  , rmse: f64
  , max_positive: f64
  , max_negative: f64
  , in_tolerance_pct: f64
  }

-- ----------------------------------------------------------------------------
-- Exact Point-to-Triangle Distance via Voronoi Feature Regions
--
-- David Eberly / Christer Ericson 7-region Voronoi classification.
-- Evaluates the exact closest point on triangle (A, B, C) to query point P.
-- ----------------------------------------------------------------------------

def point_triangle_closest (p: vec3) ((a, b, c): triangle): (vec3, f64) =
  let ab = v3_sub b a
  let ac = v3_sub c a
  let ap = v3_sub p a

  let d1 = v3_dot ab ap
  let d2 = v3_dot ac ap
  in if d1 <= 0.0 && d2 <= 0.0 then
       -- Vertex region A
       (a, v3_dist p a)
     else
       let bp = v3_sub p b
       let d3 = v3_dot ab bp
       let d4 = v3_dot ac bp
       in if d3 >= 0.0 && d4 <= d3 then
            -- Vertex region B
            (b, v3_dist p b)
          else
            let vc = d1*d4 - d3*d2
            in if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 then
                 -- Edge region AB
                 let v = d1 / (d1 - d3)
                 let q = v3_add a (v3_scale v ab)
                 in (q, v3_dist p q)
               else
                 let cp = v3_sub p c
                 let d5 = v3_dot ab cp
                 let d6 = v3_dot ac cp
                 in if d6 >= 0.0 && d5 <= d6 then
                      -- Vertex region C
                      (c, v3_dist p c)
                    else
                      let vb = d5*d2 - d1*d6
                      in if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 then
                           -- Edge region AC
                           let w = d2 / (d2 - d6)
                           let q = v3_add a (v3_scale w ac)
                           in (q, v3_dist p q)
                         else
                           let va = d3*d6 - d5*d4
                           in if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 then
                                -- Edge region BC
                                let w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
                                let bc = v3_sub c b
                                let q = v3_add b (v3_scale w bc)
                                in (q, v3_dist p q)
                              else
                                -- Face interior region
                                let denom = 1.0 / (va + vb + vc)
                                let v = vb * denom
                                let w = vc * denom
                                let q = v3_add a (v3_add (v3_scale v ab) (v3_scale w ac))
                                in (q, v3_dist p q)

-- ----------------------------------------------------------------------------
-- Signed Distance to Triangle Mesh
--
-- Evaluates closest triangle in mesh, projects displacement along surface normal:
-- sign = sgn( dot(P - Q, n_tri) )
-- ----------------------------------------------------------------------------

def point_mesh_signed_distance [m] (p: vec3) (triangles: [m]triangle): closest_result =
  if m == 0 then
    { closest_point = p, distance = 0.0, signed_distance = 0.0, triangle_index = -1 }
  else
    let (best_q, best_dist, best_idx) =
      loop (cur_q, cur_dist, cur_idx) = point_triangle_closest p triangles[0] |> (\(q, d) -> (q, d, 0i64))
      for i in 1..<m do
        let (q_i, d_i) = point_triangle_closest p triangles[i]
        in if d_i < cur_dist
           then (q_i, d_i, i)
           else (cur_q, cur_dist, cur_idx)

    let tri = triangles[best_idx]
    let n = triangle_normal tri
    let diff = v3_sub p best_q
    let sign = if v3_dot diff n >= 0.0 then 1.0 else -1.0
    in { closest_point = best_q
       , distance = best_dist
       , signed_distance = sign * best_dist
       , triangle_index = best_idx
       }

-- ----------------------------------------------------------------------------
-- Parallel Point Cloud to CAD Mesh Deviation Map
-- ----------------------------------------------------------------------------

def compute_deviations [n] [m]
    (pts: [n]vec3)
    (triangles: [m]triangle): [n]closest_result =
  map (\p -> point_mesh_signed_distance p triangles) pts

-- ----------------------------------------------------------------------------
-- Statistical Quality Analysis of Manufacturing Deviations
-- ----------------------------------------------------------------------------

def compute_deviation_stats [n]
    (results: [n]closest_result)
    (tolerance_mm: f64): deviation_stats =
  if n == 0 then
    { mean_deviation = 0.0
    , std_dev = 0.0
    , rmse = 0.0
    , max_positive = 0.0
    , max_negative = 0.0
    , in_tolerance_pct = 100.0
    }
  else
    let s_dists = map (\r -> r.signed_distance) results
    let sum_dev = reduce (+) 0.0 s_dists
    let mean_dev = sum_dev / f64.i64 n
    let sum_sq_dev = reduce (+) 0.0 (map (\d -> d * d) s_dists)
    let rmse = f64.sqrt (sum_sq_dev / f64.i64 n)
    let var_sum = reduce (+) 0.0 (map (\d -> (d - mean_dev)**2.0) s_dists)
    let std_dev = f64.sqrt (var_sum / f64.i64 n)
    let max_pos = reduce f64.max (-1e30) s_dists
    let max_neg = reduce f64.min (1e30) s_dists
    let in_tol_count =
      reduce (+) 0i64 (map (\d -> if f64.abs d <= tolerance_mm then 1i64 else 0i64) s_dists)
    let in_tol_pct = (f64.i64 in_tol_count / f64.i64 n) * 100.0
    in { mean_deviation = mean_dev
       , std_dev = std_dev
       , rmse = rmse
       , max_positive = max_pos
       , max_negative = max_neg
       , in_tolerance_pct = in_tol_pct
       }

-- ----------------------------------------------------------------------------
-- False-Color RGB Scalar Field Mapping
--
-- Maps signed deviation delta to RGB false color:
-- [-tol, +tol] -> Blue -> Cyan -> Green -> Yellow -> Red
-- ----------------------------------------------------------------------------

def deviation_to_rgb (delta: f64) (tolerance: f64): (f64, f64, f64) =
  let t = f64.max (-1.0) (f64.min 1.0 (delta / tolerance))
  -- Normalized to [0, 1]
  let s = (t + 1.0) * 0.5
  in if s < 0.25 then
       -- Blue (0, 0, 1) to Cyan (0, 1, 1)
       let w = s / 0.25
       in (0.0, w, 1.0)
     else if s < 0.5 then
       -- Cyan (0, 1, 1) to Green (0, 1, 0)
       let w = (s - 0.25) / 0.25
       in (0.0, 1.0, 1.0 - w)
     else if s < 0.75 then
       -- Green (0, 1, 0) to Yellow (1, 1, 0)
       let w = (s - 0.5) / 0.25
       in (w, 1.0, 0.0)
     else
       -- Yellow (1, 1, 0) to Red (1, 0, 0)
       let w = (s - 0.75) / 0.25
       in (1.0, 1.0 - w, 0.0)

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_compute_deviations [n] [m]
    (px: [n]f64, py: [n]f64, pz: [n]f64)
    (v0x: [m]f64, v0y: [m]f64, v0z: [m]f64)
    (v1x: [m]f64, v1y: [m]f64, v1z: [m]f64)
    (v2x: [m]f64, v2y: [m]f64, v2z: [m]f64): ([n]f64, [n]f64, [n]i64) =
  let pts = tabulate n (\i -> (px[i], py[i], pz[i]))
  let tris = tabulate m (\i ->
    ((v0x[i], v0y[i], v0z[i]), (v1x[i], v1y[i], v1z[i]), (v2x[i], v2y[i], v2z[i]))
  )
  let devs = compute_deviations pts tris
  let dists = map (\d -> d.distance) devs
  let s_dists = map (\d -> d.signed_distance) devs
  let indices = map (\d -> d.triangle_index) devs
  in (dists, s_dists, indices)
