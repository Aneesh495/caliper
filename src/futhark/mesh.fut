-- ============================================================================
-- Caliper Metrology Engine: Surface Polygonizer & Mesh Extraction
-- Module: mesh.fut
--
-- Data-parallel Marching Cubes / Tetrahedra isosurface extraction from 3D scalar
-- fields, edge linear interpolation, triangle face generation, outward normal
-- orientation, and mesh metrics (watertight volume, surface area, and bounding boxes).
-- ============================================================================

import "math"

type triangle = (vec3, vec3, vec3)

type mesh_stats =
  { num_triangles: i64
  , surface_area: f64
  , volume: f64
  , bbox_min: vec3
  , bbox_max: vec3
  }

-- ----------------------------------------------------------------------------
-- Triangle Area and Volume Primitives
-- ----------------------------------------------------------------------------

def triangle_normal ((v0, v1, v2): triangle): vec3 =
  let e1 = v3_sub v1 v0
  let e2 = v3_sub v2 v0
  in v3_normalize (v3_cross e1 e2)

def triangle_area ((v0, v1, v2): triangle): f64 =
  let e1 = v3_sub v1 v0
  let e2 = v3_sub v2 v0
  in 0.5 * v3_norm (v3_cross e1 e2)

def triangle_signed_volume ((v0, v1, v2): triangle): f64 =
  -- Divergence theorem contribution: (1/6) * dot(v0, cross(v1, v2))
  (1.0 / 6.0) * v3_dot v0 (v3_cross v1 v2)

-- ----------------------------------------------------------------------------
-- Edge Linear Interpolation for Isosurface Crossing
--
-- Finds crossing vertex p = p1 + ((isoval - val1) / (val2 - val1)) * (p2 - p1)
-- ----------------------------------------------------------------------------

def interpolate_edge
    (isoval: f64)
    (p1: vec3) (val1: f64)
    (p2: vec3) (val2: f64): vec3 =
  let denom = val2 - val1
  in if f64.abs denom < 1e-12
     then v3_lerp p1 p2 0.5
     else
       let t = (isoval - val1) / denom
       let t_clamped = f64.max 0.0 (f64.min 1.0 t)
       in v3_lerp p1 p2 t_clamped

-- ----------------------------------------------------------------------------
-- Outward Normal Alignment
--
-- Orients vertices so normal points from inside vertex towards triangle centroid.
-- ----------------------------------------------------------------------------

def orient_triangle_outward (tri: triangle) (p_inside: vec3): triangle =
  let (v0, v1, v2) = tri
  let centroid = v3_scale (1.0 / 3.0) (v3_add v0 (v3_add v1 v2))
  let out_dir = v3_sub centroid p_inside
  let n = v3_cross (v3_sub v1 v0) (v3_sub v2 v0)
  in if v3_dot n out_dir < 0.0
     then (v0, v2, v1)
     else (v0, v1, v2)

-- ----------------------------------------------------------------------------
-- Marching Tetrahedra Decomposition of a Cubic Voxel
--
-- Decomposes each cubic cell into 6 tetrahedra. Unlike classic Marching Cubes,
-- Marching Tetrahedra has only 16 cases and zero topological ambiguity,
-- guaranteeing completely manifold, watertight closed surfaces.
-- Emits fixed [2]triangle with valid count (0, 1, or 2).
-- ----------------------------------------------------------------------------

def polygonize_tetrahedron
    (isoval: f64)
    (p0: vec3) (v0: f64)
    (p1: vec3) (v1: f64)
    (p2: vec3) (v2: f64)
    (p3: vec3) (v3: f64): ([2]triangle, i32) =
  let b0 = if v0 < isoval then 1 else 0
  let b1 = if v1 < isoval then 2 else 0
  let b2 = if v2 < isoval then 4 else 0
  let b3 = if v3 < isoval then 8 else 0
  let index = b0 | b1 | b2 | b3
  -- Edge intersections
  let e01 = interpolate_edge isoval p0 v0 p1 v1
  let e02 = interpolate_edge isoval p0 v0 p2 v2
  let e03 = interpolate_edge isoval p0 v0 p3 v3
  let e12 = interpolate_edge isoval p1 v1 p2 v2
  let e13 = interpolate_edge isoval p1 v1 p3 v3
  let e23 = interpolate_edge isoval p2 v2 p3 v3
  let dummy_tri = (p0, p0, p0)
  in match index
     case 0  -> ([dummy_tri, dummy_tri], 0)
     case 15 -> ([dummy_tri, dummy_tri], 0)
     -- 1 vertex inside
     case 1  -> ([orient_triangle_outward (e01, e02, e03) p0, dummy_tri], 1)
     case 14 -> ([orient_triangle_outward (e01, e03, e02) p0, dummy_tri], 1)
     case 2  -> ([orient_triangle_outward (e01, e13, e12) p1, dummy_tri], 1)
     case 13 -> ([orient_triangle_outward (e01, e12, e13) p1, dummy_tri], 1)
     case 4  -> ([orient_triangle_outward (e02, e12, e23) p2, dummy_tri], 1)
     case 11 -> ([orient_triangle_outward (e02, e23, e12) p2, dummy_tri], 1)
     case 8  -> ([orient_triangle_outward (e03, e23, e13) p3, dummy_tri], 1)
     case 7  -> ([orient_triangle_outward (e03, e13, e23) p3, dummy_tri], 1)
     -- 2 vertices inside (quadrilateral split into 2 triangles)
     case 3  ->
       let p_in = v3_lerp p0 p1 0.5
       in ([orient_triangle_outward (e02, e12, e03) p_in, orient_triangle_outward (e12, e13, e03) p_in], 2)
     case 12 ->
       let p_in = v3_lerp p2 p3 0.5
       in ([orient_triangle_outward (e02, e03, e12) p_in, orient_triangle_outward (e12, e03, e13) p_in], 2)
     case 5  ->
       let p_in = v3_lerp p0 p2 0.5
       in ([orient_triangle_outward (e01, e12, e03) p_in, orient_triangle_outward (e12, e23, e03) p_in], 2)
     case 10 ->
       let p_in = v3_lerp p1 p3 0.5
       in ([orient_triangle_outward (e01, e03, e12) p_in, orient_triangle_outward (e12, e03, e23) p_in], 2)
     case 6  ->
       let p_in = v3_lerp p1 p2 0.5
       in ([orient_triangle_outward (e01, e02, e13) p_in, orient_triangle_outward (e02, e23, e13) p_in], 2)
     case 9  ->
       let p_in = v3_lerp p0 p3 0.5
       in ([orient_triangle_outward (e01, e13, e02) p_in, orient_triangle_outward (e02, e13, e23) p_in], 2)
     case _  -> ([dummy_tri, dummy_tri], 0)

-- ----------------------------------------------------------------------------
-- Marching Tetrahedra on a Regular 3D Voxel Grid
-- ----------------------------------------------------------------------------

def extract_isosurface [nx] [ny] [nz]
    (field: [nx][ny][nz]f64)
    (origin: vec3)
    (spacing: vec3)
    (isoval: f64): ([]triangle) =
  let n_cells = (nx - 1) * (ny - 1) * (nz - 1)
  let cell_candidates =
    tabulate n_cells (\cell_idx ->
      let iz = cell_idx % (nz - 1)
      let rem = cell_idx / (nz - 1)
      let iy = rem % (ny - 1)
      let ix = rem / (ny - 1)
      let get_pt (dx: i64) (dy: i64) (dz: i64): (vec3, f64) =
        let x = origin.0 + f64.i64 (ix + dx) * spacing.0
        let y = origin.1 + f64.i64 (iy + dy) * spacing.1
        let z = origin.2 + f64.i64 (iz + dz) * spacing.2
        in ((x, y, z), field[ix + dx, iy + dy, iz + dz])
      let (p0, v0) = get_pt 0 0 0
      let (p1, v1) = get_pt 1 0 0
      let (p2, v2) = get_pt 1 1 0
      let (p3, v3) = get_pt 0 1 0
      let (p4, v4) = get_pt 0 0 1
      let (p5, v5) = get_pt 1 0 1
      let (p6, v6) = get_pt 1 1 1
      let (p7, v7) = get_pt 0 1 1
      -- 5 tetrahedra decompose a cube
      let (t1, c1) = polygonize_tetrahedron isoval p0 v0 p1 v1 p3 v3 p4 v4
      let (t2, c2) = polygonize_tetrahedron isoval p1 v1 p2 v2 p3 v3 p6 v6
      let (t3, c3) = polygonize_tetrahedron isoval p1 v1 p4 v4 p5 v5 p6 v6
      let (t4, c4) = polygonize_tetrahedron isoval p3 v3 p4 v4 p6 v6 p7 v7
      let (t5, c5) = polygonize_tetrahedron isoval p1 v1 p3 v3 p4 v4 p6 v6
      in [ (t1[0], c1 >= 1), (t1[1], c1 >= 2)
         , (t2[0], c2 >= 1), (t2[1], c2 >= 2)
         , (t3[0], c3 >= 1), (t3[1], c3 >= 2)
         , (t4[0], c4 >= 1), (t4[1], c4 >= 2)
         , (t5[0], c5 >= 1), (t5[1], c5 >= 2)
         ]
    )
  let flat = flatten cell_candidates
  in map (\(t, _) -> t) (filter (\(_, valid) -> valid) flat)

-- ----------------------------------------------------------------------------
-- Comprehensive Mesh Statistics Computation
-- ----------------------------------------------------------------------------

def compute_mesh_stats [m] (triangles: [m]triangle): mesh_stats =
  if m == 0 then
    { num_triangles = 0
    , surface_area = 0.0
    , volume = 0.0
    , bbox_min = (0.0, 0.0, 0.0)
    , bbox_max = (0.0, 0.0, 0.0)
    }
  else
    let areas = map triangle_area triangles
    let total_area = reduce (+) 0.0 areas
    let vols = map triangle_signed_volume triangles
    let total_vol = f64.abs (reduce (+) 0.0 vols)
    let init_box = (triangles[0].0, triangles[0].0)
    let (bmin, bmax) =
      reduce (\(min_a, max_a) (min_b, max_b) ->
        (v3_min min_a min_b, v3_max max_a max_b)
      ) init_box
      (map (\(v0, v1, v2) ->
        let tmin = v3_min v0 (v3_min v1 v2)
        let tmax = v3_max v0 (v3_max v1 v2)
        in (tmin, tmax)
      ) triangles)
    in { num_triangles = m
       , surface_area = total_area
       , volume = total_vol
       , bbox_min = bmin
       , bbox_max = bmax
       }

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_mesh_stats [m]
    (v0x: [m]f64, v0y: [m]f64, v0z: [m]f64)
    (v1x: [m]f64, v1y: [m]f64, v1z: [m]f64)
    (v2x: [m]f64, v2y: [m]f64, v2z: [m]f64): (i64, f64, f64, f64, f64, f64, f64, f64, f64) =
  let tris =
    tabulate m (\i ->
      ((v0x[i], v0y[i], v0z[i]), (v1x[i], v1y[i], v1z[i]), (v2x[i], v2y[i], v2z[i]))
    )
  let stats = compute_mesh_stats tris
  in ( stats.num_triangles
     , stats.surface_area
     , stats.volume
     , stats.bbox_min.0, stats.bbox_min.1, stats.bbox_min.2
     , stats.bbox_max.0, stats.bbox_max.1, stats.bbox_max.2
     )
