-- ============================================================================
-- Caliper Metrology Engine: Linear Bounding Volume Hierarchy (LBVH)
-- Module: bvh.fut
--
-- Morton code spatial key generation (Z-order curve), Axis-Aligned Bounding Box
-- (AABB) hierarchy, and data-parallel spatial bounding queries for fast SDF evaluation.
-- ============================================================================

import "math"
import "mesh"

type aabb =
  { min_pt: vec3
  , max_pt: vec3
  }

-- ----------------------------------------------------------------------------
-- AABB Primitives & Distance Metrics
-- ----------------------------------------------------------------------------

def aabb_make (min_p: vec3) (max_p: vec3): aabb =
  { min_pt = min_p, max_pt = max_p }

def aabb_empty: aabb =
  { min_pt = ( 1e30,  1e30,  1e30)
  , max_pt = (-1e30, -1e30, -1e30)
  }

def aabb_from_point (p: vec3): aabb =
  { min_pt = p, max_pt = p }

def aabb_from_triangle ((v0, v1, v2): triangle): aabb =
  { min_pt = v3_min v0 (v3_min v1 v2)
  , max_pt = v3_max v0 (v3_max v1 v2)
  }

def aabb_union (b1: aabb) (b2: aabb): aabb =
  { min_pt = v3_min b1.min_pt b2.min_pt
  , max_pt = v3_max b1.max_pt b2.max_pt
  }

def aabb_center (box: aabb): vec3 =
  v3_scale 0.5 (v3_add box.min_pt box.max_pt)

def aabb_extents (box: aabb): vec3 =
  v3_sub box.max_pt box.min_pt

-- ----------------------------------------------------------------------------
-- Minimum Euclidean Distance from Point to AABB
--
-- Used to aggressively prune distant BVH subtrees during nearest triangle queries.
-- ----------------------------------------------------------------------------

def point_aabb_dist_sq (p: vec3) (box: aabb): f64 =
  let dx = f64.max 0.0 (f64.max (box.min_pt.0 - p.0) (p.0 - box.max_pt.0))
  let dy = f64.max 0.0 (f64.max (box.min_pt.1 - p.1) (p.1 - box.max_pt.1))
  let dz = f64.max 0.0 (f64.max (box.min_pt.2 - p.2) (p.2 - box.max_pt.2))
  in dx*dx + dy*dy + dz*dz

def point_aabb_dist (p: vec3) (box: aabb): f64 =
  f64.sqrt (point_aabb_dist_sq p box)

def point_inside_aabb (p: vec3) (box: aabb): bool =
  p.0 >= box.min_pt.0 && p.0 <= box.max_pt.0 &&
  p.1 >= box.min_pt.1 && p.1 <= box.max_pt.1 &&
  p.2 >= box.min_pt.2 && p.2 <= box.max_pt.2

-- ----------------------------------------------------------------------------
-- 3D Morton Code (Z-Order Curve) Bit Interleaving
--
-- Maps 3D coordinates normalized in [0, 1] to a 30-bit spatial locality key.
-- Preserves spatial clusters along recursive octree space partitioning.
-- ----------------------------------------------------------------------------

def expand_bits_10 (v: u32): u32 =
  let v1 = (v * 0x00010001u32) & 0xFF0000FFu32
  let v2 = (v1 ^ (v1 << 8u32)) & 0x0F00F00Fu32
  let v3 = (v2 ^ (v2 << 4u32)) & 0xC30C30C3u32
  let v4 = (v3 ^ (v3 << 2u32)) & 0x49249249u32
  in v4

def morton_code_3d (p_norm: vec3): u32 =
  -- Clamp normalized coordinates into [0.0, 1.0]
  let cx = f64.max 0.0 (f64.min 1.0 p_norm.0)
  let cy = f64.max 0.0 (f64.min 1.0 p_norm.1)
  let cz = f64.max 0.0 (f64.min 1.0 p_norm.2)
  -- Discretize to 10-bit integer [0, 1023]
  let ix = u32.i64 (i64.f64 (cx * 1023.0))
  let iy = u32.i64 (i64.f64 (cy * 1023.0))
  let iz = u32.i64 (i64.f64 (cz * 1023.0))
  let ex = expand_bits_10 ix
  let ey = expand_bits_10 iy
  let ez = expand_bits_10 iz
  in (ex << 2u32) | (ey << 1u32) | ez

-- ----------------------------------------------------------------------------
-- Parallel Morton Key Generation for Mesh Triangles
-- ----------------------------------------------------------------------------

def compute_triangle_morton_keys [m] (triangles: [m]triangle): ([m]u32, aabb) =
  let boxes = map aabb_from_triangle triangles
  let global_box = reduce aabb_union aabb_empty boxes
  let extents = aabb_extents global_box
  let safe_ext =
    ( f64.max 1e-9 extents.0
    , f64.max 1e-9 extents.1
    , f64.max 1e-9 extents.2
    )
  let keys =
    map (\tri ->
      let c = aabb_center (aabb_from_triangle tri)
      let norm_c: vec3 =
        ( (c.0 - global_box.min_pt.0) / safe_ext.0
        , (c.1 - global_box.min_pt.1) / safe_ext.1
        , (c.2 - global_box.min_pt.2) / safe_ext.2
        )
      in morton_code_3d norm_c
    ) triangles
  in (keys, global_box)

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_morton_keys [m]
    (v0x: [m]f64, v0y: [m]f64, v0z: [m]f64)
    (v1x: [m]f64, v1y: [m]f64, v1z: [m]f64)
    (v2x: [m]f64, v2y: [m]f64, v2z: [m]f64): [m]u32 =
  let tris =
    tabulate m (\i ->
      ((v0x[i], v0y[i], v0z[i]), (v1x[i], v1y[i], v1z[i]), (v2x[i], v2y[i], v2z[i]))
    )
  let (keys, _) = compute_triangle_morton_keys tris
  in keys

entry fut_point_aabb_dist
    (px: f64, py: f64, pz: f64)
    (bmin_x: f64, bmin_y: f64, bmin_z: f64)
    (bmax_x: f64, bmax_y: f64, bmax_z: f64): f64 =
  let p = (px, py, pz)
  let box = { min_pt = (bmin_x, bmin_y, bmin_z), max_pt = (bmax_x, bmax_y, bmax_z) }
  in point_aabb_dist p box
