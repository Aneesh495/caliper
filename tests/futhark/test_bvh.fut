-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Bounding Volume Hierarchy & Morton Codes
-- Module: tests/futhark/test_bvh.fut
--
-- Tests AABB metrics, point-to-box Euclidean distances, and 3D Morton code
-- space-filling curve ordering.
-- ============================================================================

import "../../src/futhark/bvh"
import "../../src/futhark/math"
import "../../src/futhark/mesh"

def test_aabb_metrics: bool =
  let b1 = aabb_make (0.0, 0.0, 0.0) (2.0, 2.0, 2.0)
  let b2 = aabb_make (1.0, 1.0, 1.0) (3.0, 4.0, 5.0)
  let u = aabb_union b1 b2
  let ok_min = f64.abs u.min_pt.0 < 1e-9 && f64.abs u.min_pt.1 < 1e-9 && f64.abs u.min_pt.2 < 1e-9
  let ok_max = f64.abs (u.max_pt.0 - 3.0) < 1e-9 && f64.abs (u.max_pt.1 - 4.0) < 1e-9 && f64.abs (u.max_pt.2 - 5.0) < 1e-9
  let p_inside: vec3 = (1.5, 1.5, 1.5)
  let p_outside: vec3 = (5.0, 2.0, 2.0)
  let d_inside = point_aabb_dist p_inside b1
  let d_outside = point_aabb_dist p_outside b1 -- Distance along X is 5.0 - 2.0 = 3.0
  let ok_dist_in = f64.abs d_inside < 1e-9
  let ok_dist_out = f64.abs (d_outside - 3.0) < 1e-9
  in ok_min && ok_max && ok_dist_in && ok_dist_out

def test_morton_encoding: bool =
  -- Points along the diagonal (0,0,0) -> (1,1,1)
  let p0: vec3 = (0.0, 0.0, 0.0)
  let p1: vec3 = (0.5, 0.5, 0.5)
  let p2: vec3 = (1.0, 1.0, 1.0)
  let k0 = morton_code_3d p0
  let k1 = morton_code_3d p1
  let k2 = morton_code_3d p2
  let ok_monotonic = k0 < k1 && k1 < k2
  let ok_k0 = k0 == 0u32
  in ok_monotonic && ok_k0

entry main: bool =
  test_aabb_metrics && test_morton_encoding
