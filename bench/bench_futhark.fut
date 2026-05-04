-- Caliper Metrology Acceleration Microbenchmarks in Futhark

import "../src/futhark/math"
import "../src/futhark/mesh"
import "../src/futhark/census"
import "../src/futhark/distance"

entry bench_census [h][w] (left: [h][w]f64) (right: [h][w]f64): [h][w][]i32 =
  let max_disp = 16i64
  let l_cen = compute_census_map left
  let r_cen = compute_census_map right
  in compute_cost_volume l_cen r_cen max_disp

entry bench_distance [n] (pts: [n]vec3) (tri_a: vec3) (tri_b: vec3) (tri_c: vec3): [n]f64 =
  map (\p ->
    let (_, dist) = point_triangle_closest p (tri_a, tri_b, tri_c)
    in dist
  ) pts
