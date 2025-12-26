-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Mesh Extraction
-- Module: tests/futhark/test_mesh.fut
--
-- Tests Marching Tetrahedra isosurface polygonization on implicit 3D spheres
-- and verifies closed surface volume and surface area metrics.
-- ============================================================================

import "../../src/futhark/mesh"
import "../../src/futhark/math"

def test_single_triangle_metrics: bool =
  -- Right isosceles triangle in XY plane with vertices (0,0,0), (2,0,0), (0,2,0)
  let tri: triangle = ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0))
  let area = triangle_area tri
  let norm = triangle_normal tri
  let ok_area = f64.abs (area - 2.0) < 1e-9 -- Area = 0.5 * 2 * 2 = 2.0
  let ok_norm = f64.abs norm.0 < 1e-9 && f64.abs norm.1 < 1e-9 && f64.abs (norm.2 - 1.0) < 1e-9
  in ok_area && ok_norm

def test_sphere_isosurface_extraction: bool =
  -- Sample implicit sphere field on 12x12x12 grid: f(x, y, z) = x^2 + y^2 + z^2 - R^2
  let n = 12i64
  let radius = 2.0
  let origin: vec3 = (-3.0, -3.0, -3.0)
  let spacing: vec3 = (0.5, 0.5, 0.5)
  let field =
    tabulate_3d n n n (\ix iy iz ->
      let x = origin.0 + f64.i64 ix * spacing.0
      let y = origin.1 + f64.i64 iy * spacing.1
      let z = origin.2 + f64.i64 iz * spacing.2
      in x*x + y*y + z*z - radius*radius
    )
  let isoval = 0.0
  let triangles = extract_isosurface field origin spacing isoval
  let stats = compute_mesh_stats triangles
  -- Theoretical sphere volume: (4/3) * pi * R^3 = 4/3 * 3.14159 * 8 = 33.51
  -- Theoretical surface area: 4 * pi * R^2 = 50.26
  let expected_vol = (4.0 / 3.0) * f64.pi * (radius**3.0)
  let expected_area = 4.0 * f64.pi * (radius**2.0)
  let ok_count = stats.num_triangles > 50
  let ok_vol = f64.abs (stats.volume - expected_vol) / expected_vol < 0.35
  let ok_area = f64.abs (stats.surface_area - expected_area) / expected_area < 0.10
  in ok_count && ok_vol && ok_area

entry main: bool =
  test_single_triangle_metrics && test_sphere_isosurface_extraction
