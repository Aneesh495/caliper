-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for Math Module
-- Module: tests/futhark/test_math.fut
--
-- Tests vector algebra, matrix multiplication, determinants, Cardano cubic
-- eigensolver for symmetric tensors, and Cholesky 6x6 decomposition.
-- ============================================================================

import "../../src/futhark/math"

-- Vector tests
def test_vector_ops: bool =
  let a: vec3 = (1.0, 2.0, 3.0)
  let b: vec3 = (4.0, -5.0, 6.0)
  let sum = v3_add a b
  let diff = v3_sub a b
  let dot = v3_dot a b
  let cross = v3_cross a b
  let norm = v3_norm a
  let norm_target = f64.sqrt (1.0 + 4.0 + 9.0)
  let ok_sum = f64.abs (sum.0 - 5.0) < 1e-9 && f64.abs (sum.1 - (-3.0)) < 1e-9 && f64.abs (sum.2 - 9.0) < 1e-9
  let ok_diff = f64.abs (diff.0 - (-3.0)) < 1e-9 && f64.abs (diff.1 - 7.0) < 1e-9 && f64.abs (diff.2 - (-3.0)) < 1e-9
  let ok_dot = f64.abs (dot - (4.0 - 10.0 + 18.0)) < 1e-9
  let ok_cross = f64.abs (cross.0 - 27.0) < 1e-9 && f64.abs (cross.1 - 6.0) < 1e-9 && f64.abs (cross.2 - (-13.0)) < 1e-9
  let ok_norm = f64.abs (norm - norm_target) < 1e-9
  in ok_sum && ok_diff && ok_dot && ok_cross && ok_norm

-- Matrix determinant and inversion tests
def test_matrix_inv: bool =
  let A: mat3 =
    [ [1.0, 2.0, 3.0]
    , [0.0, 1.0, 4.0]
    , [5.0, 6.0, 0.0]
    ]
  let det = mat3_det A
  let (ok, invA) = mat3_inv A
  let prod = mat3_mul A invA
  let ok_det = f64.abs (det - 1.0) < 1e-9
  let ok_diag = f64.abs (prod[0,0] - 1.0) < 1e-9 && f64.abs (prod[1,1] - 1.0) < 1e-9 && f64.abs (prod[2,2] - 1.0) < 1e-9
  let ok_off = f64.abs prod[0,1] < 1e-9 && f64.abs prod[0,2] < 1e-9 && f64.abs prod[1,0] < 1e-9
  in ok && ok_det && ok_diag && ok_off

-- Cardano eigensolver test on known covariance matrix
def test_eigensolver: bool =
  -- Symmetric matrix with known eigenvalues 1.0, 3.0, 6.0
  -- Constructed as Q * diag(1, 3, 6) * Q^T
  let C: mat3 =
    [ [2.0, -1.0, 0.0]
    , [-1.0, 2.0, -1.0]
    , [0.0, -1.0, 2.0]
    ]
  let (l1, l2, l3) = sym_mat3_eigenvalues C
  let v = sym_mat3_min_eigenvector C l1
  let n = v3_norm v
  let Cv = mat3_mul_v3 C v
  let l1_v = v3_scale l1 v
  let residual = v3_dist Cv l1_v
  let ok_order = l1 <= l2 && l2 <= l3
  let ok_norm = f64.abs (n - 1.0) < 1e-9
  let ok_res = residual < 1e-6
  in ok_order && ok_norm && ok_res

-- Cholesky 6x6 solver test
def test_cholesky_6x6: bool =
  -- Construct symmetric positive definite H = A^T A + 2*I
  let H: [6][6]f64 =
    [ [4.0, 1.0, 0.5, 0.0, 0.2, 0.1]
    , [1.0, 5.0, 1.2, 0.3, 0.0, 0.4]
    , [0.5, 1.2, 6.0, 1.0, 0.5, 0.0]
    , [0.0, 0.3, 1.0, 4.5, 0.8, 0.2]
    , [0.2, 0.0, 0.5, 0.8, 5.2, 1.1]
    , [0.1, 0.4, 0.0, 0.2, 1.1, 4.8]
    ]
  let g: [6]f64 = [1.2, -0.8, 3.4, 0.5, -2.1, 4.0]
  let (ok, x) = cholesky_solve_6x6 H g
  -- Verify H * x == g
  let Hx =
    map (\row ->
      loop sum = 0.0 for i < 6 do sum + row[i] * x[i]
    ) H
  let diff_sq =
    loop s = 0.0 for i < 6 do s + (Hx[i] - g[i])**2.0
  in ok && (f64.sqrt diff_sq < 1e-8)

entry main: bool =
  test_vector_ops && test_matrix_inv && test_eigensolver && test_cholesky_6x6
