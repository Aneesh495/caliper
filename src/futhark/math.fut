-- ============================================================================
-- Caliper Metrology Engine: High-Performance Geometric & Linear Algebra
-- Module: math.fut
--
-- Pure data-parallel linear algebra, 3D Euclidean geometry, analytic Cardano
-- eigensolvers for symmetric covariance tensors, and Cholesky 6x6 solvers
-- for linearized point-to-plane ICP normal equations.
-- ============================================================================

-- ----------------------------------------------------------------------------
-- Vector 3D Primitives and Vectorized Kernels
-- ----------------------------------------------------------------------------

type vec3 = (f64, f64, f64)
type vec4 = (f64, f64, f64, f64)
type mat3 = [3][3]f64
type mat4 = [4][4]f64
type quat = (f64, f64, f64, f64) -- (w, x, y, z)

def v3_make (x: f64) (y: f64) (z: f64): vec3 = (x, y, z)

def v3_to_array ((x, y, z): vec3): [3]f64 = [x, y, z]

def v3_from_array (arr: [3]f64): vec3 = (arr[0], arr[1], arr[2])

def v3_add ((x1, y1, z1): vec3) ((x2, y2, z2): vec3): vec3 =
  (x1 + x2, y1 + y2, z1 + z2)

def v3_sub ((x1, y1, z1): vec3) ((x2, y2, z2): vec3): vec3 =
  (x1 - x2, y1 - y2, z1 - z2)

def v3_scale (s: f64) ((x, y, z): vec3): vec3 =
  (s * x, s * y, s * z)

def v3_dot ((x1, y1, z1): vec3) ((x2, y2, z2): vec3): f64 =
  x1 * x2 + y1 * y2 + z1 * z2

def v3_cross ((x1, y1, z1): vec3) ((x2, y2, z2): vec3): vec3 =
  ( y1 * z2 - z1 * y2
  , z1 * x2 - x1 * z2
  , x1 * y2 - y1 * x2
  )

def v3_norm_sq (v: vec3): f64 =
  v3_dot v v

def v3_norm (v: vec3): f64 =
  f64.sqrt (v3_norm_sq v)

def v3_normalize (v: vec3): vec3 =
  let n = v3_norm v
  in if n > 1e-12
     then v3_scale (1.0 / n) v
     else (0.0, 0.0, 0.0)

def v3_dist (a: vec3) (b: vec3): f64 =
  v3_norm (v3_sub a b)

def v3_dist_sq (a: vec3) (b: vec3): f64 =
  v3_norm_sq (v3_sub a b)

def v3_lerp (a: vec3) (b: vec3) (t: f64): vec3 =
  v3_add (v3_scale (1.0 - t) a) (v3_scale t b)

def v3_min ((x1, y1, z1): vec3) ((x2, y2, z2): vec3): vec3 =
  (f64.min x1 x2, f64.min y1 y2, f64.min z1 z2)

def v3_max ((x1, y1, z1): vec3) ((x2, y2, z2): vec3): vec3 =
  (f64.max x1 x2, f64.max y1 y2, f64.max z1 z2)

-- ----------------------------------------------------------------------------
-- Matrix 3x3 Primitives and Vectorized Algebra
-- ----------------------------------------------------------------------------

def mat3_identity: mat3 =
  [ [1.0, 0.0, 0.0]
  , [0.0, 1.0, 0.0]
  , [0.0, 0.0, 1.0]
  ]

def mat3_zeros: mat3 =
  [ [0.0, 0.0, 0.0]
  , [0.0, 0.0, 0.0]
  , [0.0, 0.0, 0.0]
  ]

def mat3_add (A: mat3) (B: mat3): mat3 =
  map2 (\rowA rowB -> map2 (+) rowA rowB) A B

def mat3_scale (s: f64) (A: mat3): mat3 =
  map (\row -> map (\elem -> s * elem) row) A

def mat3_transpose (A: mat3): mat3 =
  [ [A[0,0], A[1,0], A[2,0]]
  , [A[0,1], A[1,1], A[2,1]]
  , [A[0,2], A[1,2], A[2,2]]
  ]

def mat3_mul_v3 (A: mat3) ((x, y, z): vec3): vec3 =
  ( A[0,0]*x + A[0,1]*y + A[0,2]*z
  , A[1,0]*x + A[1,1]*y + A[1,2]*z
  , A[2,0]*x + A[2,1]*y + A[2,2]*z
  )

def mat3_mul (A: mat3) (B: mat3): mat3 =
  let BT = mat3_transpose B
  in map (\rowA ->
       map (\colB ->
         rowA[0]*colB[0] + rowA[1]*colB[1] + rowA[2]*colB[2]
       ) BT
     ) A

def mat3_det (A: mat3): f64 =
  A[0,0] * (A[1,1]*A[2,2] - A[1,2]*A[2,1])
  - A[0,1] * (A[1,0]*A[2,2] - A[1,2]*A[2,0])
  + A[0,2] * (A[1,0]*A[2,1] - A[1,1]*A[2,0])

def mat3_trace (A: mat3): f64 =
  A[0,0] + A[1,1] + A[2,2]

def mat3_inv (A: mat3): (bool, mat3) =
  let det = mat3_det A
  in if f64.abs det < 1e-15
     then (false, mat3_zeros)
     else
       let inv_det = 1.0 / det
       let c00 =  (A[1,1]*A[2,2] - A[1,2]*A[2,1]) * inv_det
       let c01 = -(A[0,1]*A[2,2] - A[0,2]*A[2,1]) * inv_det
       let c02 =  (A[0,1]*A[1,2] - A[0,2]*A[1,1]) * inv_det
       let c10 = -(A[1,0]*A[2,2] - A[1,2]*A[2,0]) * inv_det
       let c11 =  (A[0,0]*A[2,2] - A[0,2]*A[2,0]) * inv_det
       let c12 = -(A[0,0]*A[1,2] - A[0,2]*A[1,0]) * inv_det
       let c20 =  (A[1,0]*A[2,1] - A[1,1]*A[2,0]) * inv_det
       let c21 = -(A[0,0]*A[2,1] - A[0,1]*A[2,0]) * inv_det
       let c22 =  (A[0,0]*A[1,1] - A[0,1]*A[1,0]) * inv_det
       in (true, [ [c00, c01, c02], [c10, c11, c12], [c20, c21, c22] ])

-- ----------------------------------------------------------------------------
-- Analytic Cardano Eigensolver for 3x3 Real Symmetric Matrices
--
-- Computes exact closed-form eigenvalues lambda1 <= lambda2 <= lambda3
-- of a real symmetric 3x3 matrix without iterative Jacobi sweeps.
-- Highly efficient for point cloud normal and curvature estimation.
-- ----------------------------------------------------------------------------

def sym_mat3_eigenvalues (A: mat3): (f64, f64, f64) =
  let a00 = A[0,0]
  let a11 = A[1,1]
  let a22 = A[2,2]
  let a01 = A[0,1]
  let a02 = A[0,2]
  let a12 = A[1,2]
  let p1 = a01*a01 + a02*a02 + a12*a12
  in if p1 < 1e-18
     then
       -- Diagonal matrix case: sort diagonal elements directly
       let vals = [a00, a11, a22]
       let l1 = f64.min vals[0] (f64.min vals[1] vals[2])
       let l3 = f64.max vals[0] (f64.max vals[1] vals[2])
       let l2 = vals[0] + vals[1] + vals[2] - l1 - l3
       in (l1, l2, l3)
     else
       let q = (a00 + a11 + a22) / 3.0
       let p2 = (a00 - q)**2.0 + (a11 - q)**2.0 + (a22 - q)**2.0 + 2.0 * p1
       let p = f64.sqrt (p2 / 6.0)
       let inv_p = 1.0 / p
       let B = [ [(a00 - q)*inv_p, a01*inv_p, a02*inv_p]
               , [a01*inv_p, (a11 - q)*inv_p, a12*inv_p]
               , [a02*inv_p, a12*inv_p, (a22 - q)*inv_p]
               ]
       let r = mat3_det B / 2.0
       let r_clamped = f64.max (-1.0) (f64.min 1.0 r)
       let phi = (f64.acos r_clamped) / 3.0
       let pi = f64.pi
       let eig1 = q + 2.0 * p * f64.cos (phi + (2.0 * pi / 3.0))
       let eig3 = q + 2.0 * p * f64.cos phi
       let eig2 = 3.0 * q - eig1 - eig3
       in (eig1, eig2, eig3)

-- ----------------------------------------------------------------------------
-- Eigenvector for Minimum Eigenvalue of Symmetric 3x3 Matrix
--
-- Uses cross-product formulation to extract the principal normal vector
-- corresponding to lambda1 (minimum eigenvalue).
-- ----------------------------------------------------------------------------

def sym_mat3_min_eigenvector (A: mat3) (l1: f64): vec3 =
  let a00 = A[0,0] - l1
  let a11 = A[1,1] - l1
  let a22 = A[2,2] - l1
  let a01 = A[0,1]
  let a02 = A[0,2]
  let a12 = A[1,2]
  -- Rows of (A - l1 * I)
  let r0: vec3 = (a00, a01, a02)
  let r1: vec3 = (a01, a11, a12)
  let r2: vec3 = (a02, a12, a22)
  -- Compute cross products between pairs of rows
  let c01 = v3_cross r0 r1
  let c02 = v3_cross r0 r2
  let c12 = v3_cross r1 r2
  let n01 = v3_norm_sq c01
  let n02 = v3_norm_sq c02
  let n12 = v3_norm_sq c12
  let best =
    if n01 >= n02 && n01 >= n12 then c01
    else if n02 >= n01 && n02 >= n12 then c02
    else c12
  in v3_normalize best

-- ----------------------------------------------------------------------------
-- Quaternions and Rigid Transformations SE(3)
-- ----------------------------------------------------------------------------

def quat_identity: quat = (1.0, 0.0, 0.0, 0.0)

def quat_normalize ((w, x, y, z): quat): quat =
  let n = f64.sqrt (w*w + x*x + y*y + z*z)
  in if n > 1e-12
     then (w / n, x / n, y / n, z / n)
     else (1.0, 0.0, 0.0, 0.0)

def quat_to_mat3 ((w, x, y, z): quat): mat3 =
  let (w_n, x_n, y_n, z_n) = quat_normalize (w, x, y, z)
  let xx = x_n * x_n
  let yy = y_n * y_n
  let zz = z_n * z_n
  let xy = x_n * y_n
  let xz = x_n * z_n
  let yz = y_n * z_n
  let wx = w_n * x_n
  let wy = w_n * y_n
  let wz = w_n * z_n
  in [ [1.0 - 2.0*(yy + zz),       2.0*(xy - wz),       2.0*(xz + wy)]
     , [      2.0*(xy + wz), 1.0 - 2.0*(xx + zz),       2.0*(yz - wx)]
     , [      2.0*(xz - wy),       2.0*(yz + wx), 1.0 - 2.0*(xx + yy)]
     ]

def transform_point (R: mat3) (t: vec3) (p: vec3): vec3 =
  v3_add (mat3_mul_v3 R p) t

def transform_normal (R: mat3) (n: vec3): vec3 =
  v3_normalize (mat3_mul_v3 R n)

-- ----------------------------------------------------------------------------
-- 6x6 Cholesky Linear System Solver for Point-to-Plane ICP
--
-- Solves (A^T A + lambda * I) x = A^T b where dimension is exactly 6x6.
-- Uses numerical Cholesky decomposition L * L^T = H with diagonal regularization.
-- ----------------------------------------------------------------------------

def cholesky_solve_6x6 (H: [6][6]f64) (g: [6]f64): (bool, [6]f64) =
  let (L, ok) =
    loop (L, ok) = (replicate 6 (replicate 6 0.0), true) for i < 6 do
      if !ok then (L, false)
      else
        loop (L_sub, ok_sub) = (L, ok) for j < i + 1 do
          if !ok_sub then (L_sub, false)
          else
            let s_prev = loop s = 0.0 for k < j do s + L_sub[i, k] * L_sub[j, k]
            let residual = H[i, j] - s_prev
            in if i == j then
                 if residual > 1e-12 then (L_sub with [i, j] = f64.sqrt residual, true)
                 else (L_sub, false)
               else
                 let diag = L_sub[j, j]
                 in if f64.abs diag > 1e-12 then (L_sub with [i, j] = residual / diag, true)
                 else (L_sub, false)
  in if !ok then (false, replicate 6 0.0)
     else
       -- Forward substitution: L * y = g
       let y =
         loop y = replicate 6 0.0 for i < 6 do
           let s = loop sum = 0.0 for k < i do sum + L[i, k] * y[k]
           in y with [i] = (g[i] - s) / L[i, i]
       -- Backward substitution: L^T * x = y
       let x =
         loop x = replicate 6 0.0 for step < 6 do
           let i = 5 - step
           let s = loop sum = 0.0 for step_k < step do
                     let k = 5 - step_k
                     in sum + L[k, i] * x[k]
           in x with [i] = (y[i] - s) / L[i, i]
       in (true, x)

-- ----------------------------------------------------------------------------
-- Entry points for Futhark compiler exports
-- ----------------------------------------------------------------------------

entry fut_v3_dot (x1: f64, y1: f64, z1: f64, x2: f64, y2: f64, z2: f64): f64 =
  v3_dot (x1, y1, z1) (x2, y2, z2)

entry fut_v3_cross (x1: f64, y1: f64, z1: f64, x2: f64, y2: f64, z2: f64): (f64, f64, f64) =
  v3_cross (x1, y1, z1) (x2, y2, z2)

entry fut_mat3_mul (A: [3][3]f64, B: [3][3]f64): [3][3]f64 =
  mat3_mul A B

entry fut_sym_eigen (A: [3][3]f64): (f64, f64, f64, f64, f64, f64) =
  let (l1, l2, l3) = sym_mat3_eigenvalues A
  let (nx, ny, nz) = sym_mat3_min_eigenvector A l1
  in (l1, l2, l3, nx, ny, nz)

entry fut_cholesky_6x6 (H: [6][6]f64, g: [6]f64): (bool, [6]f64) =
  cholesky_solve_6x6 H g
