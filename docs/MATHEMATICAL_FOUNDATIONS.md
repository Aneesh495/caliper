# Mathematical Foundations of Caliper

This document provides formal derivations, algorithmic proofs, and numerical formulations underpinning the Caliper metrology engine.

---

## 1. Optical Lens Distortion & Contraction Inversion

### 1.1 Forward Brown-Conrady Model
Let $\mathbf{x}_u = (x_u, y_u)^T$ be normalized undistorted image coordinates and $r^2 = x_u^2 + y_u^2$. The forward distorted coordinates $\mathbf{x}_d = (x_d, y_d)^T$ are given by:

$$\mathbf{x}_d = \mathbf{x}_u (1 + k_1 r^2 + k_2 r^4 + k_3 r^6) + \mathbf{d}_{\text{tan}}(\mathbf{x}_u)$$

where tangential distortion $\mathbf{d}_{\text{tan}}$ is parameterized by coefficients $(p_1, p_2)$:

$$\mathbf{d}_{\text{tan}}(\mathbf{x}_u) = \begin{bmatrix} 2 p_1 x_u y_u + p_2 (r^2 + 2 x_u^2) \\ p_1 (r^2 + 2 y_u^2) + 2 p_2 x_u y_u \end{bmatrix}$$

### 1.2 Inverse Mapping via Newton-Raphson Contraction
Given measured distorted pixel coordinates $\mathbf{x}_d$, finding undistorted coordinates $\mathbf{x}_u$ requires solving the non-linear root equation:

$$\mathbf{F}(\mathbf{x}_u) = \mathbf{x}_d - \mathbf{g}(\mathbf{x}_u) = \mathbf{0}$$

where $\mathbf{g}(\mathbf{x}_u) = \mathbf{x}_u (1 + k_1 r^2 + k_2 r^4 + k_3 r^6) + \mathbf{d}_{\text{tan}}(\mathbf{x}_u)$.

The Jacobian matrix $J(\mathbf{x}_u) = \frac{\partial \mathbf{g}}{\partial \mathbf{x}_u}$ evaluates to:

$$J_{00} = 1 + k_1(3 x_u^2 + y_u^2) + k_2(5 x_u^4 + 6 x_u^2 y_u^2 + y_u^4) + 2 p_1 y_u + 6 p_2 x_u$$
$$J_{01} = 2 k_1 x_u y_u + 4 k_2 x_u y_u r^2 + 2 p_1 x_u + 2 p_2 y_u$$
$$J_{10} = 2 k_1 x_u y_u + 4 k_2 x_u y_u r^2 + 2 p_1 x_u + 2 p_2 y_u$$
$$J_{11} = 1 + k_1(x_u^2 + 3 y_u^2) + k_2(x_u^4 + 6 x_u^2 y_u^2 + 5 y_u^4) + 6 p_1 y_u + 2 p_2 x_u$$

The iterative contraction update follows:

$$\mathbf{x}_u^{(k+1)} = \mathbf{x}_u^{(k)} - J\left(\mathbf{x}_u^{(k)}\right)^{-1} \left( \mathbf{g}\left(\mathbf{x}_u^{(k)}\right) - \mathbf{x}_d \right)$$

On modern industrial sensors, $\|J - I\| < 0.15$, guaranteeing superlinear quadratic convergence within 5 iterations.

---

## 2. Semi-Global Matching Dynamic Programming

Stereo depth estimation minimizes global energy functional $E(D)$:

$$E(D) = \sum_{\mathbf{p}} C(\mathbf{p}, D(\mathbf{p})) + \sum_{\mathbf{q} \in N_{\mathbf{p}}} P_1 \cdot \mathbf{1}_{|D(\mathbf{p}) - D(\mathbf{q})| = 1} + \sum_{\mathbf{q} \in N_{\mathbf{p}}} P_2 \cdot \mathbf{1}_{|D(\mathbf{p}) - D(\mathbf{q})| > 1}$$

Because 2D MRF energy minimization is NP-hard, Semi-Global Matching computes dynamic programming approximations along 1D paths $\mathbf{r} \in \{(0, 1), (0, -1), (1, 0), (-1, 0)\}$.

The path cost recurrence $L_{\mathbf{r}}(\mathbf{p}, d)$ at pixel $\mathbf{p}$ for candidate disparity $d$ is:

$$L_{\mathbf{r}}(\mathbf{p}, d) = C(\mathbf{p}, d) + \min \begin{cases} L_{\mathbf{r}}(\mathbf{p} - \mathbf{r}, d) \\ L_{\mathbf{r}}(\mathbf{p} - \mathbf{r}, d - 1) + P_1 \\ L_{\mathbf{r}}(\mathbf{p} - \mathbf{r}, d + 1) + P_1 \\ \min_k L_{\mathbf{r}}(\mathbf{p} - \mathbf{r}, k) + P_2 \end{cases} - \min_k L_{\mathbf{r}}(\mathbf{p} - \mathbf{r}, k)$$

Subtracting $\min_k L_{\mathbf{r}}(\mathbf{p} - \mathbf{r}, k)$ bounds intermediate accumulators within 16-bit integers, preventing numerical overflow during continuous streaming.

### 2.1 Sub-Pixel Parabolic Peak Interpolation
Let $C(d)$ denote aggregated cost at discrete integer disparity $d$. The minimum cost is located at $d^* = \arg\min_d C(d)$. Taylor expansion of $C(d)$ around $d^*$ yields:

$$\delta d = \frac{C(d^* - 1) - C(d^* + 1)}{2 \left( C(d^* - 1) - 2 C(d^*) + C(d^* + 1) \right)}$$

The refined sub-pixel disparity is $d_{\text{refined}} = d^* + \delta d$, with $\delta d \in [-0.5, 0.5]$.

---

## 3. Closed-Form Cardano Eigensolver for $3 \times 3$ Tensors

For surface normal estimation, the local point covariance tensor $C \in \mathbb{R}^{3 \times 3}$ is symmetric positive semi-definite. The characteristic polynomial is:

$$\det(C - \lambda I) = -\lambda^3 + c_2 \lambda^2 - c_1 \lambda + c_0 = 0$$

where:
$$c_2 = \text{tr}(C) = C_{00} + C_{11} + C_{22}$$
$$c_1 = (C_{00} C_{11} - C_{01}^2) + (C_{00} C_{22} - C_{02}^2) + (C_{11} C_{22} - C_{12}^2)$$
$$c_0 = \det(C) = C_{00}(C_{11} C_{22} - C_{12}^2) - C_{01}(C_{01} C_{22} - C_{12} C_{02}) + C_{02}(C_{01} C_{12} - C_{11} C_{02})$$

Applying Tschirnhaus transformation $\lambda = \mu + \frac{c_2}{3}$ eliminates the quadratic term:

$$\mu^3 - 3 p \mu - 2 q = 0$$

where:
$$p = \frac{c_2^2 - 3 c_1}{9}, \quad q = \frac{2 c_2^3 - 9 c_2 c_1 + 27 c_0}{54}$$

Because $C$ is real symmetric, all roots are real, ensuring discriminant $D = p^3 - q^2 \ge 0$. Using the trigonometric identity $\cos(3\theta) = 4\cos^3(\theta) - 3\cos(\theta)$:

$$\phi = \frac{1}{3} \arccos\left( \text{clamp}\left( \frac{q}{\sqrt{p^3}}, -1.0, 1.0 \right) \right)$$

The eigenvalues in descending order are:

$$\lambda_1 = \frac{c_2}{3} + 2 \sqrt{p} \cos(\phi)$$
$$\lambda_3 = \frac{c_2}{3} + 2 \sqrt{p} \cos\left(\phi + \frac{2\pi}{3}\right)$$
$$\lambda_2 = c_2 - \lambda_1 - \lambda_3$$

The minimal eigenvalue $\lambda_3$ represents surface variation along the normal direction. The eigenvector $\mathbf{n}$ is obtained by computing the cross product of two linearly independent rows of $(C - \lambda_3 I)$:

$$\mathbf{n} = (C - \lambda_3 I)_{0,:} \times (C - \lambda_3 I)_{1,:}$$

Normalizing $\mathbf{n} = \frac{\mathbf{n}}{\|\mathbf{n}\|}$ guarantees an orthonormal surface normal vector in $\mathcal{O}(1)$ parallel time without iterative branching.

---

## 4. Exact 7-Region Voronoi Point-to-Triangle Distance

Let a triangle be defined by vertices $A, B, C \in \mathbb{R}^3$, and let $P \in \mathbb{R}^3$ be a query point. The closest point $Q^* \in \triangle ABC$ satisfies:

$$Q^* = \arg\min_{Q \in \triangle ABC} \|P - Q\|^2$$

Let edge vectors be $\mathbf{u} = B - A$, $\mathbf{v} = C - A$, and displacement $\mathbf{w} = P - A$.
Define dot products:
$$d_1 = \mathbf{u} \cdot \mathbf{w}, \quad d_2 = \mathbf{v} \cdot \mathbf{w}$$
$$d_3 = \mathbf{u} \cdot (P - B), \quad d_4 = \mathbf{v} \cdot (P - B)$$
$$d_5 = \mathbf{u} \cdot (P - C), \quad d_6 = \mathbf{v} \cdot (P - C)$$

The 7 mutually exclusive geometric regions are:
1. **Vertex A**: If $d_1 \le 0$ and $d_2 \le 0 \implies Q^* = A$.
2. **Vertex B**: If $d_3 \ge 0$ and $d_4 \le d_3 \implies Q^* = B$.
3. **Vertex C**: If $d_6 \ge 0$ and $d_5 \le d_6 \implies Q^* = C$.
4. **Edge AB**: If $v_c \le 0$, $d_1 \ge 0$, and $d_3 \le 0 \implies Q^* = A + \left(\frac{d_1}{d_1 - d_3}\right) \mathbf{u}$.
5. **Edge AC**: If $v_b \le 0$, $d_2 \ge 0$, and $d_6 \le 0 \implies Q^* = A + \left(\frac{d_2}{d_2 - d_6}\right) \mathbf{v}$.
6. **Edge BC**: If $v_a \le 0$, $(d_4 - d_3) \ge 0$, and $(d_5 - d_6) \ge 0 \implies Q^* = B + \left(\frac{d_4 - d_3}{(d_4 - d_3) + (d_5 - d_6)}\right) (C - B)$.
7. **Face Interior**: Otherwise, with barycentric weights $(u_0, u_1, u_2) = \left(1 - \frac{v_b}{v_a + v_b + v_c} - \frac{v_c}{v_a + v_b + v_c}, \frac{v_b}{v_a + v_b + v_c}, \frac{v_c}{v_a + v_b + v_c}\right)$:
   $$Q^* = u_0 A + u_1 B + u_2 C$$

This formulation guarantees exact Euclidean projection with zero singular divisions.

---

## 5. Watertight Volume via Gauss's Divergence Theorem

For a closed, orientable, watertight polygonal surface $S = \partial V$, Gauss's Divergence Theorem states:

$$\iiint_V (\nabla \cdot \mathbf{F}) \, dV = \iint_S (\mathbf{F} \cdot \hat{\mathbf{n}}) \, dS$$

Selecting the identity vector field $\mathbf{F}(\mathbf{x}) = \mathbf{x} = (x, y, z)^T$:

$$\nabla \cdot \mathbf{F} = \frac{\partial x}{\partial x} + \frac{\partial y}{\partial y} + \frac{\partial z}{\partial z} = 1 + 1 + 1 = 3$$

Therefore:

$$3 \iiint_V dV = 3 V = \sum_{i=1}^M \iint_{T_i} (\mathbf{x} \cdot \hat{\mathbf{n}}_i) \, dS$$

For a planar triangle $T_i$ with vertices $\mathbf{v}_{0,i}, \mathbf{v}_{1,i}, \mathbf{v}_{2,i}$, the facet normal area vector is:

$$\mathbf{A}_i = \frac{1}{2} (\mathbf{v}_{1,i} - \mathbf{v}_{0,i}) \times (\mathbf{v}_{2,i} - \mathbf{v}_{0,i}) = \hat{\mathbf{n}}_i \cdot \text{Area}(T_i)$$

The inner product with triangle centroid $\bar{\mathbf{x}}_i = \frac{1}{3}(\mathbf{v}_{0,i} + \mathbf{v}_{1,i} + \mathbf{v}_{2,i})$ evaluates to the scalar triple product:

$$\iint_{T_i} (\mathbf{x} \cdot \hat{\mathbf{n}}_i) \, dS = \bar{\mathbf{x}}_i \cdot \mathbf{A}_i = \frac{1}{2} \mathbf{v}_{0,i} \cdot (\mathbf{v}_{1,i} \times \mathbf{v}_{2,i})$$

Summing over all $M$ surface facets yields the closed-form volume:

$$V = \frac{1}{6} \sum_{i=1}^M \mathbf{v}_{0,i} \cdot (\mathbf{v}_{1,i} \times \mathbf{v}_{2,i})$$

---

## 6. Point-to-Plane ICP Normal Equations

Given source points $\mathbf{p}_i$ and target planar points $\mathbf{q}_i$ with surface normals $\mathbf{n}_i$, the point-to-plane objective minimizes:

$$E(R, \mathbf{t}) = \sum_{i=1}^N \left( (R \mathbf{p}_i + \mathbf{t} - \mathbf{q}_i) \cdot \mathbf{n}_i \right)^2$$

Parameterizing small rigid rotation $\Delta \mathbf{R} \approx I + [\boldsymbol{\omega}]_\times$ where $\boldsymbol{\omega} = (\alpha, \beta, \gamma)^T$, the linearized residual is:

$$r_i(\boldsymbol{\xi}) = (\mathbf{p}_i \times \mathbf{n}_i) \cdot \boldsymbol{\omega} + \mathbf{n}_i \cdot \mathbf{t} - (\mathbf{q}_i - \mathbf{p}_i) \cdot \mathbf{n}_i$$

where $\boldsymbol{\xi} = (\alpha, \beta, \gamma, t_x, t_y, t_z)^T \in \mathbb{R}^6$.

Writing in matrix form $A \boldsymbol{\xi} = \mathbf{b}$, row $i$ of $A$ is:

$$\mathbf{a}_i^T = \begin{bmatrix} (\mathbf{p}_i \times \mathbf{n}_i)^T & \mathbf{n}_i^T \end{bmatrix} \in \mathbb{R}^{1 \times 6}$$
$$b_i = (\mathbf{q}_i - \mathbf{p}_i) \cdot \mathbf{n}_i \in \mathbb{R}$$

The normal equations are:

$$\left( \sum_{i=1}^N \mathbf{a}_i \mathbf{a}_i^T + \lambda I \right) \boldsymbol{\xi} = \sum_{i=1}^N b_i \mathbf{a}_i$$

Because $A^T A \in \mathbb{R}^{6 \times 6}$ is symmetric positive semi-definite, Caliper solves for $\boldsymbol{\xi}$ via damped Cholesky factorization ($L L^T = A^T A + \lambda I$), guaranteeing numerical stability even in degenerate cylindrical or planar geometries.
