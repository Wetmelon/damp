# API Reference

Auto-generated from `@brief` doc comments in `inc/damp/`. Regenerate with `python tools/gen_reference.py`. Flat A→Z view: [REFERENCE_INDEX.md](REFERENCE_INDEX.md).

Compile-time (or init-time) control design in the firmware tree (variant gains as `constexpr`); same runtime objects in SIL. Scope: [docs/known_limitations.md](docs/known_limitations.md).

`damp::design::` is the design-time shelf (not for PWM-rate): heavy solvers and Result structs. Domain packs (motor / power / full motion kits) are not assumed present in the public core tree.


- [API Reference](#api-reference)
  - [Core, configuration \& backend vocabulary](#core-configuration--backend-vocabulary)
  - [Scalar math, complex \& frames](#scalar-math-complex--frames)
  - [Linear algebra](#linear-algebra)
  - [LTI systems (SS / TF / ZPK / discretize)](#lti-systems-ss--tf--zpk--discretize)
  - [Design-time synthesis (not PWM-rate)](#design-time-synthesis-not-pwm-rate)
  - [Runtime controllers](#runtime-controllers)
  - [Observers \& estimators](#observers--estimators)
  - [Filters \& signal conditioning](#filters--signal-conditioning)
  - [Trajectory value types](#trajectory-value-types)
  - [Kinematics / pose](#kinematics--pose)
  - [Motor control pack (if present)](#motor-control-pack-if-present)
  - [Embedded helpers (controls-adjacent utilities)](#embedded-helpers-controls-adjacent-utilities)
  - [Frequency-domain analysis (host)](#frequency-domain-analysis-host)
  - [Simulation / SIL harness (host)](#simulation--sil-harness-host)
  - [MATLAB®-style aliases (host)](#matlab-style-aliases-host)
  - [Math backends](#math-backends)
  - [Examples](#examples)


> **Substrate — types every path uses**

## Core, configuration & backend vocabulary

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AnalyticDynamics`](inc/damp/concepts.hpp#L233) | Dynamics that provide their own Jacobians: f.jacobian(x, u) → DynamicsJacobian |
| [`AnalyticMeasurement`](inc/damp/concepts.hpp#L239) | A Measurement that provides its own Jacobians: h.jacobian(x, u) → MeasurementJacobian |
| [`Dynamics`](inc/damp/concepts.hpp#L195) | Continuous-time nonlinear dynamics ẋ = f(x, u): a callable returning ẋ as a ColVec<NX, T> |
| [`DynamicsJacobian`](inc/damp/concepts.hpp#L214) | Analytic Jacobians of a Dynamics model at (x, u), returned by its jacobian(x, u) |
| [`FeedbackLaw`](inc/damp/concepts.hpp#L118) | An LQ or MPC law that holds outputs at r: control(r) and step(r, m) |
| [`FeedbackRuntime`](inc/damp/concepts.hpp#L104) | Shared lifecycle of an LQ or MPC runtime: measurement, time update, applied input |
| [`Measurement`](inc/damp/concepts.hpp#L208) | Nonlinear measurement y = h(x, u): a callable returning y as a ColVec<NY, T> |
| [`MeasurementJacobian`](inc/damp/concepts.hpp#L221) | Analytic Jacobians of a Measurement model at (x, u), returned by its jacobian(x, u) |
| [`ParameterEstimator`](inc/damp/concepts.hpp#L159) | Online grey-box parameter estimator: physical parameters + gating |
| [`SignalSource`](inc/damp/concepts.hpp#L144) | Self-clocked signal source: u = step(), finished when done() |
| [`SisoController`](inc/damp/concepts.hpp#L87) | Scalar output-feedback controller: u = control(r, y), fixed rate |
| [`StateEstimator`](inc/damp/concepts.hpp#L131) | State estimator: x̂ = step(y, u) — the fused per-tick form |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`array_as`](inc/damp/backend.hpp#L195) | Element-wise static_cast of a damp::array to scalar type U (precision rebind) |
| [`minmax`](inc/damp/backend.hpp#L164) | Ordered {min, max} pair returned by value (+1 more overload) |

## Scalar math, complex & frames

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AlphaBeta`](inc/damp/math/transforms.hpp#L164) | Alpha-beta (stationary-frame) component pair |
| [`complex`](inc/damp/math/complex.hpp#L38) | Constexpr complex number class for compile-time computations |
| [`Convention`](inc/damp/math/transforms.hpp#L69) | Scaling convention for the Clarke/Park family |
| [`DCM`](inc/damp/math/geometry.hpp#L86) | Rotation matrix, active body→parent (SO(3) wrapper over Mat3) |
| [`DirectQuadrature`](inc/damp/math/transforms.hpp#L84) | Direct-quadrature (rotor-frame) component pair |
| [`Euler`](inc/damp/math/geometry.hpp#L226) | Intrinsic Euler angles for sequence Order (default ZYX aerospace YPR) |
| [`EulerOrder`](inc/damp/math/geometry.hpp#L55) | Intrinsic Euler angle sequence (axis order of successive principal rotations) |
| [`InstantaneousPower`](inc/damp/math/transforms.hpp#L457) | Instantaneous active and reactive power |
| [`Quaternion`](inc/damp/math/geometry.hpp#L362) | Unit quaternion rotation (w, x, y, z) (Hamilton product) |
| [`SequenceComponents`](inc/damp/math/transforms.hpp#L514) | Symmetrical (sequence) components of a three-phase phasor set |
| [`SplitReal`](inc/damp/math/split_real.hpp#L48) | A real number split into an exact integer part and a floating-point remainder |
| [`Transform4`](inc/damp/math/geometry.hpp#L785) | 4×4 homogeneous transform (SE(3) interop / DH export) |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`abs`](inc/damp/math/complex.hpp#L362) | Compute magnitude (absolute value) of a complex number (+1 more overload) |
| [`acos`](inc/damp/math/math.hpp#L152) | Arccosine ∈ [0, π]. Input is clamped to [−1, 1] in both paths |
| [`arg`](inc/damp/math/complex.hpp#L373) | Compute argument (phase angle) of a complex number |
| [`asin`](inc/damp/math/math.hpp#L133) | Arcsine ∈ [−π/2, π/2]. Input is clamped to [−1, 1] in both paths so behavior matches at compile and run time (std::asin would return NaN for \|x\| > 1) |
| [`atan`](inc/damp/math/math.hpp#L119) | Single-argument arctangent ∈ (−π/2, π/2) |
| [`atan2`](inc/damp/math/math.hpp#L107) | Two-argument arctangent, atan2(y, x) ∈ [−π, π] |
| [`cbrt`](inc/damp/math/math.hpp#L94) | Cube root (preserves sign for negative x) |
| [`ceil`](inc/damp/math/math.hpp#L383) | Ceiling — smallest integer ≥ x |
| [`clarke_park_transform`](inc/damp/math/transforms.hpp#L404) | Fused Clarke-Park transform (abc → dq) |
| [`clarke_transform`](inc/damp/math/transforms.hpp#L277) | Clarke transform (abc → αβ) |
| [`copysign`](inc/damp/math/math.hpp#L484) | Copy sign — magnitude of mag with the sign of sgn_src |
| [`cos`](inc/damp/math/math.hpp#L171) | Cosine |
| [`db2mag`](inc/damp/math/math.hpp#L552) | Decibels to magnitude, 10^(db/20) |
| [`deg2rad`](inc/damp/math/math.hpp#L574) | Degrees to radians, deg·π/180 |
| [`exp`](inc/damp/math/math.hpp#L234) | Exponential function |
| [`finite_non_negative`](inc/damp/math/math.hpp#L526) | True if x is finite and non-negative (x ≥ 0) |
| [`finite_positive`](inc/damp/math/math.hpp#L516) | True if x is finite and strictly greater than zero |
| [`floor`](inc/damp/math/math.hpp#L371) | Floor — largest integer ≤ x |
| [`fmod`](inc/damp/math/math.hpp#L415) | Floating-point remainder, x − y·trunc(x/y) (sign of x), matching std::fmod's truncated-quotient convention |
| [`hypot`](inc/damp/math/math.hpp#L78) | Euclidean distance hypot(x, y) = √(x² + y²), without overflow |
| [`instantaneous_power`](inc/damp/math/transforms.hpp#L496) | Instantaneous active and reactive power from dq quantities (+1 more overload) |
| [`inverse_clarke_transform`](inc/damp/math/transforms.hpp#L308) | Inverse Clarke transform (αβ → abc) |
| [`inverse_park_clarke_transform`](inc/damp/math/transforms.hpp#L423) | Fused inverse Park-Clarke transform (dq → abc) |
| [`inverse_park_transform`](inc/damp/math/transforms.hpp#L382) | Inverse Park transform (dq → αβ) |
| [`inverse_symmetrical_components`](inc/damp/math/transforms.hpp#L571) | Inverse symmetrical-component transform (012 → abc) |
| [`isfinite`](inc/damp/math/math.hpp#L503) | Finiteness test — false for NaN and ±∞ |
| [`log`](inc/damp/math/math.hpp#L252) | Natural logarithm |
| [`log10`](inc/damp/math/math.hpp#L434) | Base-10 logarithm, log10(x) = ln(x) / ln(10) |
| [`mag2db`](inc/damp/math/math.hpp#L538) | Magnitude to decibels, 20·log10(mag) |
| [`nearbyint`](inc/damp/math/math.hpp#L399) | Round to nearest integer; ties to even (IEEE default / `FE_TONEAREST`) |
| [`park_transform`](inc/damp/math/transforms.hpp#L355) | Park transform (αβ → dq) |
| [`pow`](inc/damp/math/math.hpp#L282) | Power function with a real exponent, base^exponent (+1 more overload) |
| [`rad2deg`](inc/damp/math/math.hpp#L563) | Radians to degrees, rad·180/π |
| [`sgn`](inc/damp/math/math.hpp#L451) | Sign function — −1 if val < 0, 1 if val > 0, 0 if val == 0 |
| [`sin`](inc/damp/math/math.hpp#L184) | Sine |
| [`sincos`](inc/damp/math/math.hpp#L202) | Combined sine and cosine, {sin(x), cos(x)} |
| [`smooth_sign`](inc/damp/math/math.hpp#L464) | Smoothed sign: an odd quintic through ±width, exactly ±1 outside it |
| [`sqrt`](inc/damp/math/complex.hpp#L337) | Compute complex square root (constexpr) (+1 more overload) |
| [`symmetrical_components`](inc/damp/math/transforms.hpp#L542) | Forward symmetrical-component (Fortescue) transform (abc → 012) |
| [`tan`](inc/damp/math/math.hpp#L217) | Tangent |
| [`wrap`](inc/damp/math/math.hpp#L596) | Wrap x into the half-open interval [min, max) (period max − min) |
| [`wrap_pi`](inc/damp/math/math.hpp#L639) | Wrap an angle to [−π, π) |
| [`wrap_two_pi`](inc/damp/math/math.hpp#L657) | Wrap an angle to [0, 2π) |
| [`zero_sequence`](inc/damp/math/transforms.hpp#L250) | Zero-sequence (common-mode) scalar of a three-phase set |

## Linear algebra

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`Block`](inc/damp/matrix/block.hpp#L38) | Block view (non-owning) into a parent matrix |
| [`ColVec`](inc/damp/matrix/colvec.hpp#L30) | Concrete Column vector specialization of Matrix<N, 1, T> |
| [`ColView`](inc/damp/matrix/views.hpp#L332) | Non-owning column view of a matrix |
| [`ConvertsExactly`](inc/damp/matrix/matrix_traits.hpp#L80) | Matrix elements of type From convert to To without rounding: the same type, or a real into the complex of its own precision |
| [`Diagonal`](inc/damp/matrix/views.hpp#L31) | Diagonal view of a square matrix |
| [`EigenResult`](inc/damp/matrix/eigen.hpp#L39) | Eigenvalue computation result |
| [`FullQR`](inc/damp/matrix/decomposition.hpp#L515) | Result of a full (complete) QR factorization |
| [`GeneralizedEigenResult`](inc/damp/matrix/generalized_eigen.hpp#L61) | Generalized eigenvalues of A − λB: the finite ones, and how many there are |
| [`is_matrix_element`](inc/damp/matrix/matrix_traits.hpp#L42) | True if T may be a Matrix element type |
| [`LDLT`](inc/damp/matrix/decomposition.hpp#L163) | LDLᴴ factors: A = L D Lᴴ with L unit lower-triangular and D real diagonal |
| [`LowerTriangle`](inc/damp/matrix/views.hpp#L165) | Lower triangular view of a square matrix |
| [`Matrix`](inc/damp/matrix/core.hpp#L205) | Fixed-size, stack-allocated matrix for linear algebra operations |
| [`MatrixLike`](inc/damp/matrix/matrix_traits.hpp#L98) | Concept for any type that provides 2D matrix-like element access |
| [`MatrixLikeOf`](inc/damp/matrix/matrix_traits.hpp#L108) | Concept for a MatrixLike type with specific dimensions |
| [`NullSpace`](inc/damp/matrix/svd.hpp#L437) | Orthonormal basis for the null space (kernel) of a matrix |
| [`QRDecomposition`](inc/damp/matrix/decomposition.hpp#L397) | Thin QR factorization A = QR (Householder) |
| [`RowVec`](inc/damp/matrix/rowvec.hpp#L31) | Row vector specialization of Matrix<1, N, T> |
| [`RowView`](inc/damp/matrix/views.hpp#L203) | Non-owning row view of a matrix |
| [`SamePrecision`](inc/damp/matrix/matrix_traits.hpp#L87) | Two matrix element types mix in one operator without rounding: equal, or a real and the complex of its own precision |
| [`SVDResult`](inc/damp/matrix/svd.hpp#L49) | Result of a full singular value decomposition A = U·Σ·Vᴴ |
| [`TransposeView`](inc/damp/matrix/views.hpp#L456) | Non-owning transpose view of a matrix (zero-copy) |
| [`UpperTriangle`](inc/damp/matrix/views.hpp#L124) | Upper triangular view of a square matrix |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`all_of`](inc/damp/matrix/core.hpp#L1421) | Whether p holds for every entry of A: the matrix counterpart of std::all_of |
| [`any_of`](inc/damp/matrix/core.hpp#L1445) | Whether p holds for some entry of A: the matrix counterpart of std::any_of |
| [`backward_substitute_transpose`](inc/damp/matrix/solve.hpp#L61) | Backward substitution for Lᴴx = b (Lᵀx = b for real T) |
| [`blkdiag`](inc/damp/matrix/concat.hpp#L135) | Block-diagonal assembly of MatrixLike arguments |
| [`cholesky`](inc/damp/matrix/decomposition.hpp#L124) | Cholesky decomposition for positive-definite matrices |
| [`cholesky_solve`](inc/damp/matrix/solve.hpp#L159) | Solve AX = B via Cholesky (A = LLᴴ) |
| [`compute_eigenvalues`](inc/damp/matrix/eigen.hpp#L763) | Compute the eigenvalues (and Schur vectors) of a real square matrix |
| [`cos`](inc/damp/matrix/functions.hpp#L1353) | Matrix cosine via scaling and double-angle reconstruction |
| [`cosh`](inc/damp/matrix/functions.hpp#L1387) | Matrix hyperbolic cosine cosh(A) = (exp(A) + exp(−A))/2 |
| [`default_tol`](inc/damp/matrix/matrix_traits.hpp#L125) | Type-appropriate default tolerance for floating-point comparisons. float  ~7 decimal digits  → 1e-6 double ~15 decimal digits → 1e-12 |
| [`det`](inc/damp/matrix/functions.hpp#L172) | Matrix determinant det(A) |
| [`diag`](inc/damp/matrix/concat.hpp#L55) | Square diagonal matrix from scalars: diag(10.0, 0.1, 100.0, 10.0) (+4 more overloads) |
| [`eps`](inc/damp/matrix/matrix_traits.hpp#L115) | Machine epsilon of T's real scalar (MATLAB® eps). Scale it by problem size and magnitude to build a relative tolerance |
| [`expm`](inc/damp/matrix/functions.hpp#L696) | Matrix exponential by scaling and squaring (Al-Mohy & Higham 2009) |
| [`forward_substitute`](inc/damp/matrix/solve.hpp#L37) | Forward substitution for Lx = b |
| [`full_qr`](inc/damp/matrix/decomposition.hpp#L533) | Full QR factorization via Householder reflections (real or complex T) |
| [`generalized_eigenvalues`](inc/damp/matrix/generalized_eigen.hpp#L509) | Generalized eigenvalues of the real pencil A − λB: the λ with det(A − λB) = 0 |
| [`hcat`](inc/damp/matrix/concat.hpp#L159) | Horizontal concatenation (same row count) |
| [`infinity_norm`](inc/damp/matrix/functions.hpp#L42) | Infinity norm ‖A‖∞: maximum absolute row sum |
| [`inverse`](inc/damp/matrix/solve.hpp#L294) | Matrix inverse A⁻¹ via mat::solve(A, I) |
| [`is_symmetric_or_hermitian`](inc/damp/matrix/decomposition.hpp#L91) | True if A equals its transpose (real) or conjugate transpose (complex) |
| [`ldlt`](inc/damp/matrix/decomposition.hpp#L191) | LDLᴴ decomposition of a symmetric (Hermitian) matrix, no pivoting |
| [`ldlt_solve`](inc/damp/matrix/solve.hpp#L189) | Solve AX = B via LDLᴴ (A = L D Lᴴ, symmetric indefinite) |
| [`log`](inc/damp/matrix/functions.hpp#L1031) | Principal matrix logarithm via inverse scaling and squaring |
| [`lu_decomposition`](inc/damp/matrix/decomposition.hpp#L254) | LU decomposition with partial pivoting |
| [`lu_solve`](inc/damp/matrix/solve.hpp#L220) | Solve AX = B via LU with partial pivoting |
| [`null_space`](inc/damp/matrix/svd.hpp#L458) | Orthonormal basis for the null space {x : A·x = 0} via SVD |
| [`one_norm`](inc/damp/matrix/functions.hpp#L71) | One-norm ‖A‖₁: maximum absolute column sum |
| [`pow`](inc/damp/matrix/functions.hpp#L1133) | Integer matrix power via binary exponentiation (+1 more overload) |
| [`pseudo_inverse`](inc/damp/matrix/svd.hpp#L404) | Moore–Penrose pseudoinverse A⁺ via SVD |
| [`qr_decompose`](inc/damp/matrix/decomposition.hpp#L444) | Thin QR via Householder reflections |
| [`quadratic_form`](inc/damp/matrix/core.hpp#L1246) | Symmetric congruence (quadratic) form  S = M X Mᵀ |
| [`rank`](inc/damp/matrix/functions.hpp#L234) | Numerical rank: the number of singular values above rel_tol·σₘₐₓ |
| [`rank_from_svd`](inc/damp/matrix/svd.hpp#L368) | Numerical rank from a precomputed SVD result |
| [`sin`](inc/damp/matrix/functions.hpp#L1336) | Matrix sine via scaling and double-angle reconstruction |
| [`sincos`](inc/damp/matrix/functions.hpp#L1320) | Compute sin(A) and cos(A) together via scaling and double-angle reconstruction |
| [`sinh`](inc/damp/matrix/functions.hpp#L1371) | Matrix hyperbolic sine sinh(A) = (exp(A) − exp(−A))/2 |
| [`solve`](inc/damp/matrix/solve.hpp#L266) | Solve AX = B — Cholesky when Hermitian, else LU (+2 more overloads) |
| [`sqrt`](inc/damp/matrix/functions.hpp#L986) | Principal matrix square root by the real Schur method |
| [`svd`](inc/damp/matrix/svd.hpp#L334) | Full singular value decomposition A = U·Σ·Vᴴ (one-sided Jacobi) |
| [`transform`](inc/damp/matrix/core.hpp#L1397) | Apply f to every entry, B(i, j) = f(A(i, j)): the matrix counterpart of std::transform |
| [`two_norm`](inc/damp/matrix/functions.hpp#L107) | Spectral norm ‖A‖₂ = σₘₐₓ(A) |
| [`vcat`](inc/damp/matrix/concat.hpp#L182) | Vertical concatenation (same column count) |

## LTI systems (SS / TF / ZPK / discretize)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`DiscretizationMethod`](inc/damp/systems/discretization.hpp#L31) | Discretization methods for continuous-time state-space systems |
| [`DiscretizationResult`](inc/damp/systems/discretization.hpp#L54) | A discretized or resampled model with its success / status, read like an optional |
| [`InputDelay`](inc/damp/systems/discretization.hpp#L557) | An input delay of θ = seconds, in up to N registers (0 ≤ θ ≤ N·Ts): each register state holds an input already sent, v(k−1) … v(k−N). `discretize(sys, Ts, InputDelay<N>{theta})` returns that model; design on it like any other. `InputDelay{theta}` is N = 1 (+1 more overload) |
| [`LeverrierResult`](inc/damp/systems/zpk.hpp#L519) | Faddeev–LeVerrier characteristic polynomial and adjoint coefficient matrices |
| [`PolyRootsResult`](inc/damp/systems/zpk.hpp#L93) | Roots of a real polynomial given highest power first |
| [`SsFeed`](inc/damp/systems/state_space.hpp#L646) | One term of a linear output mix: k · y_iy |
| [`StaircaseFormResult`](inc/damp/systems/realization.hpp#L646) | A model in staircase coordinates, the orthogonal transform that took it there, and its rank |
| [`StateScaling`](inc/damp/systems/realization.hpp#L40) | A model in balanced state coordinates, and the diagonal scaling that took it there |
| [`StateSpace`](inc/damp/systems/state_space.hpp#L139) | State-space representation for linear time-invariant systems (discrete or continuous) (+10 more overloads) |
| [`StateSpaceZPKResult`](inc/damp/systems/zpk.hpp#L663) | SISO state-space → ZPK conversion result with runtime zero count |
| [`StateTransformResult`](inc/damp/systems/realization.hpp#L124) | A model in new state coordinates z = transform·x, and the transform that took it there |
| [`TransferFunction`](inc/damp/systems/transfer_function.hpp#L102) | SISO polynomial transfer function G(s) = num(s)/den(s) |
| [`TransmissionZerosResult`](inc/damp/systems/zeros.hpp#L58) | The finite invariant zeros of a state-space model and their count |
| [`ZPK`](inc/damp/systems/zpk.hpp#L199) | Zero-pole-gain (ZPK) representation of a SISO LTI system |
| [`ZPKResult`](inc/damp/systems/zpk.hpp#L434) | Result of converting a transfer function to zero-pole-gain form |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`accumulate_poly`](inc/damp/systems/transfer_function.hpp#L38) | dst += scale·src, aligned on the constant term (the last entry: coefficients are highest power first) |
| [`balance_states`](inc/damp/systems/realization.hpp#L82) | Balance a model by a diagonal state scaling, so A, B and C have comparable magnitudes |
| [`companion_form`](inc/damp/systems/realization.hpp#L617) | Companion form of a SISO model, MATLAB®'s canon(sys, 'companion'): B = e₁, −α in A's last column |
| [`controllability_staircase`](inc/damp/systems/realization.hpp#L792) | Controllability staircase form: the controllable part split off by an orthogonal transform |
| [`controllable_canonical_form`](inc/damp/systems/realization.hpp#L519) | Controllable canonical form of a SISO model: ones on A's superdiagonal, −α in its last row, B = eₙ |
| [`discretize`](inc/damp/systems/discretization.hpp#L334) | Discretize a continuous-time state-space system (+1 more overload) |
| [`eval_frf`](inc/damp/systems/state_space.hpp#L234) | Evaluate frequency response of a state-space system |
| [`feedback`](inc/damp/systems/state_space.hpp#L392) | Negative feedback: y = sys1(u − sys2(y)) (+1 more overload) |
| [`feedback_read`](inc/damp/systems/state_space.hpp#L829) | Unity (or weighted) negative feedback around G_loop, read G_out |
| [`fractional_delay`](inc/damp/systems/discretization.hpp#L515) | Fractional-delay all-pass filter of order N (Thiran): z⁻ᴰ for a real D = delay / Ts |
| [`leverrier_num_scale`](inc/damp/systems/zpk.hpp#L585) | Magnitude scale of each strictly proper Leverrier numerator coefficient |
| [`lft`](inc/damp/systems/state_space.hpp#L506) | Lower linear fractional transformation: K closes the loop from P's last outputs to its last inputs |
| [`markov_checked_numerator`](inc/damp/systems/zpk.hpp#L615) | num, or the numerator the Markov parameters give over den when num misses it |
| [`minreal_zpk`](inc/damp/systems/zpk.hpp#L836) | Cancel matching pole-zero pairs on a ZPK model |
| [`mix`](inc/damp/systems/state_space.hpp#L711) | Fold plant outputs into one input: v = g u + Σ kᵢ yᵢ (+1 more overload) |
| [`mix_outputs`](inc/damp/systems/state_space.hpp#L788) | Linear combination of outputs: y = Σ αⱼ yⱼ |
| [`modal_form`](inc/damp/systems/realization.hpp#L286) | Real modal canonical form: A block diagonal, one block per mode |
| [`mul_linear`](inc/damp/systems/zpk.hpp#L165) | Multiply a polynomial of degree by (s − root) in place; poly[k] multiplies sᵏ |
| [`observability_staircase`](inc/damp/systems/realization.hpp#L851) | Observability staircase form: the observable part split off by an orthogonal transform |
| [`observable_canonical_form`](inc/damp/systems/realization.hpp#L569) | Observable canonical form of a SISO model: ones on A's subdiagonal, −α in its last column, C = eₙᵀ |
| [`parallel`](inc/damp/systems/state_space.hpp#L325) | Parallel connection (shared input, summed outputs) |
| [`poly_roots`](inc/damp/systems/zpk.hpp#L124) | Compute roots of a real polynomial given highest power first |
| [`resample`](inc/damp/systems/discretization.hpp#L404) | Resample a ZOH-discrete model to a new sample time |
| [`schur_stable`](inc/damp/systems/transfer_function.hpp#L64) | True if every root of z^(N−1)·p(z⁻¹) lies strictly inside the unit circle (Schur–Cohn) |
| [`select`](inc/damp/systems/state_space.hpp#L670) | SISO channel slice y_iy / u_iu |
| [`series`](inc/damp/systems/state_space.hpp#L269) | Series connection: sys2 follows sys1 (u → sys1 → sys2 → y) |
| [`similarity_transform`](inc/damp/systems/realization.hpp#L375) | Change the state basis to z = T·x: (T·A·T⁻¹, T·B, C·T⁻¹, D), and T·G for the noise input |
| [`ss2tf`](inc/damp/systems/zpk.hpp#L815) | MATLAB®-style alias for SISO to_transfer_function |
| [`ss_gain`](inc/damp/systems/state_space.hpp#L657) | Static gain y = k u (no states) |
| [`subtract`](inc/damp/systems/state_space.hpp#L588) | Differencing connection: outputs y = y₁ − y₂ |
| [`to_transfer_function`](inc/damp/systems/zpk.hpp#L785) | SISO state-space → transfer function, size NX+1 / NX+1 |
| [`to_zpk`](inc/damp/systems/zpk.hpp#L470) | Convert a SISO transfer function to zero-pole-gain form (+1 more overload) |
| [`transmission_zeros`](inc/damp/systems/zeros.hpp#L428) | The finite invariant zeros of a state-space model with any number of inputs and outputs |
| [`upsample`](inc/damp/systems/discretization.hpp#L464) | Upsample a discrete model by an integer factor: sys(z) → sys(zᴸ) at Ts / L |

> **Design-time — not for the PWM interrupt**

## Design-time synthesis (not PWM-rate)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`ActiveSetSolver`](inc/damp/design/qp.hpp#L703) | Default QP solver policy: the Goldfarb–Idnani active-set solve_qp() |
| [`AdmmSettings`](inc/damp/design/qp.hpp#L804) | Tuning parameters for the ADMM QP solver |
| [`AdmmSolver`](inc/damp/design/qp.hpp#L885) | ADMM (OSQP-style) QP solver policy — warm-started, fixed-cost iterations |
| [`ADRCResult`](inc/damp/controllers/adrc.hpp#L51) | Active Disturbance Rejection Control design result |
| [`BalancedRealizationResult`](inc/damp/design/model_reduction.hpp#L77) | Balanced realization result (Moore square-root method) |
| [`BalancedReductionResult`](inc/damp/design/model_reduction.hpp#L109) | Reduced-order model from balanced truncation or residualization |
| [`CommandProjectionResult`](inc/damp/controllers/action_governor.hpp#L91) | Result of an affine command projection |
| [`DareMethod`](inc/damp/design/riccati.hpp#L1189) | DARE solver selection for dare() |
| [`DesignStatus`](inc/damp/design/status.hpp#L35) | Why a design failed, one value per cause class (not per code branch) |
| [`DisturbanceModel`](inc/damp/design/disturbance.hpp#L82) | An explicit disturbance model: d enters the plant's state through B_d and the readings through C_d |
| [`ESCResult`](inc/damp/controllers/esc.hpp#L91) | Extremum-seeking design result: the discrete realization the runtime runs |
| [`HarmonicSuppressorResult`](inc/damp/controllers/harmonic_suppression.hpp#L50) | Design result for a multi-resonant harmonic suppressor |
| [`InputDisturbance`](inc/damp/design/disturbance.hpp#L42) | An input disturbance: d adds to u, B_d = B and C_d = D, one disturbance per input |
| [`InteriorPointSolver`](inc/damp/design/qp.hpp#L1446) | Interior-point QP solver policy for damp::MPC |
| [`InvariantSetResult`](inc/damp/design/invariant_set.hpp#L41) | Maximal positively invariant polytope G e ≤ h of x⁺ = A x |
| [`JordanBlock`](inc/damp/design/pole_placement.hpp#L488) | One Jordan mini-block of a desired closed-loop spectrum |
| [`JordanObjective`](inc/damp/design/pole_placement.hpp#L1200) | How place_jordan spends the free parameter of the Klein–Moore family (Schmid et al.) |
| [`JordanOptions`](inc/damp/design/pole_placement.hpp#L1213) | Options for place_jordan |
| [`JordanPlaceResult`](inc/damp/design/pole_placement.hpp#L1228) | place_jordan result: the PlaceResult fields plus the achieved robustness metrics |
| [`Law`](inc/damp/design/model_inverse.hpp#L884) | Law<kRes, T> of a model inverse: StateFeedforward or OutputInverse, default-constructed from the design |
| [`LeadLagResult`](inc/damp/controllers/lead_lag.hpp#L64) | Lead-lag compensator design result |
| [`LeadLagSeriesResult`](inc/damp/controllers/lead_lag.hpp#L224) | Cascaded lead+lag design result (2nd-order StateSpace + success) |
| [`LinearizationResult`](inc/damp/design/linearization.hpp#L46) | Result of nonlinear operating-point linearization |
| [`LQRCost`](inc/damp/design/riccati.hpp#L1577) | A quadratic cost on (x, u): state cost Q, input cost R, cross term N |
| [`LQRResult`](inc/damp/design/riccati.hpp#L1518) | The LQR gain: design::lqr on a plant, discrete_lqr and continuous_lqr on given matrices |
| [`MinimalRealizationResult`](inc/damp/design/minreal.hpp#L89) | Minimal realization result |
| [`ModelReductionMethod`](inc/damp/design/model_reduction.hpp#L61) | How balanced_reduction and eliminate_states eliminate states |
| [`MPCResult`](inc/damp/controllers/mpc.hpp#L1220) | An MPC design: the condensed QP, the estimator when there is one, and the steady-state target the runtime poses the QP about; consumed by damp::MPC |
| [`NoDisturbance`](inc/damp/design/disturbance.hpp#L33) | No disturbance model (the default): the plain filter on x |
| [`option_set`](inc/damp/design/options.hpp#L83) | The kinds of option a design call accepts: check<Opts...>() refuses, at compile time, an option of no listed kind and a kind given twice |
| [`option_type`](inc/damp/design/options.hpp#L42) | The type of the first option of the kind Is among Opts, else Default (+1 more overload) |
| [`OutputDisturbance`](inc/damp/design/disturbance.hpp#L53) | An output disturbance: a bias on each reading, B_d = 0 and C_d = I, one per output |
| [`OutputInverse`](inc/damp/design/model_inverse.hpp#L804) | The stable inverse at run time: the reference point (x , u ) from the output reference over the preview, seven matrix–vector products and the zero dynamics' state (+1 more overload) |
| [`OutputInverseResult`](inc/damp/design/model_inverse.hpp#L142) | The stable inverse of a discrete model over a preview of NP output samples |
| [`PIDPolePlacementSpec`](inc/damp/design/pid_design.hpp#L912) | Closed-loop targets for the continuous pid_pole_placement forms |
| [`PIDTuneSpec`](inc/damp/design/pid_design.hpp#L449) | Knobs for plant-aware pidtune |
| [`PlaceResult`](inc/damp/design/pole_placement.hpp#L1085) | Pole-placement design result: the gain K, the achieved closed-loop poles e, and Ts |
| [`PlantDisturbance`](inc/damp/design/disturbance.hpp#L67) | The plant's own disturbance channel: B_d = G and C_d = H, one disturbance per column of G (H needs as many columns) |
| [`PRResult`](inc/damp/controllers/pr.hpp#L111) | Proportional-Resonant controller design result |
| [`QPResult`](inc/damp/design/qp.hpp#L84) | Result of a dense QP solve |
| [`QPStatus`](inc/damp/design/qp.hpp#L51) | Termination status of a QP solve |
| [`RepetitiveResult`](inc/damp/controllers/repetitive.hpp#L140) | Design result for the repetitive controller |
| [`SMCResult`](inc/damp/controllers/smc.hpp#L37) | Tuning parameters for a first-order sliding-mode controller |
| [`SmithPredictorResult`](inc/damp/controllers/smith_predictor.hpp#L97) | Smith predictor design result (discrete PID + FO model + delay samples) |
| [`StateFeedbackResult`](inc/damp/controllers/state_feedback_core.hpp#L672) | A state-feedback design: the gain, and the parts the law's options add |
| [`StateFeedforward`](inc/damp/design/model_inverse.hpp#L762) | The plant-inversion feedforward at run time: u = B⁺·(x ₖ₊₁ − x ₖ) − B⁺(A − I)·x ₖ, two matrix–vector products (+1 more overload) |
| [`StateFeedforwardResult`](inc/damp/design/model_inverse.hpp#L62) | The plant-inversion feedforward of a discrete model: u ₖ = B⁺·(x ₖ₊₁ − x ₖ) − B⁺(A − I)·x ₖ |
| [`StateSpaceResult`](inc/damp/design/model_result.hpp#L59) | A designed state-space model with its success / status |
| [`SteadyStateTarget`](inc/damp/design/steady_state_target.hpp#L372) | The steady-state target at run time: (x_t, u_t) from the reference and the estimated disturbance, two matrix–vector products, and the re-anchoring of a Luenberger observer on it (+1 more overload) |
| [`SteadyStateTargetResult`](inc/damp/design/steady_state_target.hpp#L75) | The steady-state target maps: [x_t; u_t] per unit disturbance and per unit reference |
| [`STSMCResult`](inc/damp/controllers/stsmc.hpp#L34) | Super-twisting (second-order sliding-mode) controller design result |
| [`TargetPoint`](inc/damp/design/steady_state_target.hpp#L356) | A state and the input at it: a steady state (SteadyStateTarget), or a point of a reference trajectory (OutputInverse) |
| [`TrackedOutputs`](inc/damp/design/steady_state_target.hpp#L62) | The outputs a law without integrators or a disturbance model holds at r through its steady-state target (design::lqr, design::lqg): `TrackedOutputs<1>{{0}}`, in that order (the references r follow it) |
| [`TransferFunctionResult`](inc/damp/design/model_result.hpp#L35) | A designed transfer function with its success / status |
| [`WarmStartActiveSetSolver`](inc/damp/design/qp.hpp#L739) | Warm-started active-set solver policy — the damp::MPC default |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`adrc`](inc/damp/controllers/adrc.hpp#L206) | Active Disturbance Rejection Control design |
| [`ADRCController`](inc/damp/controllers/adrc.hpp#L376) | Deduce the runtime from its design, keeping the result's scalar |
| [`amigo_kappa_tau`](inc/damp/design/pid_design.hpp#L253) | AMIGO PI from ultimate gain/period and static gain Kₛ |
| [`balanced_realization`](inc/damp/design/model_reduction.hpp#L478) | Balanced realization via the square-root (Moore / Laub) method |
| [`balanced_reduction`](inc/damp/design/model_reduction.hpp#L555) | Balanced reduction to NR states, by truncation or DC-matched residualization |
| [`care`](inc/damp/design/riccati.hpp#L1491) | Solve the Continuous-time Algebraic Riccati Equation (CARE) |
| [`cbf_relative_degree_1`](inc/damp/controllers/action_governor.hpp#L205) | Build a relative-degree-1 CBF inequality row |
| [`closed_loop_poles`](inc/damp/design/stability.hpp#L564) | Compute closed-loop poles (eigenvalues) with state feedback |
| [`closed_loop_stable`](inc/damp/design/pid_design.hpp#L684) | Closed-loop certificate: den_P·den_C + num_P·num_C is Hurwitz for C_fb = Kp + Kᵢ/s + Kd s/(1 + Tf s) |
| [`cohen_coon`](inc/damp/design/pid_design.hpp#L300) | Cohen-Coon tuning from first-order-plus-dead-time model |
| [`continuous_lqr`](inc/damp/design/riccati.hpp#L1744) | Continuous-time Linear-Quadratic Regulator design |
| [`ContinuousPID`](inc/damp/controllers/pid.hpp#L1071) | Deduce the runtime from its design, keeping the result's scalar |
| [`controllability_gramian`](inc/damp/design/stability.hpp#L382) | Continuous/discrete controllability Gramian W_c |
| [`controllability_matrix`](inc/damp/design/stability.hpp#L55) | Compute the controllability matrix [B, AB, A²B, ..., A^(N-1)B] |
| [`current_loop_pi`](inc/damp/design/pid_design.hpp#L1298) | PI gains for an R–L current loop 1/(L s + R) at any bandwidth ω |
| [`damping_ratio_from_overshoot_percent`](inc/damp/design/pid_design.hpp#L858) | Map percent overshoot target to equivalent damping ratio |
| [`dare`](inc/damp/design/riccati.hpp#L1322) | Solve the Discrete Algebraic Riccati Equation (DARE) |
| [`decltype`](inc/damp/controllers/state_feedback_core.hpp#L510) | Law<kRes, T> of an LQ or gain design: the runtime it runs as (LQR, LQI, LQG, LQGI, OFLQG), default-constructed from the design |
| [`discrete_lqr`](inc/damp/design/riccati.hpp#L1610) | Discrete-time Linear-Quadratic Regulator design (+1 more overload) |
| [`discrete_lyapunov_cholesky`](inc/damp/design/lyapunov.hpp#L510) | Cholesky factor of the discrete Lyapunov solution: A X Aᵀ − X + B Bᵀ = 0, X = RᵀR |
| [`discretize_lqr_cost`](inc/damp/design/riccati.hpp#L1679) | Discretize a continuous LQR cost integral over one sample (Van Loan) |
| [`dlyap`](inc/damp/design/lyapunov.hpp#L470) | Solve the discrete-time Lyapunov (Stein) equation A X Aᵀ − X + Q = 0 |
| [`eliminate_states`](inc/damp/design/model_reduction.hpp#L627) | Eliminate the listed states, by DC-matched residualization or truncation |
| [`esc`](inc/damp/controllers/esc.hpp#L156) | Synthesize an extremum-seeking controller |
| [`esc_mppt`](inc/damp/controllers/esc.hpp#L200) | MPPT-flavored ESC: maximize a power measurement by perturbing the operating point (e.g. converter duty or reference voltage) |
| [`ExtremumSeekingController`](inc/damp/controllers/esc.hpp#L374) | Deduce the runtime from its design, keeping the result's scalar |
| [`fit_at_crossover`](inc/damp/design/pid_design.hpp#L603) | fit_to_response at spec.wc with every gain of sign σ = ±1 |
| [`fit_to_response`](inc/damp/design/pid_design.hpp#L549) | The requested law with C(jω) = e^{j(−π + PM)}/P(jω), from the plant response P at ω = spec.wc |
| [`hankel_singular_values`](inc/damp/design/model_reduction.hpp#L508) | Hankel singular values of a stable state-space system |
| [`harmonic_suppressor`](inc/damp/controllers/harmonic_suppression.hpp#L128) | Synthesize a multi-resonant harmonic suppressor |
| [`HarmonicSuppressor`](inc/damp/controllers/harmonic_suppression.hpp#L273) | Deduce the runtime from its design, keeping the result's scalar |
| [`hurwitz`](inc/damp/design/pid_design.hpp#L628) | Routh–Hurwitz: every root of c[0] sⁿ + … + c[n] in the open left half plane |
| [`is_controllable`](inc/damp/design/stability.hpp#L447) | Check if a system is controllable |
| [`is_observable`](inc/damp/design/stability.hpp#L461) | Check if a system is observable |
| [`is_stabilizable`](inc/damp/design/stability.hpp#L484) | Check if (A, B) of sys is a stabilizable pair, in the domain of sys |
| [`is_stable`](inc/damp/design/stability.hpp#L502) | Stability of a state-space model, in the domain of sys |
| [`lag`](inc/damp/controllers/lead_lag.hpp#L197) | Design a lag compensator from desired low-frequency gain boost |
| [`lambda_tuning`](inc/damp/design/pid_design.hpp#L429) | Lambda tuning for FOPDT model |
| [`lead`](inc/damp/controllers/lead_lag.hpp#L144) | Design a lead compensator from desired phase boost at a target frequency |
| [`lead_lag`](inc/damp/controllers/lead_lag.hpp#L258) | Design a lead-lag compensator (cascade of lead + lag sections) |
| [`lead_lag_direct`](inc/damp/controllers/lead_lag.hpp#L293) | Direct lead-lag specification from zero/pole locations |
| [`LeadLagController`](inc/damp/controllers/lead_lag.hpp#L436) | Deduce the runtime from its design, keeping the result's scalar |
| [`linearize`](inc/damp/design/linearization.hpp#L118) | Linearize a Dynamics / Measurement model about an operating point (+1 more overload) |
| [`lqg`](inc/damp/controllers/lqg.hpp#L42) | The linear-quadratic Gaussian law: LQR on a steady-state Kalman filter's estimate (design::kalman(sys, noise)): runs as LQG (+3 more overloads) |
| [`LQG`](inc/damp/controllers/lqg.hpp#L230) | LQG law{res}: the sizes, scalar and filter of design::lqg's result |
| [`lqg_from_parts`](inc/damp/controllers/lqg.hpp#L110) | Assemble an LQG design from separately computed Kalman and LQR designs |
| [`lqgi`](inc/damp/controllers/lqgi.hpp#L39) | The linear-quadratic Gaussian law with integral action on every output: runs as LQGI (+3 more overloads) |
| [`LQGI`](inc/damp/controllers/lqgi.hpp#L147) | LQGI law{res}: the sizes, scalar and filter of design::lqgi's result |
| [`lqi`](inc/damp/controllers/lqi.hpp#L39) | The linear-quadratic regulator with integral action on every output (xi += Ts·(r − Ci·y)): runs as LQI (+3 more overloads) |
| [`LQI`](inc/damp/controllers/lqi.hpp#L136) | LQI law{res}: the sizes and scalar of design::lqi's result |
| [`lqr`](inc/damp/controllers/lqr.hpp#L45) | The linear-quadratic regulator on the plant (DARE on a discrete plant, CARE on a continuous one): runs as LQR (+3 more overloads) |
| [`LQR`](inc/damp/controllers/lqr.hpp#L204) | LQR law{res}: the sizes and scalar of design::lqr's result (+2 more overloads) |
| [`lqr_gain`](inc/damp/design/riccati.hpp#L1359) | Optimal LQR state-feedback gain from a Riccati solution |
| [`LQRWeights`](inc/damp/controllers/state_feedback_core.hpp#L623) | Sizes from the matrices: `LQRWeights{Q, R}` (+2 more overloads) |
| [`lyap`](inc/damp/design/lyapunov.hpp#L444) | Solve the continuous-time Lyapunov equation A X + X Aᵀ + Q = 0 |
| [`lyapunov_cholesky`](inc/damp/design/lyapunov.hpp#L492) | Cholesky factor of the continuous Lyapunov solution: A X + X Aᵀ + B Bᵀ = 0, X = RᵀR |
| [`maximal_invariant_set`](inc/damp/design/invariant_set.hpp#L251) | Maximal positively invariant set of x⁺ = A x under G₀ x ≤ h₀ |
| [`MeasuredDisturbance`](inc/damp/controllers/mpc.hpp#L1207) | Sizes from the matrix: `MeasuredDisturbance{Bd}` |
| [`minreal`](inc/damp/design/minreal.hpp#L311) | Minimal realization — cancel uncontrollable and unobservable modes |
| [`mpc`](inc/damp/controllers/mpc.hpp#L1609) | Synthesize a constrained linear MPC (condensed dense QP, Δu form) on a discrete plant (+3 more overloads) |
| [`MPC`](inc/damp/controllers/mpc.hpp#L2466) | MPC law{res}: the sizes, scalar and terminal capacity of design::mpc's result (+1 more overload) |
| [`observability_gramian`](inc/damp/design/stability.hpp#L410) | Continuous/discrete observability Gramian W_o |
| [`observability_matrix`](inc/damp/design/stability.hpp#L85) | Compute the observability matrix [C; CA; CA²; ...; CA^(N-1)] |
| [`oflqg`](inc/damp/controllers/oflqg.hpp#L64) | The offset-free linear-quadratic Gaussian law: a Kalman filter on [x; d] and the steady-state target that holds the controlled outputs at r under d̂: runs as OFLQG (+1 more overload) |
| [`OFLQG`](inc/damp/controllers/oflqg.hpp#L152) | OFLQG law{res}: the sizes, scalar, filter and target of design::oflqg's result |
| [`ofmpc`](inc/damp/controllers/ofmpc.hpp#L43) | The offset-free MPC on a discrete plant: a Kalman filter on [x; d] (noise, disturbance) and the QP posed about the steady-state target under d̂: runs as OFMPC (+1 more overload) |
| [`OFMPC`](inc/damp/controllers/ofmpc.hpp#L138) | OFMPC law{res}: the sizes, scalar, filter and target of design::ofmpc's result (+1 more overload) |
| [`option`](inc/damp/design/options.hpp#L58) | No option left: the fallback (+1 more overload) |
| [`output_inverse`](inc/damp/design/model_inverse.hpp#L685) | The stable inverse of a discrete model: the state and input of its trajectory through a reference y on the controlled outputs, from y ₖ … y ₖ₊NP₋₁ (+1 more overload) |
| [`phase_margin_from_damping_ratio`](inc/damp/design/pid_design.hpp#L888) | Approximate phase margin from damping ratio |
| [`pid`](inc/damp/controllers/pid.hpp#L363) | Continuous PID gains Kₚ, Kᵢ, K_d with every other knob at its default |
| [`pid_pole_placement`](inc/damp/design/pid_design.hpp#L960) | PI gains that place the closed-loop poles of a first-order plant b₀/(a₁ s + a₀) (+4 more overloads) |
| [`PIDController`](inc/damp/controllers/pid.hpp#L904) | Deduce the runtime from its gains, keeping their scalar: PIDController pid{gains} is a PID; through an alias the mode comes from the alias: PIController pi{gains}; |
| [`pidtune`](inc/damp/design/pid_design.hpp#L807) | Plant-aware PID tune at a crossover |
| [`place`](inc/damp/design/pole_placement.hpp#L1164) | Robust pole placement (KNV): the gain K with eig(A − B·K) at poles, in the domain of sys (+1 more overload) |
| [`place_jordan`](inc/damp/design/pole_placement.hpp#L1685) | Exact pole placement with an arbitrary Jordan structure (Schmid–Ntogramatzidis–Nguyen–Pandey / Klein–Moore parametric form) (+1 more overload) |
| [`pr`](inc/damp/controllers/pr.hpp#L205) | Design a Proportional-Resonant controller |
| [`pr_harmonics`](inc/damp/controllers/pr.hpp#L243) | Design multiple-harmonic PR controller gains |
| [`PRController`](inc/damp/controllers/pr.hpp#L438) | Deduce the runtime from its design, keeping the result's scalar |
| [`project_affine`](inc/damp/controllers/action_governor.hpp#L137) | Project u_des onto the polyhedron A u ≤ b (Euclidean) |
| [`project_box`](inc/damp/controllers/action_governor.hpp#L76) | Euclidean projection of u onto a Bounds box |
| [`raise_loop_for_tracking`](inc/damp/design/pid_design.hpp#L713) | Raise the loop crossover so T −3 dB matches the 1-DOF target |
| [`rank`](inc/damp/design/stability.hpp#L425) | Numerical rank of a matrix (forwards to mat::rank) (+1 more overload) |
| [`repetitive`](inc/damp/controllers/repetitive.hpp#L214) | Synthesize a repetitive controller with a scalar robustness filter Q |
| [`repetitive_binomial`](inc/damp/controllers/repetitive.hpp#L247) | Synthesize a repetitive controller with a binomial zero-phase FIR Q |
| [`simc`](inc/damp/design/pid_design.hpp#L368) | SIMC (Skogestad Internal Model Control) tuning for FOPDT models |
| [`smc`](inc/damp/controllers/smc.hpp#L90) | Bundle hand-picked SMC parameters into an SMCResult |
| [`SMCController`](inc/damp/controllers/smc.hpp#L268) | Deduce the runtime from its design, keeping the result's scalar |
| [`smith_predictor_from_fopdt`](inc/damp/controllers/smith_predictor.hpp#L254) | Design a Smith predictor from FOPDT parameters |
| [`smith_predictor_from_pid`](inc/damp/controllers/smith_predictor.hpp#L287) | Pack a Smith predictor from an existing continuous PID and FOPDT plant |
| [`solve_miqp`](inc/damp/design/qp.hpp#L1499) | Mixed-integer QP by branch and bound over the active-set relaxation |
| [`solve_qp`](inc/damp/design/qp.hpp#L610) | Solve a strictly convex inequality-constrained QP (dual active-set) (+1 more overload) |
| [`solve_qp_admm`](inc/damp/design/qp.hpp#L1174) | One-shot ADMM QP solve (cold start) |
| [`solve_qp_interior_point`](inc/damp/design/qp.hpp#L1218) | Interior-point QP solver (primal-dual path following) |
| [`stability_margin`](inc/damp/design/stability.hpp#L527) | Distance from the poles of sys to its domain's stability boundary |
| [`state_feedforward`](inc/damp/design/model_inverse.hpp#L98) | The plant-inversion feedforward of a discrete model: the input that moves x⁺ = A·x + B·u along a state reference, u ₖ the least-squares solution of B·u = x ₖ₊₁ − A·x ₖ |
| [`steady_state_target`](inc/damp/design/steady_state_target.hpp#L328) | The steady-state target of a law on a discrete plant: the maps from the estimated disturbance and the reference to [x_t; u_t] (+1 more overload) |
| [`stsmc`](inc/damp/controllers/stsmc.hpp#L99) | Synthesize super-twisting gains from a disturbance-derivative bound |
| [`stsmc_gains`](inc/damp/controllers/stsmc.hpp#L143) | Super-twisting controller from gains you specify directly |
| [`STSMCController`](inc/damp/controllers/stsmc.hpp#L321) | Deduce the runtime from its design, keeping the result's scalar |
| [`suggest_mpc_horizon`](inc/damp/controllers/mpc.hpp#L1722) | Suggest MPC horizons from an explicit settling-time target (+1 more overload) |
| [`tyreus_luyben`](inc/damp/design/pid_design.hpp#L212) | Tyreus-Luyben tuning from ultimate gain and ultimate period |
| [`unbounded_bound`](inc/damp/design/qp.hpp#L69) | Sentinel bound treated as "no constraint" on that row |
| [`valid_spec`](inc/damp/design/pid_design.hpp#L477) | The spec's scalar domain: ωc > 0, 0 < PM < 180°, Tf ≥ 0, all finite |
| [`ziegler_nichols`](inc/damp/design/pid_design.hpp#L110) | Ziegler-Nichols tuning from ultimate gain and ultimate period |
| [`ziegler_nichols_step`](inc/damp/design/pid_design.hpp#L161) | Ziegler-Nichols step response method (reaction curve) |

> **Runtime — tick objects and signal processing**

## Runtime controllers

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`ActionGovernor`](inc/damp/controllers/action_governor.hpp#L237) | Runtime action governor — QP projection of u_des onto A u ≤ b |
| [`ADRCController`](inc/damp/controllers/adrc.hpp#L257) | Active Disturbance Rejection Control (ADRC) |
| [`Cascade`](inc/damp/controllers/composition.hpp#L50) | Series cascade of two SISO controllers: outer → inner reference |
| [`CBFConstraint`](inc/damp/controllers/action_governor.hpp#L189) | One affine row from a relative-degree-1 CBF condition |
| [`ContinuousPID`](inc/damp/controllers/pid.hpp#L941) | Continuous-gain PID with per-tick sample time (variable-rate secondary form) |
| [`DiscretePIDGains`](inc/damp/controllers/pid.hpp#L76) | Fixed-rate discrete PID coefficients (canonical deploy form) |
| [`ESCSpec`](inc/damp/controllers/esc.hpp#L71) | The extremum-seeking design input: dither, gain, filter corners and clamp; Ts is passed beside it |
| [`ExtremumSeekingController`](inc/damp/controllers/esc.hpp#L236) | Extremum-seeking controller runtime (model-free online optimizer) |
| [`ExtremumType`](inc/damp/controllers/esc.hpp#L49) | Whether ESC climbs to a maximum or descends to a minimum of the objective |
| [`HarmonicSuppressor`](inc/damp/controllers/harmonic_suppression.hpp#L166) | Multi-resonant harmonic suppressor — a parallel bank of PR resonators |
| [`IntegralAction`](inc/damp/controllers/state_feedback_core.hpp#L583) | Offset-free by integral action: integrators on r − y, every measured output (`IntegralAction{}`) or the listed ones (`IntegralAction<2>{{0, 3}}`, in that order: the references r follow it) |
| [`IntegralDesign`](inc/damp/controllers/state_feedback_core.hpp#L635) | The integral action of a design: the integrated outputs, their maps of the state and the input, and the per-sample step |
| [`Law`](inc/damp/controllers/law_detail.hpp#L203) | The law a design result deploys as, named by the result: `Law<kRes, float> law;` (+1 more overload) |
| [`LeadLagController`](inc/damp/controllers/lead_lag.hpp#L316) | Discrete lead-lag compensator |
| [`LQG`](inc/damp/controllers/lqg.hpp#L142) | Linear-quadratic Gaussian law: LQR on the Kalman filter's estimate |
| [`LQGI`](inc/damp/controllers/lqgi.hpp#L98) | Linear-quadratic Gaussian law with integral action on the measured outputs |
| [`LQI`](inc/damp/controllers/lqi.hpp#L89) | Linear-quadratic regulator with integral action, on the measured state |
| [`LQR`](inc/damp/controllers/lqr.hpp#L110) | Linear-quadratic regulator on the measured state: u = u − K·(x − x ) |
| [`LQRWeights`](inc/damp/controllers/state_feedback_core.hpp#L614) | The cost of a state-feedback design: what to penalize and how much |
| [`MeasuredDisturbance`](inc/damp/controllers/mpc.hpp#L1201) | A measured disturbance: a known input the controller does not command, entering the state update through Bd (a load torque, a feed rate, the ambient temperature) |
| [`MPC`](inc/damp/controllers/mpc.hpp#L2396) | Constrained linear MPC on the measured state (design::mpc) |
| [`MPCArtifacts`](inc/damp/controllers/mpc.hpp#L351) | Condensed-QP data produced by mpc() (MPCResult::mpc), consumed by damp::MPC |
| [`MPCConstraints`](inc/damp/controllers/mpc.hpp#L294) | Box constraints for mpc() |
| [`MPCHorizonSuggestion`](inc/damp/controllers/mpc.hpp#L1698) | Advisory NP/NC horizon values from suggest_mpc_horizon() |
| [`MPCLimits`](inc/damp/controllers/mpc.hpp#L328) | The limits an MPC takes at run time (MATLAB®'s online constraints): u, Δu and the controlled outputs, in physical units. A field left out is unbounded |
| [`MPCTerminal`](inc/damp/controllers/mpc.hpp#L230) | Opt-in terminal cost / set for `design::mpc` (default: none) |
| [`MPCWeights`](inc/damp/controllers/mpc.hpp#L256) | Cost weights for mpc() |
| [`NoEstimator`](inc/damp/controllers/state_feedback_core.hpp#L53) | The estimator of a law that measures its state: none, an empty type, so a `[[no_unique_address]]` member of it takes no storage |
| [`OFLQG`](inc/damp/controllers/oflqg.hpp#L102) | Offset-free linear-quadratic Gaussian law: a Kalman filter on [x; d] and the steady-state target |
| [`OFMPC`](inc/damp/controllers/ofmpc.hpp#L83) | Offset-free constrained linear MPC: a Kalman filter on [x; d] and the QP about the steady-state target under d̂ (design::ofmpc) |
| [`PIDController`](inc/damp/controllers/pid.hpp#L515) | Fixed-rate discrete 2-DOF PID (canonical runtime) (+3 more overloads) |
| [`PIDGains`](inc/damp/controllers/pid.hpp#L137) | 2-DOF continuous-time PID gains, limits and setpoint weights |
| [`PIDMode`](inc/damp/controllers/pid.hpp#L376) | Compile-time selection of the PID control-law structure |
| [`PIDRuntimeMode`](inc/damp/controllers/pid.hpp#L401) | Runtime operating mode for PIDController / ContinuousPID |
| [`PRController`](inc/damp/controllers/pr.hpp#L291) | Discrete Proportional-Resonant Controller |
| [`PRSpec`](inc/damp/controllers/pr.hpp#L98) | The gains and frequencies of a proportional-resonant design; Ts is passed beside it |
| [`RepetitiveConfig`](inc/damp/controllers/repetitive.hpp#L77) | Repetitive-controller tuning + period (with optional zero-phase FIR Q) |
| [`RepetitiveController`](inc/damp/controllers/repetitive.hpp#L291) | Plug-in repetitive controller runtime (fixed-size internal model) |
| [`SMCController`](inc/damp/controllers/smc.hpp#L182) | First-order sliding-mode controller (SMC) for a SISO plant |
| [`SmithPredictor`](inc/damp/controllers/smith_predictor.hpp#L329) | Smith predictor runtime: SISO primary + FO model + pure delay |
| [`STSMCController`](inc/damp/controllers/stsmc.hpp#L205) | Super-twisting controller (second-order sliding mode) |
| [`SwitchedController`](inc/damp/controllers/composition.hpp#L115) | Bumpless-ish switch between normal, experiment, and backup SISO laws |
| [`SwitchMode`](inc/damp/controllers/composition.hpp#L92) | Which path owns the plant command |

## Observers & estimators

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AdaptationGate`](inc/damp/estimation/parameter_estimation.hpp#L339) | Gated-adaptation policy: decide whether a parameter update is safe to apply |
| [`Chirp`](inc/damp/estimation/excitation/chirp.hpp#L141) | Linear or logarithmic chirp runtime generator (+1 more overload) |
| [`ChirpConfig`](inc/damp/estimation/excitation/chirp.hpp#L41) | Configuration for a sine chirp excitation |
| [`ChirpMode`](inc/damp/estimation/excitation/chirp.hpp#L27) | Chirp sweep law |
| [`ChirpResult`](inc/damp/estimation/excitation/chirp.hpp#L84) | Chirp design payload |
| [`ClassicalDOB`](inc/damp/estimation/dob.hpp#L439) | Classical Pn^-1·Q disturbance observer runtime (bolt-on compensator) (+1 more overload) |
| [`ClassicalDobResult`](inc/damp/estimation/dob.hpp#L339) | Design result for the classical Pn^-1·Q disturbance observer |
| [`ComplementaryFilter`](inc/damp/estimation/sensor_fusion.hpp#L61) | Simple complementary filter for orientation estimation |
| [`DOB`](inc/damp/estimation/dob.hpp#L168) | Lightweight SISO disturbance observer runtime (+1 more overload) |
| [`DOBConfig`](inc/damp/estimation/dob.hpp#L32) | Configuration for a first-order disturbance observer |
| [`ErrorStateJacobian`](inc/damp/estimation/eskf.hpp#L287) | ESKF prediction Jacobians F and G (nominal state updated outside) |
| [`ErrorStateKalmanFilter`](inc/damp/estimation/eskf.hpp#L341) | Runtime error-state KF: tracks δx and P, not the full nominal state (+1 more overload) |
| [`ErrorStateKalmanFilterConfig`](inc/damp/estimation/eskf.hpp#L312) | Configuration of an ErrorStateKalmanFilter: the process noise and the initial covariance |
| [`ESKFOrientationFilter`](inc/damp/estimation/sensor_fusion.hpp#L380) | Turnkey 6-state attitude ESKF: owns q, b_g, and the error filter |
| [`ESKFPredictFn`](inc/damp/estimation/eskf.hpp#L299) | Callable that returns ErrorStateJacobian for one predict step |
| [`ESKFResult`](inc/damp/estimation/eskf.hpp#L85) | Design payload for an ESKF: Q, R, P₀, and success |
| [`ESO`](inc/damp/estimation/eso.hpp#L273) | Extended-state observer (ESO) (+1 more overload) |
| [`ESOResult`](inc/damp/estimation/eso.hpp#L55) | Extended-state observer design: continuous β and, for Ts > 0, the discrete ESO |
| [`ExtendedKalmanFilter`](inc/damp/estimation/ekf.hpp#L118) | Extended Kalman filter: continuous nonlinear dynamics, sampled measurements (+1 more overload) |
| [`ExtendedKalmanFilterConfig`](inc/damp/estimation/ekf.hpp#L42) | Configuration of an ExtendedKalmanFilter: its sample time, the process noise and the initial estimate |
| [`Filtered`](inc/damp/estimation/state_variable_filter.hpp#L67) | One signal after the state-variable filter: its filtered value and first two derivatives |
| [`FirstOrderPlantEstimator`](inc/damp/estimation/parameter_estimation.hpp#L126) | Online estimator for a first-order plant's gain and time constant |
| [`FirstOrderPlantEstimatorConfig`](inc/damp/estimation/parameter_estimation.hpp#L85) | Configuration for the first-order grey-box estimator |
| [`FrequencyResponseEstimator`](inc/damp/estimation/frequency_response.hpp#L168) | Lock-in FRF estimator over a fixed table of NFreq bins |
| [`FrfFeatureKind`](inc/damp/estimation/frequency_response.hpp#L437) | Peak (resonance) or valley (anti-resonance) on an FRF magnitude curve |
| [`FrfFeatureResult`](inc/damp/estimation/frequency_response.hpp#L520) | Combined peak + valley reduction of an FRF table |
| [`FrfMargins`](inc/damp/estimation/frequency_response.hpp#L555) | Open-loop margins read from a discrete FRF table (on-target) |
| [`FrfPoint`](inc/damp/estimation/frequency_response.hpp#L57) | One measured frequency-response point (on-target FRF table entry) |
| [`Impulse`](inc/damp/estimation/excitation/impulse.hpp#L89) | Rectangular impulse (finite-width Dirac stand-in) |
| [`ImpulseConfig`](inc/damp/estimation/excitation/impulse.hpp#L32) | Configuration for a one-shot rectangular impulse (Dirac stand-in) |
| [`ImpulseResult`](inc/damp/estimation/excitation/impulse.hpp#L52) | Impulse design payload |
| [`ImuSample`](inc/damp/estimation/ins_mechanization.hpp#L125) | IMU sample in the body frame |
| [`InsNavigator`](inc/damp/estimation/ins_eskf.hpp#L695) | Runtime strapdown navigator: nominal InsState + 15-state ESKF (+2 more overloads) |
| [`InsState`](inc/damp/estimation/ins_mechanization.hpp#L139) | Nominal strapdown navigation state (external to the ESKF error state) |
| [`KalmanFilter`](inc/damp/estimation/kalman.hpp#L667) | The recursive Kalman filter: the covariance propagates each tick, and the gain is formed from it in Joseph form (+2 more overloads) |
| [`KalmanFilterConfig`](inc/damp/estimation/kalman.hpp#L618) | Configuration of a KalmanFilter built from a model rather than from design::kalman |
| [`KalmanNoise`](inc/damp/estimation/kalman.hpp#L80) | The estimator's noise: process noise on the plant's noise inputs (through G), measurement noise on the sensors (through H) (+2 more overloads) |
| [`KalmanResult`](inc/damp/estimation/kalman.hpp#L138) | Steady-state Kalman filter design result |
| [`Luenberger`](inc/damp/estimation/luenberger.hpp#L583) | The fixed-gain state observer: the innovation gain M held, the current form, and the optional disturbance anchor (+4 more overloads) |
| [`LuenbergerResult`](inc/damp/estimation/luenberger.hpp#L117) | Luenberger observer design result: the plant, the gains L and M, and the achieved error poles |
| [`MadgwickFilter`](inc/damp/estimation/sensor_fusion.hpp#L135) | Madgwick gradient-descent AHRS filter |
| [`MahonyFilter`](inc/damp/estimation/sensor_fusion.hpp#L177) | Mahony nonlinear complementary filter with PI correction |
| [`MHE`](inc/damp/estimation/mhe.hpp#L346) | Runtime moving-horizon estimator (fixed per-tick iteration budget) (+1 more overload) |
| [`MHEConstraints`](inc/damp/estimation/mhe.hpp#L97) | Box constraints on the estimated states (applied to every window state) and the QP budget |
| [`MHEResult`](inc/damp/estimation/mhe.hpp#L117) | Window QP data produced by mhe(), consumed by damp::MHE |
| [`ModeExtractorConfig`](inc/damp/estimation/frequency_response.hpp#L475) | Configuration for FRF peak / valley extraction |
| [`ModeExtractorResult`](inc/damp/estimation/frequency_response.hpp#L492) | Result of reducing an FRF table to at most MaxModes peaks (resonances) |
| [`MultiSine`](inc/damp/estimation/excitation/multi_sine.hpp#L213) | Sum-of-tones multi-sine runtime generator (+1 more overload) |
| [`MultiSineConfig`](inc/damp/estimation/excitation/multi_sine.hpp#L62) | Configuration for fixed-component multi-sine excitation |
| [`MultiSineResult`](inc/damp/estimation/excitation/multi_sine.hpp#L104) | Multi-sine design payload |
| [`NavFrame`](inc/damp/estimation/ins_mechanization.hpp#L77) | Local-level navigation frame for strapdown mechanization |
| [`ParameterDriftMonitor`](inc/damp/estimation/parameter_estimation.hpp#L252) | Residual-based drift detector: has the plant moved away from the model? |
| [`PRBS`](inc/damp/estimation/excitation/prbs.hpp#L189) | Maximal-length PRBS runtime generator (+1 more overload) |
| [`PRBSConfig`](inc/damp/estimation/excitation/prbs.hpp#L82) | Configuration for maximal-length pseudo-random binary excitation |
| [`PRBSResult`](inc/damp/estimation/excitation/prbs.hpp#L126) | PRBS design payload |
| [`Ramp`](inc/damp/estimation/excitation/ramp.hpp#L117) | Rate-limited ramp runtime generator (+1 more overload) |
| [`RampConfig`](inc/damp/estimation/excitation/ramp.hpp#L33) | Configuration for a slew-rate-limited ramp excitation |
| [`RampResult`](inc/damp/estimation/excitation/ramp.hpp#L65) | Ramp design payload |
| [`ReducedLuenberger`](inc/damp/estimation/luenberger.hpp#L866) | Reduced-order (Gopinath) state observer (runtime) (+1 more overload) |
| [`ReducedLuenbergerResult`](inc/damp/estimation/luenberger.hpp#L342) | Reduced-order (Gopinath) observer design result |
| [`RelayAutotuneConfig`](inc/damp/estimation/relay_autotune.hpp#L128) | Configuration for the relay-feedback autotuning experiment |
| [`RelayAutotuneOutput`](inc/damp/estimation/relay_autotune.hpp#L260) | Per-tick output of RelayAutotuner::step |
| [`RelayAutotuner`](inc/damp/estimation/relay_autotune.hpp#L279) | Runtime relay-feedback autotuner (+1 more overload) |
| [`RelayAutotuneResult`](inc/damp/estimation/relay_autotune.hpp#L187) | Relay-autotuner design payload |
| [`RelayAutotuneStatus`](inc/damp/estimation/relay_autotune.hpp#L244) | Lifecycle state of a RelayAutotuner |
| [`ResonantMode`](inc/damp/estimation/frequency_response.hpp#L448) | One extracted FRF extremum (peak or valley) |
| [`RlsConfig`](inc/damp/estimation/rls.hpp#L65) | Common RLS configuration |
| [`RlsVector`](inc/damp/estimation/rls.hpp#L150) | Vector runtime RLS estimator (NP parameters) (+1 more overload) |
| [`RlsVectorResult`](inc/damp/estimation/rls.hpp#L106) | Vector RLS design payload for N parameters |
| [`RlsVectorState`](inc/damp/estimation/rls.hpp#L92) | Vector RLS runtime state for N parameters |
| [`StateVariableFilter`](inc/damp/estimation/state_variable_filter.hpp#L129) | State-variable filter over NS signals: filtered values and derivatives, no differentiation |
| [`StateVariableFilterResult`](inc/damp/estimation/state_variable_filter.hpp#L78) | Discretized state-variable filter (shared by every channel) |
| [`Step`](inc/damp/estimation/excitation/step.hpp#L117) | Single step with optional \|dy/dt\| settle detection (+1 more overload) |
| [`StepConfig`](inc/damp/estimation/excitation/step.hpp#L36) | Configuration for a single step with optional settle detection |
| [`SteppedSine`](inc/damp/estimation/excitation/stepped_sine.hpp#L136) | Stepped-sine excitation — one pure tone at a time across a frequency table (+1 more overload) |
| [`SteppedSineConfig`](inc/damp/estimation/excitation/stepped_sine.hpp#L34) | Configuration for single-tone stepped-sine excitation |
| [`SteppedSineResult`](inc/damp/estimation/excitation/stepped_sine.hpp#L74) | Stepped-sine design payload |
| [`StepResult`](inc/damp/estimation/excitation/step.hpp#L68) | Single-step design payload |
| [`StepTrain`](inc/damp/estimation/excitation/step_train.hpp#L116) | Alternating +/- step train runtime generator (+1 more overload) |
| [`StepTrainConfig`](inc/damp/estimation/excitation/step_train.hpp#L32) | Configuration for alternating +/- step excitation |
| [`StepTrainResult`](inc/damp/estimation/excitation/step_train.hpp#L64) | Step-train design payload |
| [`Tone`](inc/damp/estimation/excitation/multi_sine.hpp#L29) | One sinusoidal component in a multi-sine excitation |
| [`UnscentedKalmanFilter`](inc/damp/estimation/ukf.hpp#L133) | Unscented (sigma-point) Kalman Filter: continuous nonlinear dynamics, sampled measurements (+1 more overload) |
| [`UnscentedKalmanFilterConfig`](inc/damp/estimation/ukf.hpp#L82) | Configuration of an UnscentedKalmanFilter: its sample time, the process noise, the initial estimate and the unscented transform's tuning |
| [`UnscentedParams`](inc/damp/estimation/ukf.hpp#L60) | Tuning parameters for the scaled unscented transform |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`chirp`](inc/damp/estimation/excitation/chirp.hpp#L123) | Build a chirp design payload from a configuration |
| [`classical_dob`](inc/damp/estimation/dob.hpp#L404) | Synthesize a classical disturbance observer from a nominal plant and Q-filter |
| [`discretize_process_noise`](inc/damp/estimation/kalman.hpp#L269) | Discretize continuous process noise over one sample (Van Loan) |
| [`dob`](inc/damp/estimation/dob.hpp#L143) | Validate and package DOB configuration into a runtime-ready design result |
| [`enu_from_ned`](inc/damp/estimation/ins_mechanization.hpp#L116) | Map a NED vector into ENU (axis permute) |
| [`eskf_imu`](inc/damp/estimation/eskf.hpp#L219) | Design Q, R, P₀ for a 6-state IMU attitude ESKF (gyro + accel) |
| [`eskf_marg`](inc/damp/estimation/eskf.hpp#L258) | Design Q, R, P₀ for a 6-state MARG attitude ESKF (gyro + accel + mag) |
| [`eso`](inc/damp/estimation/eso.hpp#L207) | Extended-state observer: poles at −ω_o, plant gain b₀, optional sample time |
| [`extract_frf_features`](inc/damp/estimation/frequency_response.hpp#L810) | Extract peaks and valleys in one pass over the FRF table |
| [`extract_modes`](inc/damp/estimation/frequency_response.hpp#L785) | Extract resonant peaks (local \|G\| maxima) from an FRF table |
| [`extract_valleys`](inc/damp/estimation/frequency_response.hpp#L799) | Extract anti-resonance valleys (local \|G\| minima) from an FRF table |
| [`finish_extremum`](inc/damp/estimation/frequency_response.hpp#L631) | Half-power (peaks) or double-power (valleys) bandwidth → ζ |
| [`gravity_nav`](inc/damp/estimation/ins_mechanization.hpp#L93) | Gravity vector in the local-level frame [m/s²] |
| [`heading_from_baseline_nav`](inc/damp/estimation/ins_eskf.hpp#L491) | Heading [rad] from a nav-frame baseline vector (dual-antenna difference) |
| [`impulse`](inc/damp/estimation/excitation/impulse.hpp#L73) | Build an impulse design payload |
| [`ins_aid_position`](inc/damp/estimation/ins_eskf.hpp#L458) | Predicted antenna / marker position with body lever-arm |
| [`ins_apply_correction`](inc/damp/estimation/ins_eskf.hpp#L444) | Inject filter correction into x and reset the error state |
| [`ins_error_jacobian`](inc/damp/estimation/ins_eskf.hpp#L262) | Discrete error-state transition and noise input for one IMU step |
| [`ins_eskf_design`](inc/damp/estimation/ins_eskf.hpp#L174) | Build Q / R / P0 for a 15-state INS ESKF from IMU noise densities |
| [`ins_heading`](inc/damp/estimation/ins_eskf.hpp#L476) | Horizontal heading [rad] of a body axis expressed in the nav frame |
| [`ins_inject`](inc/damp/estimation/ins_eskf.hpp#L392) | Inject δx into the nominal INS state (right-multiplicative q) |
| [`ins_predict`](inc/damp/estimation/ins_eskf.hpp#L421) | Mechanize nominal state and propagate the error covariance |
| [`ins_propagate_covariance`](inc/damp/estimation/ins_eskf.hpp#L358) | P ← F P Fᵀ + Q using sparse F (exact for ins_error_jacobian structure) |
| [`ins_update_heading`](inc/damp/estimation/ins_eskf.hpp#L579) | Sparse dual-antenna / yaw heading update (scalar) |
| [`ins_update_orientation`](inc/damp/estimation/ins_eskf.hpp#L617) | Sparse orientation update from a measured body→nav quaternion |
| [`ins_update_pose`](inc/damp/estimation/ins_eskf.hpp#L642) | Generic pose alias: position then orientation (two sparse updates) |
| [`ins_update_pose_heading`](inc/damp/estimation/ins_eskf.hpp#L661) | Pose alias: position + dual-antenna heading (two sparse updates) |
| [`ins_update_position`](inc/damp/estimation/ins_eskf.hpp#L506) | Sparse position update: y = p + Rℓ |
| [`ins_update_velocity`](inc/damp/estimation/ins_eskf.hpp#L538) | Sparse velocity update: y = v at the IMU origin (no lever-arm) |
| [`ins_update_zupt`](inc/damp/estimation/ins_eskf.hpp#L557) | Zero-velocity update (ZUPT): velocity measurement of zero |
| [`kalman`](inc/damp/estimation/kalman.hpp#L485) | Steady-state Kalman filter design (+3 more overloads) |
| [`luenberger`](inc/damp/estimation/luenberger.hpp#L294) | Design a Luenberger observer by robust pole placement, in the domain of sys (+1 more overload) |
| [`margins_from_frf`](inc/damp/estimation/frequency_response.hpp#L842) | Gain / phase margins from a fixed on-target FRF table |
| [`mechanize_step`](inc/damp/estimation/ins_mechanization.hpp#L210) | One dead-reckoning step from a raw IMU sample (subtracts the state's biases) |
| [`mechanize_step_from_corrected`](inc/damp/estimation/ins_mechanization.hpp#L175) | One dead-reckoning step from an IMU sample (strapdown mechanization) |
| [`mhe`](inc/damp/estimation/mhe.hpp#L227) | Synthesize a constrained moving-horizon estimator |
| [`multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L139) | Build a multi-sine design payload from a configuration |
| [`ned_from_enu`](inc/damp/estimation/ins_mechanization.hpp#L106) | Map an ENU vector into NED (axis permute) |
| [`prbs`](inc/damp/estimation/excitation/prbs.hpp#L166) | Build a PRBS design payload from a configuration |
| [`ramp`](inc/damp/estimation/excitation/ramp.hpp#L101) | Build a ramp design payload from a configuration |
| [`reduced_luenberger`](inc/damp/estimation/luenberger.hpp#L501) | Design a reduced-order (Gopinath) observer by pole placement, in the domain of sys (+1 more overload) |
| [`relay_autotune`](inc/damp/estimation/relay_autotune.hpp#L228) | Build a validated relay-autotune design payload |
| [`rls_vector`](inc/damp/estimation/rls.hpp#L135) | Build vector RLS design payload |
| [`rows`](inc/damp/estimation/ins_eskf.hpp#L303) | Sparse y = F x for the 15-state INS first-order structure (G = I path) |
| [`schroeder_multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L193) | Build a multi-sine design payload with Schroeder phases applied |
| [`specific_force_at_rest`](inc/damp/estimation/ins_mechanization.hpp#L227) | Specific force [m/s²] a stationary IMU measures for the given orientation |
| [`state_variable_filter`](inc/damp/estimation/state_variable_filter.hpp#L99) | Design the filter F(s) = ω³/(s + ω)³ sampled at Ts |
| [`step`](inc/damp/estimation/excitation/step.hpp#L100) | Build a single-step design payload |
| [`step_train`](inc/damp/estimation/excitation/step_train.hpp#L100) | Build a step-train design payload from a configuration |
| [`stepped_sine`](inc/damp/estimation/excitation/stepped_sine.hpp#L111) | Build a stepped-sine design payload from a configuration |
| [`with_schroeder_phases`](inc/damp/estimation/excitation/multi_sine.hpp#L170) | Assign Schroeder low-crest-factor phases to a multi-sine tone table |

## Filters & signal conditioning

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`Biquad`](inc/damp/filters/biquad.hpp#L538) | Second-order IIR (biquad) section runtime (+1 more overload) |
| [`BiquadCascade`](inc/damp/filters/biquad.hpp#L628) | Cascade of second-order sections (SOS) for higher-order IIR filters (+1 more overload) |
| [`Complementary`](inc/damp/filters/complementary.hpp#L32) | Scalar (1-D) complementary filter — fuse a fast rate with a slow absolute |
| [`Delay`](inc/damp/filters/delay.hpp#L27) | Discrete-time delay buffer |
| [`DsogiPll`](inc/damp/filters/pll.hpp#L306) | Dual-SOGI three-phase positive-sequence PLL (DSOGI-PLL) |
| [`FirDesignResult`](inc/damp/filters/fir.hpp#L67) | Window-method FIR design result |
| [`FirFilter`](inc/damp/filters/fir.hpp#L465) | Direct-form FIR filter runtime (tapped delay line + dot product) (+1 more overload) |
| [`FirstOrderCoeffs`](inc/damp/filters/iir_design.hpp#L31) | DSP coefficients for first-order IIR filter |
| [`FirstOrderSection`](inc/damp/filters/first_order_section.hpp#L44) | First-order IIR section runtime, Direct Form I (+1 more overload) |
| [`FirType`](inc/damp/filters/fir.hpp#L38) | Ideal frequency-selective response for fir1 / fir_window |
| [`FirWindow`](inc/damp/filters/fir.hpp#L50) | Window applied to the ideal truncated impulse response |
| [`Goertzel`](inc/damp/filters/spectral.hpp#L52) | Generalized Goertzel single-bin DFT — amplitude/phase at one frequency |
| [`HarmonicAnalyzer`](inc/damp/filters/spectral.hpp#L204) | Harmonic analyzer — a Goertzel bank over a fundamental and K−1 harmonics |
| [`LowPass`](inc/damp/filters/lowpass.hpp#L223) | Nth-order low-pass filter (direct form I) (+1 more overload) |
| [`LowPassResult`](inc/damp/filters/lowpass.hpp#L58) | Direct-form coefficients of an Nth-order low-pass, with success / status |
| [`MedianFilter`](inc/damp/filters/median.hpp#L34) | Sliding-window median filter — nonlinear spike/outlier rejection |
| [`MovingAverage`](inc/damp/filters/moving_average.hpp#L67) | Moving-average (boxcar) filter — also a DC-preserving harmonic-notch comb |
| [`MSTOGI`](inc/damp/filters/sogi.hpp#L299) | Runtime MSTOGI with exact resonator and forward-Euler washout |
| [`NlmsFilter`](inc/damp/filters/fir.hpp#L575) | Normalized LMS adaptive FIR filter |
| [`RobustExactDifferentiator`](inc/damp/filters/differentiator.hpp#L62) | First-order robust exact differentiator (super-twisting differentiator) |
| [`SecondOrderCoeffs`](inc/damp/filters/iir_design.hpp#L72) | DSP coefficients for second-order IIR filter |
| [`SecondOrderSectionsResult`](inc/damp/filters/biquad.hpp#L42) | An IIR filter as a cascade of second-order sections: the result BiquadCascade runs |
| [`SinglePhasePLL`](inc/damp/filters/pll.hpp#L60) | Single-Phase PLL |
| [`SOGI`](inc/damp/filters/sogi.hpp#L226) | Runtime SOGI quadrature-signal generator (live-frequency resonator) |
| [`SogiFll`](inc/damp/filters/sogi.hpp#L396) | SOGI with a Frequency-Locked Loop — self-tuning single-tone tracker |
| [`ThreePhasePLL`](inc/damp/filters/pll.hpp#L160) | Synchronous-reference-frame (SRF) PLL for balanced three-phase input |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`allpass`](inc/damp/filters/iir_design.hpp#L718) | Second-order all-pass filter design (continuous-time) (+1 more overload) |
| [`bandpass`](inc/damp/filters/iir_design.hpp#L609) | Second-order band-pass filter design (continuous-time) (+1 more overload) |
| [`butterworth_lowpass`](inc/damp/filters/iir_design.hpp#L776) | Butterworth low-pass filter design (+1 more overload) |
| [`comb_notch_window`](inc/damp/filters/moving_average.hpp#L34) | Window length for a moving-average comb that notches f_notch and all its harmonics: N = round(1 / (f_notch·Ts)) |
| [`continuous_lpf_exact_step`](inc/damp/filters/lowpass.hpp#L39) | Exact step of continuous first-order LPF ẏ = −ω_c (y − u) |
| [`fir1`](inc/damp/filters/fir.hpp#L246) | Window-method FIR design (normalized frequency, single cutoff) (+1 more overload) |
| [`fir1_hz`](inc/damp/filters/fir.hpp#L387) | Window-method FIR design from a cutoff in Hz and the sample time (+1 more overload) |
| [`fir_window`](inc/damp/filters/fir.hpp#L423) | Alias for fir1 (descriptive name) (+1 more overload) |
| [`highpass_1st`](inc/damp/filters/iir_design.hpp#L417) | First-order high-pass filter design (continuous-time) (+1 more overload) |
| [`highpass_2nd`](inc/damp/filters/iir_design.hpp#L553) | Second-order high-pass filter design (continuous-time) (+1 more overload) |
| [`highpass_sections`](inc/damp/filters/biquad.hpp#L500) | Nth-order Butterworth high-pass at cutoff fc as second-order sections, any order |
| [`highshelf`](inc/damp/filters/iir_design.hpp#L905) | High-shelf EQ filter: boost or cut everything above fc |
| [`lowpass`](inc/damp/filters/lowpass.hpp#L100) | Discretize a continuous low-pass TF (Tustin) into direct-form coefficients (+1 more overload) |
| [`lowpass_1st`](inc/damp/filters/iir_design.hpp#L360) | First-order low-pass filter design (continuous-time) (+1 more overload) |
| [`lowpass_2nd`](inc/damp/filters/iir_design.hpp#L474) | Second-order low-pass filter design (continuous-time) (+1 more overload) |
| [`lowpass_sections`](inc/damp/filters/biquad.hpp#L464) | Nth-order Butterworth low-pass at cutoff fc as second-order sections, any order |
| [`lowshelf`](inc/damp/filters/iir_design.hpp#L873) | Low-shelf EQ filter: boost or cut everything below fc |
| [`mstogi`](inc/damp/filters/sogi.hpp#L153) | Mixed Second/Third-Order Generalized Integrator (MSTOGI) (+1 more overload) |
| [`negative_sequence_ab`](inc/damp/filters/pll.hpp#L267) | Instantaneous negative-sequence αβ from a quadrature signal pair |
| [`notch`](inc/damp/filters/iir_design.hpp#L661) | Second-order band-reject (notch) filter design (continuous-time) (+1 more overload) |
| [`pade_delay_1st`](inc/damp/filters/pade.hpp#L33) | First-order Pade approximation of time delay (+1 more overload) |
| [`pade_delay_2nd`](inc/damp/filters/pade.hpp#L85) | Second-order Pade approximation of time delay (+1 more overload) |
| [`peaking`](inc/damp/filters/iir_design.hpp#L848) | Peaking (bell) EQ filter: boost or cut a band around f0 |
| [`positive_sequence_ab`](inc/damp/filters/pll.hpp#L246) | Instantaneous positive-sequence αβ from a quadrature signal pair |
| [`second_order_sections`](inc/damp/filters/biquad.hpp#L374) | Factor a discrete transfer function in z into second-order sections, for BiquadCascade (+1 more overload) |
| [`sogi`](inc/damp/filters/sogi.hpp#L69) | Second-Order Generalized Integrator (SOGI) design (+1 more overload) |
| [`to_coeffs`](inc/damp/filters/iir_design.hpp#L132) | Convert StateSpace system to first-order DSP coefficients (+4 more overloads) |

> **Pose / trajectory helpers (core)**

## Trajectory value types

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`Jet`](inc/damp/trajectory/trajectory_types.hpp#L98) | Truncated jet of a scalar motion sample: d[k] = s⁽ᵏ⁾ |
| [`TrajectoryBoundary`](inc/damp/trajectory/trajectory_types.hpp#L181) | Boundary conditions at one endpoint of a polynomial trajectory: a position and its time derivatives through jerk |
| [`TrajectoryLimits`](inc/damp/trajectory/trajectory_types.hpp#L48) | Asymmetric kinematic limits for a trapezoidal or S-curve motion profile |
| [`TrajectoryState`](inc/damp/trajectory/trajectory_types.hpp#L79) | A point on a motion profile: commanded position, velocity, acceleration |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`factorial`](inc/damp/trajectory/trajectory_types.hpp#L209) | k! (exact for modest k in float/double; used up through nonic BVP / jets) |
| [`falling_factorial`](inc/damp/trajectory/trajectory_types.hpp#L222) | Falling factorial i·(i−1)···(i−k+1) = i! / (i−k)! — the k-th derivative coefficient of tⁱ. Zero when i < k |

## Kinematics / pose

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`Pose`](inc/damp/kinematics/pose.hpp#L61) | Rigid-body pose: a translation and an orientation (unit quaternion) |
| [`Translation3`](inc/damp/kinematics/pose.hpp#L38) | A 3-D translation — a thin Vec3 with domain-named conveniences |

## Motor control pack (if present)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`DqCommand`](inc/damp/motor/foc.hpp#L86) | Clamped dq voltage command from FOController::current_controller |
| [`FOController`](inc/damp/motor/foc.hpp#L101) | dq current regulator: PI + decoupling FF → clamped Vdq |
| [`Modulator`](inc/damp/motor/modulation.hpp#L497) | Thin scheme-holding wrapper for FOC / deploy paths |
| [`PwmScheme`](inc/damp/motor/modulation.hpp#L79) | Carrier-based three-phase VSI modulation scheme |
| [`SvmDuties`](inc/damp/motor/modulation.hpp#L103) | Result of a duty-map: half-bridge duties plus an over-modulation flag |
| [`TorqueSpeedPoint`](inc/damp/motor/spm.hpp#L144) | One sample of the SPM torque–speed envelope (MTPA i_d=0 branch) |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`constant_power_torque_limit`](inc/damp/motor/foc.hpp#L44) | Constant-power torque ceiling Tₘₐₓ = P_rated / \|ω\| |
| [`decoupling_feedforward`](inc/damp/motor/foc.hpp#L72) | Decoupling + back-EMF feedforward at R = 0 (reference currents) |
| [`dpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L269) | DPWM zero-sequence for a chosen discontinuous scheme |
| [`duties_from_phase_voltages`](inc/damp/motor/modulation.hpp#L361) | Map phase voltages + zero-sequence to clamped half-bridge duties |
| [`flux_from_Kv`](inc/damp/motor/spm.hpp#L55) | λ from Kᵥ via torque_constant_from_Kv |
| [`flux_from_torque_constant`](inc/damp/motor/spm.hpp#L37) | λ from datasheet Kₜ (amplitude / peak-per-phase convention) |
| [`iq_from_torque`](inc/damp/motor/spm.hpp#L71) | i_q = T_e / Kₜ for id = 0 |
| [`linear_modulation_voltage`](inc/damp/motor/modulation.hpp#L127) | Peak phase voltage at the linear SVPWM hexagon limit |
| [`max_iq`](inc/damp/motor/spm.hpp#L111) | Max positive i_q for an SPM (i_d = 0) at electrical speed |
| [`max_torque_at_speed`](inc/damp/motor/spm.hpp#L172) | SPM max-torque point at a mechanical speed (i_d = 0 only) |
| [`modulation_duty_cycles`](inc/damp/motor/modulation.hpp#L400) | Carrier-based VSI duty cycles from an αβ voltage command |
| [`modulation_index`](inc/damp/motor/modulation.hpp#L165) | Modulation index relative to the linear SVPWM circle |
| [`modulation_zero_sequence`](inc/damp/motor/modulation.hpp#L327) | Zero-sequence for any PwmScheme |
| [`motor_constant`](inc/damp/motor/spm.hpp#L63) | Motor constant Kₘ = Kₜ / √R [Nm/√W] |
| [`six_step_duty_cycles`](inc/damp/motor/modulation.hpp#L466) | Classical six-step (full-wave) duty pattern from the αβ angle |
| [`six_step_fundamental_voltage`](inc/damp/motor/modulation.hpp#L146) | Six-step (square-wave) fundamental peak phase voltage |
| [`steady_state_vdq`](inc/damp/motor/spm.hpp#L79) | Steady-state Vdq (rotor frame) |
| [`steady_state_voltage_magnitude`](inc/damp/motor/spm.hpp#L92) | \|Vdq\| at a steady-state operating point |
| [`svm_duty_cycles`](inc/damp/motor/modulation.hpp#L443) | Space-vector PWM duty cycles from an αβ voltage command |
| [`svpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L200) | Min-max zero-sequence injection for space-vector PWM |
| [`thipwm_zero_sequence`](inc/damp/motor/modulation.hpp#L229) | Third-harmonic injection (THIPWM) zero-sequence — 1/6 of fundamental |
| [`torque`](inc/damp/motor/spm.hpp#L103) | Electromagnetic torque of an SPM (i_d = 0): T_e = Kₜ i_q |
| [`torque_constant_from_flux`](inc/damp/motor/spm.hpp#L29) | Kₜ = 1½ p λ [Nm/A] (amplitude-invariant) |
| [`torque_constant_from_Kv`](inc/damp/motor/spm.hpp#L46) | Kₜ from hobby Kᵥ [RPM/V] (peak line-to-line / bus-volt sense) |
| [`torque_speed_envelope`](inc/damp/motor/spm.hpp#L203) | Fixed table of SPM torque–speed envelope samples (host / LUT bake) |
| [`voltage_circle_radius`](inc/damp/motor/foc.hpp#L60) | SVPWM voltage-circle radius Vₘₐₓ = m · Vdc / √3 |

> **Embedded helpers — controls-adjacent utilities (not ETL)**

## Embedded helpers (controls-adjacent utilities)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AdaptiveLut1D`](inc/damp/toolbox/lookup.hpp#L428) | Adaptive 1-D lookup: fixed breakpoints, online-updated cell values |
| [`AdaptiveLut2D`](inc/damp/toolbox/lookup.hpp#L501) | Adaptive 2-D lookup: fixed grid, online-updated cells (nearest grid node) |
| [`AdaptOutOfRange`](inc/damp/toolbox/lookup.hpp#L48) | How AdaptiveLut1D / AdaptiveLut2D treat samples outside the breakpoint span |
| [`AffineCal`](inc/damp/toolbox/scaling.hpp#L83) | Affine sensor calibration `y = gain·x + offset` |
| [`AnalogCrossMode`](inc/damp/toolbox/conditioning.hpp#L467) | How two redundant analog channels should relate in engineering units |
| [`AnalogCrossResult`](inc/damp/toolbox/conditioning.hpp#L489) | Output of cross_check_analog |
| [`AnalogCrossStatus`](inc/damp/toolbox/conditioning.hpp#L475) | Result of comparing two redundant analog channels |
| [`AnalogInput`](inc/damp/toolbox/io.hpp#L47) | A single analog input: range/fault check on the raw reading, then affine calibration to engineering units |
| [`AxisInput`](inc/damp/toolbox/io.hpp#L148) | Operator-axis conditioning chain (joystick / RC stick → command) |
| [`AxisPosition`](inc/damp/encoder/axis_position.hpp#L47) | The position of an axis, rotary or linear: an exact count carrying its own scale |
| [`BLINK`](inc/damp/toolbox/iec61131.hpp#L542) | BLINK (free-running square-wave / flasher) |
| [`Bounds`](inc/damp/toolbox/bounds.hpp#L52) | A per-channel closed-interval box constraint |
| [`Button`](inc/damp/toolbox/io.hpp#L195) | Debounced momentary push-button with edge and long-press detection |
| [`ConstantInertiaFeedforward`](inc/damp/toolbox/actuator.hpp#L201) | Per-axis decoupled torque feedforward: `τ = J·a + b·v + τ_c·sign(v) + g` |
| [`Counter`](inc/damp/toolbox/logic.hpp#L273) | Edge-counting up/down counter: increments on each rising edge of up, decrements on each rising edge of down. Returns the running count, saturating at the limits of T (no wrap / signed overflow) |
| [`CTD`](inc/damp/toolbox/iec61131.hpp#L354) | CTD Counter (Count Down) |
| [`CTU`](inc/damp/toolbox/iec61131.hpp#L314) | CTU Counter (Count Up) |
| [`CTUD`](inc/damp/toolbox/iec61131.hpp#L394) | CTUD Counter (Count Up Down) |
| [`Debounce`](inc/damp/toolbox/logic.hpp#L212) | Debounce: the output adopts in only after in differs from the current output continuously for stable_time. Rejects contact bounce and brief glitches. (Not an IEC block — the one everyone hand-rolls.) |
| [`DFF`](inc/damp/toolbox/iec61131.hpp#L451) | D Flip-Flop (edge-triggered data latch) |
| [`DLATCH`](inc/damp/toolbox/iec61131.hpp#L482) | D Latch (level-sensitive / transparent latch) |
| [`DrivePair`](inc/damp/toolbox/io.hpp#L234) | Left/right actuator pair (differential/arcade drive output) |
| [`EncoderGate`](inc/damp/encoder/gate.hpp#L109) | Validates encoder readings before they reach the tracker |
| [`EncoderGateConfig`](inc/damp/encoder/gate.hpp#L43) | Configuration for EncoderGate. Positions in counts, as the tracker's |
| [`EncoderTracker`](inc/damp/encoder/tracker.hpp#L70) | Encoder position (exact) and speed (smoothed) from absolute counts |
| [`EncoderTrackerConfig`](inc/damp/encoder/tracker.hpp#L23) | Configuration for EncoderTracker |
| [`EncoderVerdict`](inc/damp/encoder/gate.hpp#L31) | What EncoderGate made of a reading |
| [`Extrapolation`](inc/damp/toolbox/lookup.hpp#L36) | Out-of-range behaviour for a Lut1D / Lut2D / Lut3D query beyond its breakpoints |
| [`F_TRIG`](inc/damp/toolbox/iec61131.hpp#L154) | F_TRIG (Falling Edge Trigger) |
| [`FallingEdge`](inc/damp/toolbox/logic.hpp#L48) | Falling-edge detector: true on the tick x goes true → false |
| [`HolonomicOutput`](inc/damp/toolbox/io.hpp#L355) | Four-wheel actuator set (mecanum/holonomic drive output), X-config |
| [`Hysteresis`](inc/damp/toolbox/conditioning.hpp#L284) | Hysteresis comparator (Schmitt trigger): bool output with separate on/off thresholds to reject chatter |
| [`Interpolation`](inc/damp/toolbox/lookup.hpp#L42) | In-range interpolant for Lut1D (2-D stays bilinear; use SplineSurface for C¹ grids) |
| [`Latch`](inc/damp/toolbox/logic.hpp#L66) | Set/reset latch. @tparam SetDominant which input wins when both are asserted (default: set-dominant, e.g. a trip overriding a clear) |
| [`Lut1D`](inc/damp/toolbox/lookup.hpp#L109) | 1-D interpolating lookup table over monotonic breakpoints |
| [`Lut2D`](inc/damp/toolbox/lookup.hpp#L201) | 2-D bilinear interpolating lookup table over a regular grid |
| [`Lut3D`](inc/damp/toolbox/lookup.hpp#L275) | 3-D trilinear interpolating lookup table over a regular grid |
| [`OffDelayTimer`](inc/damp/toolbox/logic.hpp#L136) | Off-delay timer: output goes true immediately when in is true and stays true until in has been false continuously for delay |
| [`OnDelayTimer`](inc/damp/toolbox/logic.hpp#L100) | On-delay timer: output goes true once in has been held true continuously for delay; drops immediately when in goes false |
| [`Periodic`](inc/damp/toolbox/timing.hpp#L170) | Periodic trigger — fires once per elapsed period |
| [`PhaseUnwrapper`](inc/damp/encoder/unwrap.hpp#L40) | Counts whole periods from successive fine phases (no hardware counter) |
| [`PulseTimer`](inc/damp/toolbox/logic.hpp#L171) | Pulse timer (non-retriggerable): a rising edge of in emits a fixed |
| [`QuadMode`](inc/damp/encoder/quadrature.hpp#L47) | Quadrature decode resolution (edges counted per A/B cycle) |
| [`QuadratureDecoder`](inc/damp/encoder/quadrature.hpp#L75) | Software A/B quadrature decoder with optional index |
| [`R_TRIG`](inc/damp/toolbox/iec61131.hpp#L119) | R_TRIG (Rising Edge Trigger) |
| [`RangeMonitor`](inc/damp/toolbox/conditioning.hpp#L403) | Analog-input range/fault monitor (NAMUR NE43 pattern) |
| [`RedundantAnalogPair`](inc/damp/toolbox/io.hpp#L91) | Dual-channel analog input with same-slope or opposite-slope cross-check |
| [`RisingEdge`](inc/damp/toolbox/logic.hpp#L34) | Rising-edge detector: true on the tick x goes false → true |
| [`RS`](inc/damp/toolbox/iec61131.hpp#L92) | RS Latch (Reset-Set Latch) |
| [`ServoAxis`](inc/damp/toolbox/actuator.hpp#L107) | One servoactuator transmission: SI joint unit ⟷ drive (motor) units |
| [`ServoBank`](inc/damp/toolbox/actuator.hpp#L233) | A bank of ServoAxis transmissions: maps a synchronized multi-axis TrajectoryState array straight to drive commands in one call |
| [`ServoCommand`](inc/damp/toolbox/actuator.hpp#L84) | A drive-native servoactuator setpoint: position, velocity, torque |
| [`SignalRange`](inc/damp/toolbox/conditioning.hpp#L352) | The measuring span and wire-fault limits of an analog input (NAMUR NE43 pattern) |
| [`SignalStatus`](inc/damp/toolbox/conditioning.hpp#L320) | Classification of an analog input against its valid/fault bands |
| [`SinCosCalibration`](inc/damp/encoder/sincos.hpp#L66) | Per-channel gain/offset trims for the sin/cos front end (hardware calibration) |
| [`SinCosSample`](inc/damp/encoder/sincos.hpp#L78) | One decoded sin/cos sample: the fine phase and the signal amplitude |
| [`SlewLimiter`](inc/damp/toolbox/conditioning.hpp#L220) | Slew-rate limiter: bound how fast the output may follow the target |
| [`SplineSurface`](inc/damp/toolbox/lookup.hpp#L394) | 2-D interpolating surface with smooth (Catmull-Rom spline) blending |
| [`SR`](inc/damp/toolbox/iec61131.hpp#L64) | SR Latch (Set-dominant Set-Reset Latch) |
| [`Stopwatch`](inc/damp/toolbox/timing.hpp#L91) | Free-running elapsed-time accumulator |
| [`TFF`](inc/damp/toolbox/iec61131.hpp#L507) | T Flip-Flop (toggle on rising edge) |
| [`Thermistor`](inc/damp/toolbox/thermistor.hpp#L184) | NTC thermistor linearization (resistance → temperature) |
| [`ThermistorResult`](inc/damp/toolbox/thermistor.hpp#L44) | Fitted NTC coefficients in Steinhart-Hart form |
| [`Timeout`](inc/damp/toolbox/timing.hpp#L126) | One-shot timeout |
| [`TOF`](inc/damp/toolbox/iec61131.hpp#L224) | TOF Timer (Timer Off Delay) |
| [`Toggle`](inc/damp/toolbox/logic.hpp#L244) | Toggle (T flip-flop): output flips on each rising edge of in |
| [`TON`](inc/damp/toolbox/iec61131.hpp#L184) | TON Timer (Timer On Delay) |
| [`TP`](inc/damp/toolbox/iec61131.hpp#L262) | TP Timer (Timer Pulse) |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`catmull_1d`](inc/damp/toolbox/lookup.hpp#L340) | 1-D non-uniform Catmull-Rom (cubic Hermite) evaluation |
| [`classify_range`](inc/damp/toolbox/conditioning.hpp#L369) | Classify x against the valid span and the wire-fault limits of range |
| [`combine_quadrature`](inc/damp/encoder/sincos.hpp#L120) | Absolute position from a quadrature edge counter and the fine phase |
| [`cross_check_analog`](inc/damp/toolbox/conditioning.hpp#L543) | Cross-check two calibrated analog readings for redundant sensors |
| [`deadband`](inc/damp/toolbox/conditioning.hpp#L49) | Dead zone over `[lower, upper]`, matching Simulink®'s Dead Zone block (+1 more overload) |
| [`decode_sincos`](inc/damp/encoder/sincos.hpp#L92) | Trim a (sin, cos) pair and decode its phase |
| [`differential_drive`](inc/damp/toolbox/io.hpp#L262) | Differential (tank) drive mixer: throttle + turn → left/right |
| [`expo`](inc/damp/toolbox/conditioning.hpp#L204) | Exponential response curve `y = (1−k)·x + k·x³` (RC "expo") |
| [`h_pattern`](inc/damp/toolbox/io.hpp#L285) | H-pattern (dual-lever / tank) drive: one lever per track, no mixing |
| [`holonomic_drive`](inc/damp/toolbox/io.hpp#L386) | Four-wheel holonomic drive mixer (mecanum / omni) |
| [`inverse_deadband`](inc/damp/toolbox/conditioning.hpp#L81) | Inverse dead zone: add an offset to overcome a physical dead zone (valve overlap, static friction, motor stiction), with independent negative/positive offsets (+1 more overload) |
| [`inverse_deadband_sine`](inc/damp/toolbox/conditioning.hpp#L164) | Sine-shaped inverse dead zone (full boost once \|x\|≥w) |
| [`inverse_deadband_soft`](inc/damp/toolbox/conditioning.hpp#L140) | Smooth inverse dead zone via soft-sign (tanh-like Coulomb boost) |
| [`inverse_lerp`](inc/damp/toolbox/scaling.hpp#L48) | Inverse of lerp: the fraction t such that `lerp(a, b, t) == x` |
| [`is_fault`](inc/damp/toolbox/conditioning.hpp#L334) | True for a wire fault (FaultLow/FaultHigh) — i.e. not a real reading at all |
| [`is_usable`](inc/damp/toolbox/conditioning.hpp#L508) | True when a channel reading may still be used (in-span or saturated) |
| [`is_valid`](inc/damp/toolbox/conditioning.hpp#L329) | True only for the in-span status |
| [`iso_c_drive`](inc/damp/toolbox/io.hpp#L304) | ISO-C (single-stick, coordinated / curvature) drive |
| [`lerp`](inc/damp/toolbox/scaling.hpp#L35) | Linear interpolation between a and b by fraction t |
| [`LIMIT`](inc/damp/toolbox/iec61131.hpp#L593) | LIMIT (IEC 61131-3 selection function): clamp in to [mn, mx] |
| [`linear_screw`](inc/damp/toolbox/actuator.hpp#L170) | Build a ServoAxis for a linear axis driven by a leadscrew/belt |
| [`lut_segment`](inc/damp/toolbox/lookup.hpp#L63) | Index of the interpolation segment containing x |
| [`MUX`](inc/damp/toolbox/iec61131.hpp#L612) | MUX (IEC 61131-3 multiplexer): select input k of N (0-based) |
| [`poly_horner`](inc/damp/toolbox/scaling.hpp#L138) | Evaluate a polynomial at x by Horner's method |
| [`rescale`](inc/damp/toolbox/scaling.hpp#L69) | Affine map of x from the input range to the output range |
| [`rotary_gearbox`](inc/damp/toolbox/actuator.hpp#L153) | Build a ServoAxis for a rotary joint behind a gearbox |
| [`scaled_deadband`](inc/damp/toolbox/conditioning.hpp#L185) | Center dead zone that rescales the surviving range back to full span |
| [`SEL`](inc/damp/toolbox/iec61131.hpp#L601) | SEL (IEC 61131-3 binary selection): g ? in1 : in0 |
| [`soft_sign`](inc/damp/toolbox/conditioning.hpp#L107) | Smooth unit direction x / √(x²+ε²) ∈ (−1,1) |
| [`steer_map`](inc/damp/toolbox/io.hpp#L344) | Curve-driven steering map — arbitrary `(throttle, turn)` geometry |
| [`steinhart_hart`](inc/damp/toolbox/thermistor.hpp#L129) | Fit the Steinhart-Hart coefficients from three calibration points |
| [`thermistor_beta`](inc/damp/toolbox/thermistor.hpp#L83) | Fit NTC coefficients from the Beta-parameter model |
| [`two_point_cal`](inc/damp/toolbox/scaling.hpp#L112) | Fit an AffineCal through two `(raw, engineering)` points |
| [`wrapped_delta`](inc/damp/encoder/quadrature.hpp#L41) | Signed difference between two unsigned counter readings, wrap-safe |

> **Host — analysis, sim, MATLAB® aliases**

## Frequency-domain analysis (host)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`BodeResult`](inc/damp/analysis/frequency.hpp#L59) | Bode plot data for a SISO system |
| [`DiskMarginResult`](inc/damp/analysis/disk_margin.hpp#L36) | Disk margin of a SISO loop: the disk size α and the gain and phase margins it guarantees |
| [`FrequencyPoint`](inc/damp/analysis/frequency.hpp#L37) | Single-point frequency response result |
| [`GangOfFourResult`](inc/damp/analysis/sensitivity.hpp#L29) | Gang of four of a loop u = C(r − y): S, T, PS and CS, each a state-space model on the closed loop's states |
| [`GangOfSevenResult`](inc/damp/analysis/sensitivity.hpp#L45) | Gang of seven of a two-degree-of-freedom loop u = K[r; y]: the gang of four of its feedback path plus the reference responses |
| [`HinfNormResult`](inc/damp/analysis/norms.hpp#L119) | H∞ norm of a stable system and the frequency of its peak |
| [`ImpedanceResult`](inc/damp/analysis/frequency.hpp#L903) | Result of impedance frequency response evaluation |
| [`LoopDiskMarginResult`](inc/damp/analysis/disk_margin.hpp#L55) | Loop-at-a-time disk margins of a square MIMO loop: one SISO disk margin per channel |
| [`LoopResponseResult`](inc/damp/analysis/frequency.hpp#L417) | Open-loop and closed-loop frequency response package |
| [`LoopSummary`](inc/damp/analysis/frequency.hpp#L450) | Compact loop summary metrics for quick stability/robustness checks |
| [`LsimInfo`](inc/damp/analysis/time_response.hpp#L513) | Transient characteristics of an arbitrary response signal |
| [`LsimResult`](inc/damp/analysis/time_response.hpp#L66) | Result of a single-trajectory simulation: time, output, and state history |
| [`MiddlebrookResult`](inc/damp/analysis/frequency.hpp#L927) | Result of Middlebrook minor loop gain analysis |
| [`NicholsPoint`](inc/damp/analysis/frequency.hpp#L685) | Single-point Nichols chart sample (open-loop phase vs magnitude) |
| [`NicholsResult`](inc/damp/analysis/frequency.hpp#L697) | Nichols chart locus across a frequency sweep |
| [`NyquistPoint`](inc/damp/analysis/frequency.hpp#L371) | Single-point Nyquist response data |
| [`NyquistResult`](inc/damp/analysis/frequency.hpp#L383) | Nyquist response data across a frequency sweep |
| [`PoleInfo`](inc/damp/analysis/poles.hpp#L56) | Natural frequency and damping ratio for each pole |
| [`PoleZeroMap`](inc/damp/analysis/poles.hpp#L102) | Poles and zeros of a system, for pole-zero plotting |
| [`RootLocusResult`](inc/damp/analysis/poles.hpp#L234) | Root-locus data: closed-loop poles along a gain grid |
| [`SigmaPoint`](inc/damp/analysis/frequency.hpp#L770) | Singular values of G(jω) at one frequency |
| [`SigmaResult`](inc/damp/analysis/frequency.hpp#L784) | Singular-value frequency response over a grid |
| [`StepInfo`](inc/damp/analysis/time_response.hpp#L381) | Step-response characteristics of a single output signal |
| [`TimeResponse`](inc/damp/analysis/time_response.hpp#L51) | Multi-channel time-domain response sampled on a time grid |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`arange`](inc/damp/analysis/linspace.hpp#L115) | Half-open arithmetic range [start, stop) with step step (+2 more overloads) |
| [`bode`](inc/damp/analysis/frequency.hpp#L267) | Compute Bode plot data for a SISO state-space system (+2 more overloads) |
| [`canonical_phase_margin`](inc/damp/analysis/frequency.hpp#L150) | Normalize phase margin to (-180, 180] |
| [`damp`](inc/damp/analysis/poles.hpp#L74) | Compute natural frequency and damping for each pole |
| [`dcgain`](inc/damp/analysis/norms.hpp#L42) | Compute DC gain of a continuous-time system |
| [`disk_margin`](inc/damp/analysis/disk_margin.hpp#L158) | Disk margin of a SISO loop L under unity negative feedback |
| [`frf_s_or_z`](inc/damp/analysis/frequency.hpp#L1005) | FRF evaluation point: jω (continuous) or e^{jωTs} (discrete) |
| [`gain_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L211) | Find the gain margin using the unwrapped phase trajectory |
| [`gang_of_four`](inc/damp/analysis/sensitivity.hpp#L207) | Gang of four of a plant P and a controller C on the error, u = C(r − y) |
| [`gang_of_seven`](inc/damp/analysis/sensitivity.hpp#L150) | Gang of seven of a plant and a two-degree-of-freedom controller u = K[r; y] |
| [`geomspace`](inc/damp/analysis/linspace.hpp#L68) | Geometric sequence from start to end (endpoint values, not exponents) (+1 more overload) |
| [`impedance`](inc/damp/analysis/frequency.hpp#L1073) | Compute impedance frequency response from a SISO impedance transfer function |
| [`impedance_from_admittance`](inc/damp/analysis/frequency.hpp#L1031) | Compute impedance frequency response from a SISO admittance system |
| [`impulse`](inc/damp/analysis/time_response.hpp#L198) | Impulse response of a (MIMO) state-space system (+1 more overload) |
| [`initial`](inc/damp/analysis/time_response.hpp#L246) | Initial-condition (free) response of a (MIMO) state-space system |
| [`linspace`](inc/damp/analysis/linspace.hpp#L36) | num evenly spaced values from start to end inclusive (+1 more overload) |
| [`loop_at_a_time_disk_margin`](inc/damp/analysis/disk_margin.hpp#L194) | Loop-at-a-time disk margins of a square MIMO loop L under unity negative feedback |
| [`loop_metrics`](inc/damp/analysis/frequency.hpp#L482) | One-call loop analysis: compute L/S/T response and return compact metrics (+1 more overload) |
| [`loop_response`](inc/damp/analysis/frequency.hpp#L584) | Compute open-loop L, sensitivity S, complementary sensitivity T, and Nyquist data (+1 more overload) |
| [`lsim`](inc/damp/analysis/time_response.hpp#L284) | Forced time response of a (MIMO) state-space system to an input signal (+1 more overload) |
| [`lsiminfo`](inc/damp/analysis/time_response.hpp#L533) | Compute transient characteristics from an output/time signal |
| [`middlebrook`](inc/damp/analysis/frequency.hpp#L1113) | Middlebrook stability analysis of a cascaded source and load from their impedances |
| [`nichols`](inc/damp/analysis/frequency.hpp#L714) | Build Nichols points from existing Bode data (+2 more overloads) |
| [`norm_h2`](inc/damp/analysis/norms.hpp#L78) | H2 norm of a state-space system |
| [`norm_hinf`](inc/damp/analysis/norms.hpp#L337) | H∞ norm of a stable state-space system, sup σ̄(G) over the imaginary axis or the unit circle, and its frequency |
| [`nyquist`](inc/damp/analysis/frequency.hpp#L499) | Compute Nyquist data for a SISO state-space system (+1 more overload) |
| [`phase_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L169) | Find the phase margin using the unwrapped phase trajectory |
| [`poles`](inc/damp/analysis/poles.hpp#L40) | Compute open-loop poles (eigenvalues of A matrix) |
| [`poly_roots`](inc/damp/analysis/poles.hpp#L150) | Roots of a polynomial given highest power first (MATLAB® `roots`) |
| [`pzmap`](inc/damp/analysis/poles.hpp#L173) | Pole-zero map of a SISO transfer function (MATLAB® `pzmap(tf)`) (+2 more overloads) |
| [`rlocus`](inc/damp/analysis/poles.hpp#L293) | Root locus of a SISO state-space plant over an explicit gain grid (+1 more overload) |
| [`sigma`](inc/damp/analysis/frequency.hpp#L827) | Singular-value frequency response of a (possibly MIMO) state-space system (+1 more overload) |
| [`step`](inc/damp/analysis/time_response.hpp#L164) | Step response of a (MIMO) state-space system (+1 more overload) |
| [`stepinfo`](inc/damp/analysis/time_response.hpp#L402) | Compute step-response characteristics from an output/time signal (+1 more overload) |
| [`summarize_loop_response`](inc/damp/analysis/frequency.hpp#L462) | Summarize loop_response() results into one compact metrics struct |
| [`unwrap_phase_deg`](inc/damp/analysis/frequency.hpp#L127) | Unwrap phase data in degrees to avoid +/-180 discontinuities |

## Simulation / SIL harness (host)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AdaptiveOptions`](inc/damp/simulation/solver.hpp#L349) | Options for adaptive-step ODE integration |
| [`AdaptiveStepIntegrator`](inc/damp/simulation/integrator.hpp#L115) | Integrator that reports a genuine embedded local-error estimate |
| [`AdaptiveStepSolver`](inc/damp/simulation/solver.hpp#L404) | Adaptive-step ODE solver (+1 more overload) |
| [`BackwardEuler`](inc/damp/simulation/integrator.hpp#L343) | Backward Euler integrator |
| [`BDF2`](inc/damp/simulation/integrator.hpp#L403) | Backward Differentiation Formula 2 (BDF2) integrator |
| [`CachedZoh`](inc/damp/simulation/cached_zoh.hpp#L63) | Cached discrete ZOH maps Ad, Bd for step size h |
| [`Discrete`](inc/damp/simulation/integrator.hpp#L133) | Discrete-time integrator (no integration, just one step) |
| [`DP45`](inc/damp/simulation/integrator.hpp#L780) | Dormand-Prince 5(4) adaptive integrator — the ode45 pair |
| [`ErrorTolerance`](inc/damp/simulation/solver.hpp#L302) | Scalar-or-vector local-error tolerance for adaptive ODE control |
| [`EventSimConfig`](inc/damp/simulation/hybrid.hpp#L67) | Configuration for event-driven hybrid / multi-rate style runs |
| [`Exact`](inc/damp/simulation/integrator.hpp#L167) | Exact integrator for LTI systems |
| [`FixedStepSolver`](inc/damp/simulation/solver.hpp#L145) | Fixed-step ODE solver (+1 more overload) |
| [`ForwardEuler`](inc/damp/simulation/integrator.hpp#L223) | Forward Euler integrator |
| [`Heun`](inc/damp/simulation/integrator.hpp#L863) | Heun's method (Improved Euler, RK2) integrator |
| [`HybridEventAction`](inc/damp/simulation/hybrid.hpp#L104) | Action applied when an event fires (after Exact advance to the event time) |
| [`HybridSimulationResult`](inc/damp/simulation/hybrid.hpp#L79) | Result of a piecewise-LTI hybrid simulation |
| [`IntegrationResult`](inc/damp/simulation/integrator.hpp#L64) | Result of an integration step |
| [`LogPolicy`](inc/damp/simulation/hybrid.hpp#L58) | How often to record samples along a hybrid trajectory |
| [`MultiRateConfig`](inc/damp/simulation/multirate.hpp#L81) | Configuration for multi-rate closed-loop runs |
| [`PiecewiseLTI`](inc/damp/simulation/hybrid.hpp#L93) | Plant: NModes continuous LTI pieces (same state/input/output dimensions) |
| [`RK23`](inc/damp/simulation/integrator.hpp#L918) | Bogacki-Shampine 2(3) adaptive integrator |
| [`RK3`](inc/damp/simulation/integrator.hpp#L985) | Classical 3rd-order Runge-Kutta (RK3) integrator |
| [`RK4`](inc/damp/simulation/integrator.hpp#L710) | Classical 4th-order Runge-Kutta (RK4) integrator |
| [`SimulationResult`](inc/damp/simulation/simulate.hpp#L128) | Result of a closed-loop simulation |
| [`SisoReferenceAdapter`](inc/damp/simulation/simulate.hpp#L96) | Adapt a SISO `control(r, y)` controller for callables that want `u = f(y)` or `u = f(t, y)` (the older simulate contracts) |
| [`SolveResult`](inc/damp/simulation/solver.hpp#L86) | Result of an ODE solve operation |
| [`SymplecticEuler`](inc/damp/simulation/integrator.hpp#L278) | Semi-implicit (symplectic) Euler for mechanical systems |
| [`Trapezoidal`](inc/damp/simulation/integrator.hpp#L652) | Trapezoidal (Tustin) integrator |
| [`TRBDF2`](inc/damp/simulation/integrator.hpp#L527) | TR-BDF2 composite integrator with an embedded error estimate (the ode23tb scheme) |
| [`TwoRateSimulationResult`](inc/damp/simulation/multirate.hpp#L108) | Result of a two-rate cascade simulation |
| [`XyPlotOpts`](inc/damp/simulation/plot_plotly.hpp#L537) | Optional axis framing for plot_xy |
| [`XySeries`](inc/damp/simulation/plot_plotly.hpp#L519) | One series for a planar (x, y) scatter / path plot |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`adaptive_solve`](inc/damp/simulation/solver.hpp#L814) | Adaptive-step solve in one call |
| [`apply_wet_plotly_shell`](inc/damp/simulation/plot_plotly.hpp#L66) | Apply damp's full-window HTML shell to plotlypp-generated markup |
| [`bodemag`](inc/damp/simulation/plot_plotly.hpp#L798) | Plot a magnitude-only Bode diagram (log frequency, dB magnitude) |
| [`bodeplot`](inc/damp/simulation/plot_plotly.hpp#L784) | Plot magnitude and phase Bode subplots |
| [`exact_lti_step`](inc/damp/simulation/hybrid.hpp#L303) | Exact open-loop step of a single continuous LTI system |
| [`extract_channel`](inc/damp/simulation/plot_plotly.hpp#L225) | Extract the i-th element from each vector entry into a std::vector<double> |
| [`fires_at`](inc/damp/simulation/multirate.hpp#L73) | True when base tick tick is a firing instant for divisor |
| [`fixed_solve`](inc/damp/simulation/solver.hpp#L286) | Fixed-step solve in one call |
| [`impulseplot`](inc/damp/simulation/plot_plotly.hpp#L735) | Plot an impulse response, one trace per input/output pair |
| [`locate_zero_crossing_exact`](inc/damp/simulation/hybrid.hpp#L333) | Locate a state zero-crossing of g(x) on an Exact LTI segment |
| [`lsimplot`](inc/damp/simulation/plot_plotly.hpp#L751) | Plot a forced (lsim) simulation, one trace per output |
| [`nicholsplot`](inc/damp/simulation/plot_plotly.hpp#L944) | Plot a Nichols chart (open-loop phase vs magnitude) with M-circle grid |
| [`nyquistplot`](inc/damp/simulation/plot_plotly.hpp#L823) | Plot a Nyquist locus with the -1 critical point marked |
| [`panel_color`](inc/damp/simulation/plot_plotly.hpp#L316) | Plotly default colorway entry for index i within a single subplot |
| [`panel_legend`](inc/damp/simulation/plot_plotly.hpp#L299) | Per-panel legend box just right of a stacked subplot band |
| [`panel_legend_id`](inc/damp/simulation/plot_plotly.hpp#L285) | Legend id for subplot row (1-based): "legend", "legend2", … |
| [`panel_x_domain`](inc/damp/simulation/plot_plotly.hpp#L275) | Horizontal domain for stacked multi-panel x-axes |
| [`plot_bode`](inc/damp/simulation/plot_plotly.hpp#L424) | Plot Bode magnitude and phase as subplots |
| [`plot_line`](inc/damp/simulation/plot_plotly.hpp#L493) | Simple line plot of time vs value |
| [`plot_simulation`](inc/damp/simulation/plot_plotly.hpp#L345) | Plot simulation results with subplots for states, outputs, and inputs |
| [`plot_step`](inc/damp/simulation/plot_plotly.hpp#L623) | Plot step response data |
| [`plot_xy`](inc/damp/simulation/plot_plotly.hpp#L558) | Multi-series planar scatter (true 2-D path / locus) (+1 more overload) |
| [`plots_root_for`](inc/damp/simulation/plot_plotly.hpp#L120) | Walk up from dir until a directory named "plots" is found |
| [`pzplot`](inc/damp/simulation/plot_plotly.hpp#L850) | Plot a pole-zero map on the complex plane (poles as ×, zeros as ○) |
| [`quantize`](inc/damp/simulation/quantize.hpp#L66) | ADC reading of x: x rounded to the nearest multiple of lsb, then clamped to [lo, hi] |
| [`rlocusplot`](inc/damp/simulation/plot_plotly.hpp#L1042) | Plot a root locus (closed-loop poles vs gain) on the complex plane |
| [`sigmaplot`](inc/damp/simulation/plot_plotly.hpp#L999) | Plot singular-value frequency response (log frequency, dB) |
| [`simulate`](inc/damp/simulation/simulate.hpp#L159) | Simulate a nonlinear plant with a controller in closed loop |
| [`simulate_cached_zoh`](inc/damp/simulation/cached_zoh.hpp#L157) | Simulate a continuous LTI plant with fixed ZOH step and held input (+1 more overload) |
| [`simulate_cached_zoh_multirate`](inc/damp/simulation/cached_zoh.hpp#L236) | Multi-rate ZOH: held control samples, dense plant samples for plotting |
| [`simulate_discrete`](inc/damp/simulation/simulate.hpp#L525) | Simulate a discrete-time system with a controller |
| [`simulate_discrete_nonlinear`](inc/damp/simulation/simulate.hpp#L478) | Simulate a discrete-time nonlinear plant with a controller |
| [`simulate_lti`](inc/damp/simulation/simulate.hpp#L359) | Simulate a continuous LTI system with a controller |
| [`simulate_lti_siso`](inc/damp/simulation/simulate.hpp#L421) | SIL: continuous SISO LTI plant under a sampled SISO controller with constant reference |
| [`simulate_multirate`](inc/damp/simulation/multirate.hpp#L142) | Multi-rate closed-loop simulation with MultiRateConfig (+1 more overload) |
| [`simulate_piecewise_lti`](inc/damp/simulation/hybrid.hpp#L185) | Simulate a piecewise-LTI plant with Exact integration between events |
| [`simulate_sampled`](inc/damp/simulation/simulate.hpp#L232) | Simulate a continuous plant under a discrete (sampled) controller — multi-rate |
| [`simulate_state_feedback`](inc/damp/simulation/simulate.hpp#L303) | Simulate a nonlinear plant with state-feedback controller |
| [`simulate_two_rate`](inc/damp/simulation/multirate.hpp#L266) | Two-rate cascade with MultiRateConfig (+1 more overload) |
| [`simulate_two_rate_lti`](inc/damp/simulation/multirate.hpp#L415) | Two-rate cascade on a continuous LTI plant (A, B, C, D) |
| [`siso_ref`](inc/damp/simulation/simulate.hpp#L115) | Build a SisoReferenceAdapter for a SISO `control(r,y)` controller |
| [`stepplot`](inc/damp/simulation/plot_plotly.hpp#L719) | Plot a step response, one trace per input/output pair |
| [`to_double_vector`](inc/damp/simulation/plot_plotly.hpp#L238) | Convert std::vector\<T\> to std::vector<double> |
| [`to_std_vector`](inc/damp/simulation/plot_plotly.hpp#L213) | Convert a ColVec<N,T> to std::vector<double> for plotlypp |
| [`write_html`](inc/damp/simulation/plot_plotly.hpp#L147) | Write a plotlypp figure to HTML with damp's shell (prefer over writeHtml) |
| [`write_hybrid_trace_npy`](inc/damp/simulation/npy_export.hpp#L101) | Pack HybridSimulationResult into (N, 2+NX+NU+NY) float64: t, mode, x..., u..., y |
| [`write_npy_f64`](inc/damp/simulation/npy_export.hpp#L42) | Write a C-order float64 array as a .npy file |

## MATLAB®-style aliases (host)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AllMarginResult`](inc/damp/matlab.hpp#L1293) | All classical margins including delay margin (superset of margin) |
| [`MarginResult`](inc/damp/matlab.hpp#L1245) | Gain/phase margins and their crossover frequencies |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`acker`](inc/damp/matlab.hpp#L849) | Single-input pole placement (MATLAB®'s acker) |
| [`allmargin`](inc/damp/matlab.hpp#L1314) | Gain, phase, and delay margins of a SISO loop over a frequency grid |
| [`balreal`](inc/damp/matlab.hpp#L343) | MATLAB® short alias for design::balanced_realization |
| [`balred`](inc/damp/matlab.hpp#L362) | MATLAB® short alias for design::balanced_reduction to NR states |
| [`bandwidth`](inc/damp/matlab.hpp#L1402) | -3 dB bandwidth of a SISO system over a frequency grid |
| [`c2d`](inc/damp/matlab.hpp#L425) | MATLAB® interface function c2d to discretize a continuous-time state-space system (+1 more overload) |
| [`canon`](inc/damp/matlab.hpp#L500) | MATLAB® canon(sys) / canon(sys, 'modal'): real modal canonical form and its transform (+1 more overload) |
| [`ctrb`](inc/damp/matlab.hpp#L386) | MATLAB® short alias for controllability_matrix (+1 more overload) |
| [`ctrbf`](inc/damp/matlab.hpp#L541) | MATLAB® ctrbf(A, B, C): controllability staircase form, its orthogonal transform and rank |
| [`d2d`](inc/damp/matlab.hpp#L458) | MATLAB® d2d: resample a ZOH-discrete model to a new sample time |
| [`diskmargin`](inc/damp/matlab.hpp#L1331) | Disk margin of a SISO loop (+1 more overload) |
| [`dlqe`](inc/damp/matlab.hpp#L959) | Discrete linear-quadratic estimator from matrices (MATLAB®'s dlqe) |
| [`dlqr`](inc/damp/matlab.hpp#L907) | Discrete-time Linear-Quadratic Regulator design (MATLAB®'s dlqr) |
| [`dlyap`](inc/damp/matlab.hpp#L1422) | MATLAB® alias for the discrete Lyapunov solve AXAᵀ−X+Q=0 |
| [`dlyapchol`](inc/damp/matlab.hpp#L1442) | MATLAB® alias for the Cholesky factor of the discrete Lyapunov solution (X = RᵀR) |
| [`eig`](inc/damp/matlab.hpp#L669) | MATLAB® short alias for the eigenvalues of a square matrix |
| [`estim`](inc/damp/matlab.hpp#L731) | Form the state estimator from a plant and an estimator gain |
| [`expm`](inc/damp/matlab.hpp#L681) | MATLAB® short alias for the matrix exponential (mat::expm) |
| [`eye`](inc/damp/matlab.hpp#L609) | Create an identity matrix of size n x n |
| [`gram`](inc/damp/matlab.hpp#L1498) | MATLAB® alias for the controllability/observability Gramian of a system |
| [`hankelsv`](inc/damp/matlab.hpp#L352) | MATLAB® short alias for design::hankel_singular_values |
| [`hinfnorm`](inc/damp/matlab.hpp#L1521) | MATLAB® alias for the H∞ system norm [ninf, fpeak] = hinfnorm(sys, tol) |
| [`isstable`](inc/damp/matlab.hpp#L1355) | Stability predicate for a state-space model (+1 more overload) |
| [`kalmd`](inc/damp/matlab.hpp#L982) | Discrete Kalman estimator for a continuous plant and continuous noise (MATLAB®'s kalmd) |
| [`linmod`](inc/damp/matlab.hpp#L588) | MATLAB®-style nonlinear linearization about an operating point |
| [`logm`](inc/damp/matlab.hpp#L707) | MATLAB® short alias for the principal matrix logarithm (mat::log) |
| [`logspace`](inc/damp/matlab.hpp#L1217) | n points from 10ᵃ to 10ᵇ, logarithmically spaced (decade exponents, as MATLAB®'s) |
| [`lqg`](inc/damp/matlab.hpp#L1097) | Linear-Quadratic-Gaussian regulator design |
| [`lqgreg`](inc/damp/matlab.hpp#L1113) | Combine separate Kalman filter and LQR designs into an LQG controller |
| [`lqgtrack`](inc/damp/matlab.hpp#L1134) | Linear-Quadratic-Gaussian design with integral action for tracking (+2 more overloads) |
| [`lqi`](inc/damp/matlab.hpp#L1058) | Linear-Quadratic Integral design for tracking (+2 more overloads) |
| [`lqr`](inc/damp/matlab.hpp#L892) | Continuous-time LQR design (MATLAB®'s lqr) |
| [`lqrd`](inc/damp/matlab.hpp#L922) | Sampled-data LQR from a continuous plant: ZOH dynamics + Van Loan cost (MATLAB®'s lqrd) (+1 more overload) |
| [`lqry`](inc/damp/matlab.hpp#L1008) | Output-weighted continuous LQR (state cost Q = Cᵀ Q_y C) (+1 more overload) |
| [`lyap`](inc/damp/matlab.hpp#L1412) | MATLAB® alias for the continuous Lyapunov solve AX+XAᵀ+Q=0 |
| [`lyapchol`](inc/damp/matlab.hpp#L1432) | MATLAB® alias for the Cholesky factor of the continuous Lyapunov solution (X = RᵀR) |
| [`make1DOF`](inc/damp/matlab.hpp#L304) | Force 1-DOF setpoint weights on a PID design result (b=c=1) |
| [`make2DOF`](inc/damp/matlab.hpp#L320) | Apply 2-DOF setpoint weights on a PID design result |
| [`margin`](inc/damp/matlab.hpp#L1266) | Gain and phase margins of a SISO loop over a frequency grid |
| [`minreal`](inc/damp/matlab.hpp#L334) | MATLAB® short alias for design::minreal, with its default tolerance √ε |
| [`modred`](inc/damp/matlab.hpp#L375) | MATLAB® short alias for design::eliminate_states (0-based indices) |
| [`norm`](inc/damp/matlab.hpp#L1511) | MATLAB® alias for the H2 system norm norm(sys,2) |
| [`null`](inc/damp/matlab.hpp#L654) | MATLAB® short alias for an orthonormal null-space basis |
| [`obsv`](inc/damp/matlab.hpp#L402) | MATLAB® short alias for observability_matrix (+1 more overload) |
| [`obsvf`](inc/damp/matlab.hpp#L551) | MATLAB® obsvf(A, B, C): observability staircase form, its orthogonal transform and rank |
| [`ones`](inc/damp/matlab.hpp#L618) | N×N matrix of ones |
| [`pade`](inc/damp/matlab.hpp#L1384) | Padé approximation of order N (1 or 2) of the pure delay e^{−sT} |
| [`pid`](inc/damp/matlab.hpp#L249) | MATLAB®-style parallel-form continuous PID constructor |
| [`pidstd`](inc/damp/matlab.hpp#L272) | Standard-form continuous PID constructor (1-DOF) |
| [`pidstd2`](inc/damp/matlab.hpp#L291) | Standard-form continuous 2-DOF PID constructor |
| [`pidtune`](inc/damp/matlab.hpp#L1182) | Plant-aware PID tune at crossover ωc (MATLAB® pidtune spelling): 60° PID |
| [`pinv`](inc/damp/matlab.hpp#L641) | MATLAB® short alias for the Moore–Penrose pseudoinverse |
| [`place`](inc/damp/matlab.hpp#L873) | Robust multi-input pole placement (MATLAB®'s place) |
| [`pole`](inc/damp/matlab.hpp#L1229) | MATLAB® short alias for the open-loop poles of a system |
| [`prescale`](inc/damp/matlab.hpp#L481) | MATLAB® prescale: diagonal state scaling that balances A, B and C |
| [`quadprog`](inc/damp/matlab.hpp#L1453) | MATLAB® alias for the dense inequality-constrained QP solve (+1 more overload) |
| [`reg`](inc/damp/matlab.hpp#L762) | Form the observer-based compensator y → u from a plant, K and L |
| [`sqrtm`](inc/damp/matlab.hpp#L694) | MATLAB® short alias for the principal matrix square root (mat::sqrt) |
| [`ss`](inc/damp/matlab.hpp#L121) | MATLAB®-style state-space model constructor |
| [`ss2ss`](inc/damp/matlab.hpp#L531) | MATLAB® ss2ss(sys, T): the model in the state z = T·x |
| [`ss2zpk`](inc/damp/matlab.hpp#L227) | MATLAB® short alias for SISO to_zpk (state space → zero-pole-gain) |
| [`ssbal`](inc/damp/matlab.hpp#L490) | MATLAB® ssbal: diagonal state balancing of [A B; C 0] |
| [`svd`](inc/damp/matlab.hpp#L632) | MATLAB® short alias for the singular value decomposition |
| [`tf`](inc/damp/matlab.hpp#L85) | MATLAB®-style transfer function constructor (+1 more overload) |
| [`tf2sos`](inc/damp/matlab.hpp#L182) | MATLAB®-style transfer function → second-order sections, from MATLAB®'s z⁻¹ vectors (+1 more overload) |
| [`tf2zpk`](inc/damp/matlab.hpp#L157) | MATLAB® short alias for to_zpk (transfer function → zero-pole-gain) |
| [`thiran`](inc/damp/matlab.hpp#L572) | MATLAB® thiran: fractional-delay all-pass of order N for a delay tau at Ts |
| [`tzero`](inc/damp/matlab.hpp#L561) | MATLAB® tzero(sys): the finite invariant zeros of a model with any number of inputs and outputs |
| [`upsample`](inc/damp/matlab.hpp#L471) | MATLAB® upsample: sys(z) → sys(zᴸ) at Ts / L |
| [`zp2sos`](inc/damp/matlab.hpp#L218) | MATLAB®-style zero-pole-gain → second-order sections, H(z) = k·Π(z − zᵢ)/Π(z − pⱼ) |
| [`zpk`](inc/damp/matlab.hpp#L144) | MATLAB®-style zero-pole-gain model constructor |
| [`zpk2tf`](inc/damp/matlab.hpp#L166) | MATLAB® short alias for ZPK::to_transfer_function |

## Math backends

Internal, compile-time-selected implementations of the `damp::` scalar-math surface (`sin`/`cos`/`sqrt`/`exp`/…), chosen via `damp/config.hpp`. Every backend exposes the same functions, so they're listed once here as files rather than repeated in the tables above. The public dispatcher is [`damp/math/math.hpp`](inc/damp/math/math.hpp).

| File | Role |
| ---- | ---- |
| [`damp/math/math_backend.hpp`](inc/damp/math/math_backend.hpp) | Pluggable runtime math backend selection (freestanding-safe) |
| [`damp/math/std_fallback.hpp`](inc/damp/math/std_fallback.hpp) | Composable std:: (\<cmath\>) base for the hosted math backends |
| [`damp/math/damp_backend.hpp`](inc/damp/math/damp_backend.hpp) | Fast float and double math backend (trig.hpp) over the std:: fallback |
| [`damp/math/series_backend.hpp`](inc/damp/math/series_backend.hpp) | Freestanding runtime math backend: routes scalar math to the constexpr series in constexpr_math.hpp, pulling no hosted headers (no \<cmath\>) |
| [`damp/math/constexpr_math.hpp`](inc/damp/math/constexpr_math.hpp) | Compile-time scalar math for the damp:: dispatch layer |
| [`damp/math/trig.hpp`](inc/damp/math/trig.hpp) | Fast float and double sine, cosine and inverse trig, with Cody–Waite reduction |

## Examples

Runnable programs in `examples/` (37 total). Build with `make` (or `tup --quiet examples`); outputs go to `examples/build/`.

| Example | Description |
| ------- | ----------- |
| [`adrc_sil.cpp`](examples/control/adrc/adrc_sil.cpp) | ADRC vs well-tuned PI-D on a rigid inertia (host SIL + plot) |
| [`adrc_sketch.cpp`](examples/control/adrc/adrc_sketch.cpp) | Position ADRC — flashable sketch (Design Is Deploy) |
| [`cart_pole_sil.cpp`](examples/control/cart_pole/cart_pole_sil.cpp) | Cart-pole LQR — host nonlinear SIL (calls cart_pole_controller.hpp) |
| [`cart_pole_sketch.cpp`](examples/control/cart_pole/cart_pole_sketch.cpp) | Cart-pole LQR — flashable sketch (Design Is Deploy) |
| [`encoder_velocity_sil.cpp`](examples/estimation/encoder_velocity/encoder_velocity_sil.cpp) | Encoder / tach velocity estimation — host comparison + plots |
| [`encoder_velocity_sketch.cpp`](examples/estimation/encoder_velocity/encoder_velocity_sketch.cpp) | Encoder PLL velocity — thin deploy-shaped smoke |
| [`eskf_sil.cpp`](examples/estimation/eskf/eskf_sil.cpp) | MARG attitude ESKF — finite-tick host smoke (calls eskf_estimator.hpp) |
| [`eskf_sketch.cpp`](examples/estimation/eskf/eskf_sketch.cpp) | MARG attitude ESKF — flashable sketch (Design Is Deploy) |
| [`example_math_backend.cpp`](examples/example_math_backend.cpp) | Pluggable math backend example |
| [`fo_plant_inertia_sil.cpp`](examples/estimation/fo_plant_inertia/fo_plant_inertia_sil.cpp) | Online J,b from FirstOrderPlantEstimator — host identification demo |
| [`fo_plant_inertia_sketch.cpp`](examples/estimation/fo_plant_inertia/fo_plant_inertia_sketch.cpp) | FO plant inertia ID — thin deploy-shaped smoke |
| [`foc_sil.cpp`](examples/motor/foc/foc_sil.cpp) | Host teaching plot: FOC current-loop PI vs I-P (calls foc_controller.hpp) |
| [`foc_sketch.cpp`](examples/motor/foc/foc_sketch.cpp) | FOC current-loop — thin smoke sketch (host teaching demo) |
| [`heater_lqg_sil.cpp`](examples/control/heater_lqg/heater_lqg_sil.cpp) | Heater with dead time — host SIL for the offset-free LQG, the dead time in one sample: the flashable tick on a hotend with the true transport delay, a colder room than assumed, and a fan at 240 s |
| [`heater_lqg_sketch.cpp`](examples/control/heater_lqg/heater_lqg_sketch.cpp) | Heater with dead time — offset-free LQG, the dead time in one sample, flashable sketch |
| [`heater_ofmpc_sil.cpp`](examples/control/heater_ofmpc/heater_ofmpc_sil.cpp) | Heater with dead time — host SIL: the delay-aware MPC vs the same MPC without the delay |
| [`heater_ofmpc_sketch.cpp`](examples/control/heater_ofmpc/heater_ofmpc_sketch.cpp) | Heater with dead time — offset-free MPC, flashable sketch (Design Is Deploy) |
| [`imu_pose_sil.cpp`](examples/estimation/imu_pose/imu_pose_sil.cpp) | Animated body-frame triad from gyro integration (IMU attitude only) |
| [`imu_pose_sketch.cpp`](examples/estimation/imu_pose/imu_pose_sketch.cpp) | Open-loop gyro attitude — thin smoke (teaching) |
| [`ins_eskf_sil.cpp`](examples/estimation/ins_eskf/ins_eskf_sil.cpp) | Animated INS: free-run bias vs InsNavigator (position + dual-antenna heading) |
| [`ins_eskf_sketch.cpp`](examples/estimation/ins_eskf/ins_eskf_sketch.cpp) | InsNavigator thin smoke — same design as ins_eskf_sil (float deploy) |
| [`ins_mechanization_sil.cpp`](examples/estimation/ins_mechanization/ins_mechanization_sil.cpp) | Animated strapdown INS path: IMU samples → p, v, q (mechanization) |
| [`ins_mechanization_sketch.cpp`](examples/estimation/ins_mechanization/ins_mechanization_sketch.cpp) | Strapdown mechanization — thin rest-frame smoke |
| [`ins_navigator_sil.cpp`](examples/estimation/ins_navigator/ins_navigator_sil.cpp) | 15-state INS navigator — finite-tick host smoke (calls ins_navigator_estimator.hpp) |
| [`ins_navigator_sketch.cpp`](examples/estimation/ins_navigator/ins_navigator_sketch.cpp) | 15-state INS navigator — flashable sketch (Design Is Deploy) |
| [`lpf_sil.cpp`](examples/control/lpf/lpf_sil.cpp) | Low-pass filter — host step-response demo (calls lpf_filter.hpp) |
| [`lpf_sketch.cpp`](examples/control/lpf/lpf_sketch.cpp) | Low-pass filter — flashable sketch |
| [`lqi_sil.cpp`](examples/control/lqi/lqi_sil.cpp) | Position LQI — host plant, same tick as the sketch; clamp-aware vs naive |
| [`lqi_sketch.cpp`](examples/control/lqi/lqi_sketch.cpp) | Position LQI with a clamp — flashable sketch |
| [`pendulum_sil.cpp`](examples/control/pendulum/pendulum_sil.cpp) | Upright pendulum LQR — host nonlinear SIL (calls pendulum_controller.hpp) |
| [`pendulum_sketch.cpp`](examples/control/pendulum/pendulum_sketch.cpp) | Upright pendulum LQR — flashable sketch (Design Is Deploy) |
| [`pid_sil.cpp`](examples/control/pid/pid_sil.cpp) | PI — thin host smoke (calls pid_controller.hpp) |
| [`pid_sketch.cpp`](examples/control/pid/pid_sketch.cpp) | PI — flashable sketch (Design Is Deploy) |
| [`reaction_wheel_sil.cpp`](examples/control/reaction_wheel/reaction_wheel_sil.cpp) | Reaction-wheel gallery — print gains + one tick each (host) |
| [`reaction_wheel_sketch.cpp`](examples/control/reaction_wheel/reaction_wheel_sketch.cpp) | Reaction-wheel cascade PI — flashable sketch (Design Is Deploy) |
| [`workflow_end_to_end_sil.cpp`](examples/control/workflow_end_to_end/workflow_end_to_end_sil.cpp) | Workflow end-to-end host SIL (calls workflow_end_to_end_controller.hpp) |
| [`workflow_end_to_end_sketch.cpp`](examples/control/workflow_end_to_end/workflow_end_to_end_sketch.cpp) | Workflow e2e — flashable float runtime smoke (after host design) |
