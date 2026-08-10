# API Reference

Auto-generated from `@brief` doc comments in `inc/damp/`. Regenerate with `python tools/gen_reference.py`. Flat A→Z view: [REFERENCE_INDEX.md](REFERENCE_INDEX.md).

Compile-time (or init-time) control design in the firmware tree (variant gains as `constexpr`); same runtime objects in SIL. Scope: [docs/known_limitations.md](docs/known_limitations.md).

`damp::design::` is the design-time shelf (not for PWM-rate): heavy solvers and Result structs. Domain packs (motor / power / full motion kits) are not assumed present in the public core tree.


- [API Reference](#api-reference)
  - [Core, configuration \& backend vocabulary](#core-configuration--backend-vocabulary)
  - [Scalar math \& complex](#scalar-math--complex)
  - [Linear algebra](#linear-algebra)
  - [LTI systems (SS / TF / ZPK / discretize)](#lti-systems-ss--tf--zpk--discretize)
  - [Design-time synthesis (not PWM-rate)](#design-time-synthesis-not-pwm-rate)
  - [Runtime controllers](#runtime-controllers)
  - [Observers \& estimators](#observers--estimators)
  - [Filters \& signal conditioning](#filters--signal-conditioning)
  - [Trajectory value types](#trajectory-value-types)
  - [Kinematics / pose](#kinematics--pose)
  - [Motor control pack (if present)](#motor-control-pack-if-present)
  - [Power electronics pack (if present)](#power-electronics-pack-if-present)
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
| [`OutputFeedbackController`](inc/damp/concepts.hpp#L88) | Vector output-feedback controller: u = control(r, y), self-contained tick |
| [`ParameterEstimator`](inc/damp/concepts.hpp#L142) | Online grey-box parameter estimator: physical parameters + gating |
| [`SignalSource`](inc/damp/concepts.hpp#L127) | Self-clocked signal source: u = step(), finished when done() |
| [`SisoController`](inc/damp/concepts.hpp#L74) | Scalar output-feedback controller: u = control(r, y), fixed rate |
| [`StateEstimator`](inc/damp/concepts.hpp#L114) | State estimator: x̂ = estimate(y, u) — the fused per-tick form |
| [`StateFeedbackController`](inc/damp/concepts.hpp#L101) | State-feedback law: u = control(r, x), reference-first |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`minmax`](inc/damp/backend.hpp#L164) | Ordered {min, max} pair returned by value (+1 more overload) |

## Scalar math & complex

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`complex`](inc/damp/math/complex.hpp#L38) | Constexpr complex number class for compile-time computations |
| [`DCM`](inc/damp/math/geometry.hpp#L63) | Direction cosine matrix — 3×3 rotation (SO(3) wrapper over Mat3) |
| [`Euler`](inc/damp/math/geometry.hpp#L205) | Intrinsic Euler angles for sequence Order (default ZYX aerospace YPR) |
| [`EulerOrder`](inc/damp/math/geometry.hpp#L39) | Intrinsic Euler angle sequence (axis order of successive principal rotations) |
| [`Quaternion`](inc/damp/math/geometry.hpp#L381) | Unit quaternion rotation (w, x, y, z) (Hamilton product) |
| [`Transform4`](inc/damp/math/geometry.hpp#L826) | 4×4 homogeneous transform (SE(3) interop / DH export) |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`abs`](inc/damp/math/complex.hpp#L358) | Compute magnitude (absolute value) of a complex number (+1 more overload) |
| [`acos`](inc/damp/math/math.hpp#L150) | Arccosine ∈ [0, π]. Input is clamped to [−1, 1] in both paths |
| [`arg`](inc/damp/math/complex.hpp#L369) | Compute argument (phase angle) of a complex number |
| [`asin`](inc/damp/math/math.hpp#L131) | Arcsine ∈ [−π/2, π/2]. Input is clamped to [−1, 1] in both paths so behavior matches at compile and run time (std::asin would return NaN for \|x\| > 1) |
| [`atan`](inc/damp/math/math.hpp#L117) | Single-argument arctangent ∈ (−π/2, π/2) |
| [`atan2`](inc/damp/math/math.hpp#L105) | Two-argument arctangent, atan2(y, x) ∈ [−π, π] |
| [`cbrt`](inc/damp/math/math.hpp#L92) | Cube root (preserves sign for negative x) |
| [`ceil`](inc/damp/math/math.hpp#L329) | Ceiling — smallest integer ≥ x |
| [`copysign`](inc/damp/math/math.hpp#L406) | Copy sign — magnitude of mag with the sign of sgn_src |
| [`cos`](inc/damp/math/math.hpp#L169) | Cosine |
| [`db2mag`](inc/damp/math/math.hpp#L470) | Decibels to magnitude, 10^(db/20) |
| [`deg2rad`](inc/damp/math/math.hpp#L492) | Degrees to radians, deg·π/180 |
| [`exp`](inc/damp/math/math.hpp#L232) | Exponential function |
| [`finite_non_negative`](inc/damp/math/math.hpp#L448) | True if x is finite and non-negative (x ≥ 0) |
| [`finite_positive`](inc/damp/math/math.hpp#L438) | True if x is finite and strictly greater than zero |
| [`floor`](inc/damp/math/math.hpp#L317) | Floor — largest integer ≤ x |
| [`fmod`](inc/damp/math/math.hpp#L361) | Floating-point remainder, x − y·trunc(x/y) (sign of x), matching std::fmod's truncated-quotient convention |
| [`hypot`](inc/damp/math/math.hpp#L76) | Euclidean distance hypot(x, y) = √(x² + y²), without overflow |
| [`isfinite`](inc/damp/math/math.hpp#L425) | Finiteness test — false for NaN and ±∞ |
| [`log`](inc/damp/math/math.hpp#L249) | Natural logarithm |
| [`log10`](inc/damp/math/math.hpp#L380) | Base-10 logarithm, log10(x) = ln(x) / ln(10) |
| [`mag2db`](inc/damp/math/math.hpp#L459) | Magnitude to decibels, 20·log10(mag) |
| [`nearbyint`](inc/damp/math/math.hpp#L345) | Round to nearest integer; ties to even (IEEE default / `FE_TONEAREST`) |
| [`pow`](inc/damp/math/math.hpp#L268) | Power function, base^exponent (+1 more overload) |
| [`rad2deg`](inc/damp/math/math.hpp#L481) | Radians to degrees, rad·180/π |
| [`sgn`](inc/damp/math/math.hpp#L397) | Sign function — −1 if val < 0, 1 if val > 0, 0 if val == 0 |
| [`sin`](inc/damp/math/math.hpp#L182) | Sine |
| [`sincos`](inc/damp/math/math.hpp#L200) | Combined sine and cosine, {sin(x), cos(x)} |
| [`sqrt`](inc/damp/math/complex.hpp#L331) | Compute complex square root (constexpr) (+1 more overload) |
| [`tan`](inc/damp/math/math.hpp#L215) | Tangent |
| [`wrap`](inc/damp/math/math.hpp#L509) | Wrap x into the half-open interval [min, max) (period max − min) |
| [`wrap_pi`](inc/damp/math/math.hpp#L544) | Wrap an angle to [−π, π) |
| [`wrap_two_pi`](inc/damp/math/math.hpp#L562) | Wrap an angle to [0, 2π) |

## Linear algebra

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`Block`](inc/damp/matrix/block.hpp#L38) | Block view (non-owning) into a parent matrix |
| [`ColVec`](inc/damp/matrix/colvec.hpp#L27) | Concrete Column vector specialization of Matrix<N, 1, T> |
| [`ColView`](inc/damp/matrix/views.hpp#L257) | Non-owning column view of a matrix |
| [`Diagonal`](inc/damp/matrix/views.hpp#L49) | Diagonal view of a square matrix |
| [`EigenResult`](inc/damp/matrix/eigen.hpp#L34) | Eigenvalue computation result |
| [`FullQR`](inc/damp/matrix/decomposition.hpp#L273) | Result of a full (complete) QR factorization |
| [`is_matrix_element`](inc/damp/matrix/matrix_traits.hpp#L41) | True if T may be a Matrix element type |
| [`LowerTriangle`](inc/damp/matrix/views.hpp#L141) | Lower triangular view of a square matrix |
| [`Matrix`](inc/damp/matrix/core.hpp#L74) | Fixed-size, stack-allocated matrix for linear algebra operations |
| [`MatrixLike`](inc/damp/matrix/matrix_traits.hpp#L70) | Concept for any type that provides 2D matrix-like element access |
| [`MatrixLikeOf`](inc/damp/matrix/matrix_traits.hpp#L80) | Concept for a MatrixLike type with specific dimensions |
| [`NullSpace`](inc/damp/matrix/svd.hpp#L331) | Orthonormal basis for the null space (kernel) of a matrix |
| [`QRDecomposition`](inc/damp/matrix/decomposition.hpp#L191) | Thin QR factorization A = QR (modified Gram–Schmidt) |
| [`RowVec`](inc/damp/matrix/rowvec.hpp#L29) | Row vector specialization of Matrix<1, N, T> |
| [`RowView`](inc/damp/matrix/views.hpp#L165) | Non-owning row view of a matrix |
| [`SVDResult`](inc/damp/matrix/svd.hpp#L228) | Result of a full singular value decomposition A = U·Σ·Vᴴ |
| [`TransposeView`](inc/damp/matrix/views.hpp#L344) | Non-owning transpose view of a matrix (zero-copy) |
| [`UpperTriangle`](inc/damp/matrix/views.hpp#L114) | Upper triangular view of a square matrix |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`backward_substitute_transpose`](inc/damp/matrix/solve.hpp#L59) | Backward substitution for Lᵀx = b (real) / Lᴴx = b layout |
| [`cholesky`](inc/damp/matrix/decomposition.hpp#L82) | Cholesky decomposition for positive-definite matrices |
| [`cholesky_solve`](inc/damp/matrix/solve.hpp#L159) | Solve AX = B via Cholesky (A = LLᴴ) |
| [`compute_eigenvalues`](inc/damp/matrix/eigen.hpp#L382) | Compute the eigenvalues (and Schur vectors) of a real square matrix |
| [`cos`](inc/damp/matrix/functions.hpp#L738) | Matrix cosine via scaling and double-angle reconstruction |
| [`cosh`](inc/damp/matrix/functions.hpp#L809) | Matrix hyperbolic cosine cosh(A) = (exp(A) + exp(−A))/2 |
| [`default_tol`](inc/damp/matrix/matrix_traits.hpp#L86) | Type-appropriate default tolerance for floating-point comparisons. float  ~7 decimal digits  → 1e-6 double ~15 decimal digits → 1e-12 |
| [`det`](inc/damp/matrix/functions.hpp#L182) | Matrix determinant det(A) |
| [`expm`](inc/damp/matrix/functions.hpp#L323) | Matrix exponential via scaling and squaring with Padé approximant of degree 13 |
| [`forward_substitute`](inc/damp/matrix/solve.hpp#L37) | Forward substitution for Lx = b |
| [`frobenius_norm`](inc/damp/matrix/functions.hpp#L96) | Frobenius norm ‖A‖F = √(Σᵢⱼ \|aᵢⱼ\|²) |
| [`full_qr`](inc/damp/matrix/decomposition.hpp#L291) | Full QR factorization via Householder reflections (real or complex T) |
| [`infinity_norm`](inc/damp/matrix/functions.hpp#L35) | Infinity norm ‖A‖∞: maximum absolute row sum |
| [`inverse`](inc/damp/matrix/solve.hpp#L248) | Matrix inverse A⁻¹ via mat::solve(A, I) |
| [`log`](inc/damp/matrix/functions.hpp#L511) | Principal matrix logarithm via inverse scaling and squaring |
| [`logm`](inc/damp/matrix/functions.hpp#L550) | @brief MATLAB®-style alias for log |
| [`lu_decomposition`](inc/damp/matrix/decomposition.hpp#L131) | LU decomposition with partial pivoting |
| [`lu_solve`](inc/damp/matrix/solve.hpp#L187) | Solve AX = B via LU with partial pivoting |
| [`max`](inc/damp/matrix/core.hpp#L823) | Addition of two MatrixLike types (with broadcasting support) (+1 more overload) |
| [`null_space`](inc/damp/matrix/svd.hpp#L350) | Orthonormal basis for the null space {x : A·x = 0} via SVD |
| [`one_norm`](inc/damp/matrix/functions.hpp#L64) | One-norm ‖A‖₁: maximum absolute column sum |
| [`pow`](inc/damp/matrix/functions.hpp#L570) | Integer matrix power via binary exponentiation (+1 more overload) |
| [`pseudo_inverse`](inc/damp/matrix/svd.hpp#L303) | Moore–Penrose pseudoinverse A⁺ via SVD |
| [`qr_decompose`](inc/damp/matrix/decomposition.hpp#L223) | Thin QR via modified Gram–Schmidt |
| [`quadratic_form`](inc/damp/matrix/core.hpp#L924) | Symmetric congruence (quadratic) form  S = M X Mᵀ |
| [`rank`](inc/damp/matrix/functions.hpp#L259) | Matrix rank via Gaussian elimination with partial pivoting |
| [`rank_from_svd`](inc/damp/matrix/svd.hpp#L271) | Numerical rank from a precomputed SVD result |
| [`sin`](inc/damp/matrix/functions.hpp#L723) | Matrix sine via scaling and double-angle reconstruction |
| [`sincos`](inc/damp/matrix/functions.hpp#L653) | Compute sin(A) and cos(A) together via scaling and double-angle reconstruction |
| [`sinh`](inc/damp/matrix/functions.hpp#L793) | Matrix hyperbolic sine sinh(A) = (exp(A) − exp(−A))/2 |
| [`solve`](inc/damp/matrix/solve.hpp#L82) | Solve lower-triangular system LX = B via forward substitution (+2 more overloads) |
| [`sqrt`](inc/damp/matrix/functions.hpp#L446) | Matrix square root via Denman–Beavers iteration |
| [`sqrtm`](inc/damp/matrix/functions.hpp#L483) | @brief MATLAB®-style alias for sqrt |
| [`svd`](inc/damp/matrix/svd.hpp#L246) | Full singular value decomposition A = U·Σ·Vᴴ (one-sided Jacobi) |
| [`two_norm`](inc/damp/matrix/functions.hpp#L128) | Spectral norm ‖A‖₂ = σₘₐₓ(A) |

## LTI systems (SS / TF / ZPK / discretize)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`DiscretizationMethod`](inc/damp/systems/discretization.hpp#L26) | Discretization methods for continuous-time state-space systems |
| [`LeverrierResult`](inc/damp/systems/zpk.hpp#L500) | Faddeev–LeVerrier characteristic polynomial and adjoint coefficient matrices |
| [`PolyRootsResult`](inc/damp/systems/zpk.hpp#L89) | Roots of a real polynomial given in ascending powers |
| [`StateSpace`](inc/damp/systems/state_space.hpp#L125) | State-space representation for linear time-invariant systems (discrete or continuous) |
| [`StateSpaceZPKResult`](inc/damp/systems/zpk.hpp#L570) | SISO state-space → ZPK conversion result with runtime zero count |
| [`TransferFunction`](inc/damp/systems/transfer_function.hpp#L74) | SISO polynomial transfer function G(s) = num(s)/den(s) |
| [`ZPK`](inc/damp/systems/zpk.hpp#L177) | Zero-pole-gain (ZPK) representation of a SISO LTI system |
| [`ZPKResult`](inc/damp/systems/zpk.hpp#L399) | Result of converting a transfer function to zero-pole-gain form |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`discretize`](inc/damp/systems/discretization.hpp#L218) | Discretize a continuous-time state-space system |
| [`eval_frf`](inc/damp/systems/state_space.hpp#L166) | Evaluate frequency response of a state-space system |
| [`feedback`](inc/damp/systems/state_space.hpp#L311) | Negative feedback: y = sys1(u − sys2(y)) (+1 more overload) |
| [`minreal_zpk`](inc/damp/systems/zpk.hpp#L867) | Cancel matching pole-zero pairs on a ZPK model |
| [`parallel`](inc/damp/systems/state_space.hpp#L253) | Parallel connection (shared input, summed outputs) |
| [`poly_roots`](inc/damp/systems/zpk.hpp#L118) | Compute roots of an ascending-power real polynomial |
| [`series`](inc/damp/systems/state_space.hpp#L196) | Series connection: sys2 follows sys1 (u → sys1 → sys2 → y) (+1 more overload) |
| [`ss2tf`](inc/damp/systems/zpk.hpp#L843) | MATLAB®-style alias for SISO to_transfer_function |
| [`ss2zpk`](inc/damp/systems/zpk.hpp#L813) | MATLAB®-style alias for SISO to_zpk (state-space) |
| [`subtract`](inc/damp/systems/state_space.hpp#L394) | Differencing connection: outputs y = y₁ − y₂ |
| [`tf2zpk`](inc/damp/systems/zpk.hpp#L471) | MATLAB®-style alias for to_zpk (transfer function) |
| [`to_transfer_function`](inc/damp/systems/zpk.hpp#L826) | SISO state-space → transfer function (Leverrier), size NX+1 / NX+1 |
| [`to_zpk`](inc/damp/systems/zpk.hpp#L429) | Convert a SISO transfer function to zero-pole-gain form (+1 more overload) |
| [`zpk2tf`](inc/damp/systems/zpk.hpp#L480) | MATLAB®-style alias for ZPK::to_transfer_function |

> **Design-time — not for the PWM interrupt**

## Design-time synthesis (not PWM-rate)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`ActiveSetSolver`](inc/damp/design/qp.hpp#L618) | Default QP solver policy: the Goldfarb–Idnani active-set solve_qp() |
| [`AdmmSettings`](inc/damp/design/qp.hpp#L713) | Tuning parameters for the ADMM QP solver |
| [`AdmmSolver`](inc/damp/design/qp.hpp#L752) | ADMM (OSQP-style) QP solver policy — warm-started, fixed-cost iterations |
| [`ADRCResult`](inc/damp/controllers/adrc.hpp#L33) | Active Disturbance Rejection Control design result |
| [`BalancedRealizationResult`](inc/damp/design/model_reduction.hpp#L74) | Balanced realization result (Moore square-root method) |
| [`BalancedReductionResult`](inc/damp/design/model_reduction.hpp#L100) | Reduced-order model from balanced truncation or residualization |
| [`CommandProjectionResult`](inc/damp/controllers/action_governor.hpp#L185) | Result of an affine command projection |
| [`DiscretePIDResult`](inc/damp/controllers/pid.hpp#L43) | Fixed-rate discrete PID coefficients (canonical deploy form) |
| [`ESCResult`](inc/damp/controllers/esc.hpp#L103) | Design result for the extremum-seeking controller |
| [`HarmonicSuppressorResult`](inc/damp/controllers/harmonic_suppression.hpp#L45) | Design result for a multi-resonant harmonic suppressor |
| [`InductorCurrentPIResult`](inc/damp/design/pid_design.hpp#L835) | Fixed-rate inductor current-loop PI (topology-agnostic) |
| [`InteriorPointSolver`](inc/damp/design/qp.hpp#L1163) | Interior-point QP solver policy for damp::MPC |
| [`JordanBlock`](inc/damp/design/pole_placement.hpp#L831) | One Jordan mini-block of a desired closed-loop spectrum |
| [`JordanObjective`](inc/damp/design/pole_placement.hpp#L1165) | Robustness objective for place_jordan_optimal (the paper's two methods) |
| [`LeadLagResult`](inc/damp/controllers/lead_lag.hpp#L59) | Lead-lag compensator design result |
| [`LeadLagSeriesResult`](inc/damp/controllers/lead_lag.hpp#L204) | Cascaded lead+lag design result (2nd-order StateSpace + success) |
| [`LinearizationResult`](inc/damp/design/linearization.hpp#L41) | Result of nonlinear operating-point linearization |
| [`LQGAnalysisModels`](inc/damp/design/synthesis.hpp#L39) | Analysis-oriented models produced from an LQG design |
| [`LQGArtifacts`](inc/damp/design/synthesis.hpp#L89) | Synthesis artifact bundle: design + analysis models + runtime bundle |
| [`LQGIAnalysisModels`](inc/damp/design/synthesis.hpp#L152) | Analysis-oriented models produced from an LQGI design |
| [`LQGIArtifacts`](inc/damp/design/synthesis.hpp#L180) | Synthesis artifact bundle: LQGI servo design + analysis + ready-to-run controller |
| [`LQGIResult`](inc/damp/controllers/lqgi.hpp#L34) | LQGI design result |
| [`LQGPRArtifacts`](inc/damp/design/synthesis.hpp#L133) | Synthesis artifact bundle for SISO LQG + PR |
| [`LQGPRRuntimeBundle`](inc/damp/design/synthesis.hpp#L100) | Runtime bundle for SISO LQG + PR internal model compensation |
| [`LQGResult`](inc/damp/controllers/lqg.hpp#L46) | LQG design result |
| [`LQGRuntimeBundle`](inc/damp/design/synthesis.hpp#L50) | Runtime bundle for LQG control |
| [`LQIAnalysisModels`](inc/damp/design/synthesis.hpp#L143) | Analysis-oriented models produced from an LQI design |
| [`LQIArtifacts`](inc/damp/design/synthesis.hpp#L165) | Synthesis artifact bundle: LQI servo design + analysis + ready-to-run controller |
| [`LQIResult`](inc/damp/controllers/lqi.hpp#L48) | LQI design result |
| [`LQRResult`](inc/damp/controllers/lqr.hpp#L56) | Linear-Quadratic Regulator design result |
| [`MinimalRealizationResult`](inc/damp/design/minreal.hpp#L84) | Minimal realization result |
| [`ModelReductionMethod`](inc/damp/design/model_reduction.hpp#L58) | Method for eliminating states in modred |
| [`OptimalJordanPlacement`](inc/damp/design/pole_placement.hpp#L1172) | Result of optimized arbitrary pole placement (place_jordan_optimal) |
| [`PIDPerformanceSpec`](inc/damp/design/pid_design.hpp#L585) | Time-domain performance targets for quick PID synthesis |
| [`PIDResult`](inc/damp/controllers/pid.hpp#L97) | 2-DOF continuous-time PID controller design result |
| [`PIDType`](inc/damp/design/pid_design.hpp#L47) | PID controller type selection for tuning methods |
| [`PRResult`](inc/damp/controllers/pr.hpp#L72) | Proportional-Resonant controller design result |
| [`QPResult`](inc/damp/design/qp.hpp#L82) | Result of a dense QP solve |
| [`QPStatus`](inc/damp/design/qp.hpp#L51) | Termination status of a QP solve |
| [`RepetitiveResult`](inc/damp/controllers/repetitive.hpp#L134) | Design result for the repetitive controller |
| [`SMCResult`](inc/damp/controllers/smc.hpp#L32) | Tuning parameters for a first-order sliding-mode controller |
| [`SmithPredictorResult`](inc/damp/controllers/smith_predictor.hpp#L85) | Smith predictor design result (discrete PID + FO model + delay samples) |
| [`STSMCResult`](inc/damp/controllers/stsmc.hpp#L30) | Super-twisting (second-order sliding-mode) controller design result |
| [`WarmStartActiveSetSolver`](inc/damp/design/qp.hpp#L654) | Warm-started active-set solver policy — the damp::MPC default |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`ackermann`](inc/damp/design/pole_placement.hpp#L1582) | Single-input pole placement via Ackermann's formula |
| [`adrc`](inc/damp/controllers/adrc.hpp#L86) | Active Disturbance Rejection Control design |
| [`amigo_kappa_tau`](inc/damp/design/pid_design.hpp#L237) | AMIGO PI from ultimate gain/period and static gain Kₛ |
| [`balanced_realization`](inc/damp/design/model_reduction.hpp#L569) | Descriptive alias for balreal |
| [`balanced_truncation`](inc/damp/design/model_reduction.hpp#L586) | Descriptive alias for balred |
| [`balreal`](inc/damp/design/model_reduction.hpp#L383) | Balanced realization via the square-root (Moore / Laub) method |
| [`balred`](inc/damp/design/model_reduction.hpp#L471) | Balanced truncation to NR states |
| [`bandwidth_from_settling_time`](inc/damp/design/pid_design.hpp#L571) | Map settling-time and damping-ratio targets to a bandwidth estimate |
| [`build_lqg_analysis_models`](inc/damp/design/synthesis.hpp#L211) | Build analysis models from an LQG design |
| [`build_lqgi_analysis_models`](inc/damp/design/synthesis.hpp#L290) | Build analysis models from an LQGI servo design |
| [`build_lqi_analysis_models`](inc/damp/design/synthesis.hpp#L259) | Build analysis models from an LQI servo design |
| [`build_mpc_analysis_models`](inc/damp/controllers/mpc.hpp#L690) | Build the unconstrained-MPC LTI analysis models |
| [`care`](inc/damp/design/riccati.hpp#L763) | Solve the Continuous-time Algebraic Riccati Equation (CARE) |
| [`cbf_relative_degree_1`](inc/damp/controllers/action_governor.hpp#L299) | Build a relative-degree-1 CBF inequality row |
| [`closed_loop_poles`](inc/damp/design/stability.hpp#L335) | Compute closed-loop poles (eigenvalues) with state feedback |
| [`cohen_coon`](inc/damp/design/pid_design.hpp#L280) | Cohen-Coon tuning from first-order-plus-dead-time model |
| [`continuous_lqr`](inc/damp/controllers/lqr.hpp#L377) | Continuous-time Linear-Quadratic Regulator design (+1 more overload) |
| [`controllability_gramian`](inc/damp/design/stability.hpp#L118) | Continuous/discrete controllability Gramian W_c |
| [`controllability_matrix`](inc/damp/design/stability.hpp#L52) | Compute the controllability matrix [B, AB, A²B, ..., A^(N-1)B] |
| [`damping_ratio_from_overshoot_percent`](inc/damp/design/pid_design.hpp#L510) | Map percent overshoot target to equivalent damping ratio |
| [`dare`](inc/damp/design/riccati.hpp#L634) | Solve the Discrete Algebraic Riccati Equation (DARE) |
| [`discrete_lqg`](inc/damp/controllers/lqg.hpp#L120) | Discrete Linear-Quadratic-Gaussian regulator design |
| [`discrete_lqgi`](inc/damp/controllers/lqgi.hpp#L124) | Discrete LQG with integral action (LQI + Kalman) for output tracking |
| [`discrete_lqi`](inc/damp/controllers/lqi.hpp#L93) | Discrete Linear-Quadratic-Integral (LQI) design for output tracking |
| [`discrete_lqr`](inc/damp/controllers/lqr.hpp#L115) | Discrete-time Linear-Quadratic Regulator design |
| [`discrete_lqr_from_continuous`](inc/damp/controllers/lqr.hpp#L253) | Design discrete LQR from continuous-time system via discretization (+1 more overload) |
| [`discretize_lqr_cost`](inc/damp/controllers/lqr.hpp#L181) | Discretize a continuous LQR cost integral over one sample (Van Loan) |
| [`dlqr`](inc/damp/controllers/lqr.hpp#L298) | Discrete-time LQR design (MATLAB®-style short name) |
| [`dlyap`](inc/damp/design/lyapunov.hpp#L130) | Solve the discrete-time Lyapunov (Stein) equation A X Aᵀ − X + Q = 0 |
| [`esc`](inc/damp/controllers/esc.hpp#L153) | Synthesize an extremum-seeking controller |
| [`esc_lpf_alpha`](inc/damp/controllers/esc.hpp#L131) | First-order discrete low-pass coefficient for corner wc [rad/s] at Ts |
| [`esc_mppt`](inc/damp/controllers/esc.hpp#L200) | MPPT-flavored ESC: maximize a power measurement by perturbing the operating point (e.g. converter duty or reference voltage) |
| [`hankel_singular_values`](inc/damp/design/model_reduction.hpp#L596) | Descriptive alias for hankelsv |
| [`hankelsv`](inc/damp/design/model_reduction.hpp#L439) | Hankel singular values of a stable state-space system |
| [`harmonic_suppressor`](inc/damp/controllers/harmonic_suppression.hpp#L80) | Synthesize a multi-resonant harmonic suppressor |
| [`is_closed_loop_stable_discrete`](inc/damp/design/stability.hpp#L245) | Check closed-loop stability for discrete system with state feedback |
| [`is_controllable`](inc/damp/design/stability.hpp#L176) | Check if a system is controllable |
| [`is_observable`](inc/damp/design/stability.hpp#L194) | Check if a system is observable |
| [`is_stabilizable`](inc/damp/design/riccati.hpp#L44) | Check if (A, B) is a stabilizable pair |
| [`is_stable_discrete`](inc/damp/design/stability.hpp#L215) | Check if a discrete-time system matrix A is stable |
| [`lag`](inc/damp/controllers/lead_lag.hpp#L181) | Design a lag compensator from desired low-frequency gain boost |
| [`lambda_tuning`](inc/damp/design/pid_design.hpp#L398) | Lambda tuning for FOPDT model |
| [`lead`](inc/damp/controllers/lead_lag.hpp#L132) | Design a lead compensator from desired phase boost at a target frequency |
| [`lead_lag`](inc/damp/controllers/lead_lag.hpp#L234) | Design a lead-lag compensator (cascade of lead + lag sections) |
| [`lead_lag_direct`](inc/damp/controllers/lead_lag.hpp#L264) | Direct lead-lag specification from zero/pole locations |
| [`linearize`](inc/damp/design/linearization.hpp#L135) | Linearize nonlinear dynamics and output maps about an operating point (+1 more overload) |
| [`lqg_bundle`](inc/damp/design/synthesis.hpp#L322) | Synthesize the full LQG artifact bundle in one call |
| [`lqg_from_parts`](inc/damp/controllers/lqg.hpp#L149) | Assemble an LQG design from separately computed Kalman and LQR results |
| [`lqg_pr_bundle`](inc/damp/design/synthesis.hpp#L349) | Synthesize a SISO LQG + PR design with internal-model compensation |
| [`lqgi_bundle`](inc/damp/design/synthesis.hpp#L396) | Synthesize the full LQGI servo artifact bundle in one call |
| [`lqi_bundle`](inc/damp/design/synthesis.hpp#L378) | Synthesize the full LQI servo artifact bundle in one call |
| [`lqr_gain`](inc/damp/design/riccati.hpp#L714) | Optimal LQR state-feedback gain from a Riccati solution |
| [`lqrd`](inc/damp/controllers/lqr.hpp#L317) | Sampled-data LQR from continuous plant (MATLAB®-style short name) (+1 more overload) |
| [`lyap`](inc/damp/design/lyapunov.hpp#L105) | Solve the continuous-time Lyapunov equation A X + X Aᵀ + Q = 0 |
| [`minimal_realization`](inc/damp/design/minreal.hpp#L428) | Descriptive alias for minreal |
| [`minreal`](inc/damp/design/minreal.hpp#L284) | Minimal realization — cancel uncontrollable and unobservable modes |
| [`modred`](inc/damp/design/model_reduction.hpp#L518) | Model reduction by truncation or DC-matched residualization |
| [`MPC`](inc/damp/controllers/mpc.hpp#L907) | Deduce the runtime from its artifacts: MPC controller{art}; (+1 more overload) |
| [`mpc`](inc/damp/controllers/offset_free_mpc.hpp#L132) | Synthesize an offset-free constrained MPC (controller + estimator) |
| [`observability_gramian`](inc/damp/design/stability.hpp#L146) | Continuous/discrete observability Gramian W_o |
| [`observability_matrix`](inc/damp/design/stability.hpp#L82) | Compute the observability matrix [C; CA; CA²; ...; CA^(N-1)] |
| [`OffsetFreeMPC`](inc/damp/controllers/offset_free_mpc.hpp#L258) | Deduce the runtime from its artifacts: OffsetFreeMPC controller{art}; |
| [`phase_margin_from_damping_ratio`](inc/damp/design/pid_design.hpp#L540) | Approximate phase margin from damping ratio |
| [`pi_pole_placement_first_order`](inc/damp/design/pid_design.hpp#L769) | PI gains that place the closed-loop poles of a first-order plant (+1 more overload) |
| [`pid`](inc/damp/controllers/pid.hpp#L222) | 2-DOF continuous PID controller design |
| [`pid_from_bandwidth`](inc/damp/design/pid_design.hpp#L432) | Design PID from desired bandwidth and phase margin |
| [`pid_from_performance_spec`](inc/damp/design/pid_design.hpp#L603) | Design PID directly from settling-time and overshoot targets |
| [`pid_pole_placement`](inc/damp/design/pid_design.hpp#L640) | Direct PID pole placement for a first-order-plus-dead-time model (+1 more overload) |
| [`place`](inc/damp/design/pole_placement.hpp#L88) | Robust multi-input pole placement (Kautsky–Nichols–Van Dooren, real poles) (+5 more overloads) |
| [`place_discrete`](inc/damp/design/pole_placement.hpp#L581) | Discrete plant + z-plane poles → discrete \(K\) (+3 more overloads) |
| [`place_jordan`](inc/damp/design/pole_placement.hpp#L1138) | Exact pole placement with an arbitrary Jordan structure (Schmid–Ntogramatzidis–Nguyen–Pandey / Klein–Moore parametric form) |
| [`place_jordan_optimal`](inc/damp/design/pole_placement.hpp#L1237) | Robust / minimum-gain arbitrary pole placement (Schmid et al., Methods 1–2) |
| [`place_observer`](inc/damp/design/pole_placement.hpp#L659) | Continuous observer gain: place eigenvalues of \(A - LC\) (s-plane poles) (+5 more overloads) |
| [`place_observer_discrete`](inc/damp/design/pole_placement.hpp#L769) | Discrete plant + z-plane poles → discrete observer gain \(L\) (+3 more overloads) |
| [`pr`](inc/damp/controllers/pr.hpp#L153) | Design a Proportional-Resonant controller |
| [`pr_harmonics`](inc/damp/controllers/pr.hpp#L185) | Design multiple-harmonic PR controller gains |
| [`project_affine`](inc/damp/controllers/action_governor.hpp#L231) | Project u_des onto the polyhedron A u ≤ b (Euclidean) |
| [`project_box`](inc/damp/controllers/action_governor.hpp#L78) | Euclidean projection of u onto the axis-aligned box [umin, umax] (+1 more overload) |
| [`rank`](inc/damp/design/stability.hpp#L163) | Compute rank of a matrix via Gaussian elimination with partial pivoting |
| [`repetitive`](inc/damp/controllers/repetitive.hpp#L171) | Synthesize a repetitive controller with a scalar robustness filter Q |
| [`repetitive_binomial`](inc/damp/controllers/repetitive.hpp#L223) | Synthesize a repetitive controller with a binomial zero-phase FIR Q |
| [`simc`](inc/damp/design/pid_design.hpp#L338) | SIMC (Skogestad Internal Model Control) tuning for FOPDT models |
| [`smc`](inc/damp/controllers/smc.hpp#L59) | Bundle hand-picked SMC parameters into an SMCResult |
| [`smith_predictor_from_fopdt`](inc/damp/controllers/smith_predictor.hpp#L141) | Design a Smith predictor from FOPDT parameters |
| [`smith_predictor_from_pid`](inc/damp/controllers/smith_predictor.hpp#L207) | Pack a Smith predictor from an existing continuous PID and FOPDT plant |
| [`solve_miqp`](inc/damp/design/qp.hpp#L1215) | Mixed-integer QP by branch and bound over the active-set relaxation |
| [`solve_qp`](inc/damp/design/qp.hpp#L525) | Solve a strictly convex inequality-constrained QP (dual active-set) (+1 more overload) |
| [`solve_qp_admm`](inc/damp/design/qp.hpp#L916) | One-shot ADMM QP solve (cold start) |
| [`solve_qp_interior_point`](inc/damp/design/qp.hpp#L960) | Interior-point QP solver (primal-dual path following) |
| [`stability_margin_continuous`](inc/damp/design/stability.hpp#L268) | Compute stability margin for continuous system |
| [`stability_margin_discrete`](inc/damp/design/stability.hpp#L297) | Compute stability margin for discrete system |
| [`state_mpc`](inc/damp/controllers/mpc.hpp#L502) | Synthesize a constrained linear MPC (condensed dense QP, Δu form) (+1 more overload) |
| [`stsmc`](inc/damp/controllers/stsmc.hpp#L87) | Synthesize super-twisting gains from a disturbance-derivative bound |
| [`stsmc_gains`](inc/damp/controllers/stsmc.hpp#L123) | Super-twisting controller from gains you specify directly |
| [`suggest_mpc_horizon`](inc/damp/controllers/mpc.hpp#L585) | Suggest MPC horizons from an explicit settling-time target (+1 more overload) |
| [`tyreus_luyben`](inc/damp/design/pid_design.hpp#L175) | Tyreus-Luyben tuning from ultimate gain and ultimate period |
| [`unbounded_bound`](inc/damp/design/qp.hpp#L67) | Sentinel bound treated as "no constraint" on that row |
| [`ziegler_nichols`](inc/damp/design/pid_design.hpp#L75) | Ziegler-Nichols tuning from ultimate gain and ultimate period |
| [`ziegler_nichols_step`](inc/damp/design/pid_design.hpp#L124) | Ziegler-Nichols step response method (reaction curve) |

> **Runtime — tick objects and signal processing**

## Runtime controllers

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`ActionGovernor`](inc/damp/controllers/action_governor.hpp#L331) | Runtime action governor — QP projection of u_des onto A u ≤ b |
| [`ADRCController`](inc/damp/controllers/adrc.hpp#L148) | Active Disturbance Rejection Control (ADRC) |
| [`BoxCommandFilter`](inc/damp/controllers/action_governor.hpp#L119) | Runtime box command filter — closed-form clamp per tick |
| [`Cascade`](inc/damp/controllers/composition.hpp#L49) | Series cascade of two SISO controllers: outer → inner reference |
| [`CBFConstraint`](inc/damp/controllers/action_governor.hpp#L283) | One affine row from a relative-degree-1 CBF condition |
| [`ContinuousPID`](inc/damp/controllers/pid.hpp#L633) | Continuous-gain PID with per-tick sample time (variable-rate secondary form) |
| [`ESCConfig`](inc/damp/controllers/esc.hpp#L59) | Extremum-seeking controller configuration (discrete realization) |
| [`ExtremumSeekingController`](inc/damp/controllers/esc.hpp#L238) | Extremum-seeking controller runtime (model-free online optimizer) |
| [`ExtremumType`](inc/damp/controllers/esc.hpp#L47) | Whether ESC climbs to a maximum or descends to a minimum of the objective |
| [`HarmonicSuppressor`](inc/damp/controllers/harmonic_suppression.hpp#L121) | Multi-resonant harmonic suppressor — a parallel bank of PR resonators |
| [`LeadLagController`](inc/damp/controllers/lead_lag.hpp#L281) | Discrete lead-lag compensator |
| [`LQG`](inc/damp/controllers/lqg.hpp#L176) | Linear-Quadratic-Gaussian (LQG) controller |
| [`LQGI`](inc/damp/controllers/lqgi.hpp#L158) | Linear-Quadratic-Gaussian-Integral (LQGI) controller |
| [`LQI`](inc/damp/controllers/lqi.hpp#L155) | Linear-Quadratic-Integral (LQI) controller |
| [`LQRCost`](inc/damp/controllers/lqr.hpp#L142) | Discretized LQR cost weights (Q, R, N) for a sampled-data problem |
| [`MPC`](inc/damp/controllers/mpc.hpp#L770) | Runtime constrained MPC controller (fixed per-tick iteration budget) |
| [`MPCAnalysisModels`](inc/damp/controllers/mpc.hpp#L666) | LTI models of the unconstrained MPC loop, for margin/robustness analysis |
| [`MPCArtifacts`](inc/damp/controllers/mpc.hpp#L235) | Condensed-QP data produced by mpc(), consumed by damp::MPC |
| [`MPCConstraints`](inc/damp/controllers/mpc.hpp#L197) | Box constraints for mpc() |
| [`MPCHorizonSuggestion`](inc/damp/controllers/mpc.hpp#L561) | Advisory NP/NC horizon values from suggest_mpc_horizon() |
| [`MPCWeights`](inc/damp/controllers/mpc.hpp#L162) | Cost weights for mpc() |
| [`MultiPRController`](inc/damp/controllers/pr.hpp#L369) | Multi-harmonic PR Controller |
| [`OffsetFreeMPC`](inc/damp/controllers/offset_free_mpc.hpp#L196) | Runtime offset-free MPC: constrained MPC + disturbance-augmented Kalman filter |
| [`OffsetFreeMPCArtifacts`](inc/damp/controllers/offset_free_mpc.hpp#L84) | Combined MPC + disturbance-augmented Kalman design, consumed by damp::OffsetFreeMPC |
| [`PIDController`](inc/damp/controllers/pid.hpp#L312) | Fixed-rate discrete 2-DOF PID (canonical runtime) (+2 more overloads) |
| [`PIDMode`](inc/damp/controllers/pid.hpp#L247) | Compile-time selection of the PID control-law structure |
| [`PIDRuntimeMode`](inc/damp/controllers/pid.hpp#L271) | Runtime operating mode for PIDController / ContinuousPID |
| [`PRController`](inc/damp/controllers/pr.hpp#L222) | Discrete Proportional-Resonant Controller |
| [`RepetitiveConfig`](inc/damp/controllers/repetitive.hpp#L70) | Repetitive-controller tuning + period (with optional zero-phase FIR Q) |
| [`RepetitiveController`](inc/damp/controllers/repetitive.hpp#L273) | Plug-in repetitive controller runtime (fixed-size internal model) |
| [`SMCController`](inc/damp/controllers/smc.hpp#L136) | First-order sliding-mode controller (SMC) for a SISO plant |
| [`SmithPredictor`](inc/damp/controllers/smith_predictor.hpp#L267) | Smith predictor runtime: SISO primary + FO model + pure delay |
| [`StateFeedback`](inc/damp/controllers/lqr.hpp#L431) | Runtime full-state feedback law u = −Kx |
| [`STSMCController`](inc/damp/controllers/stsmc.hpp#L179) | Super-twisting controller (second-order sliding mode) |
| [`SwitchedController`](inc/damp/controllers/composition.hpp#L117) | Bumpless-ish switch between normal, experiment, and backup SISO laws |
| [`SwitchMode`](inc/damp/controllers/composition.hpp#L94) | Which path owns the plant command |

## Observers & estimators

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AdaptationGate`](inc/damp/estimation/parameter_estimation.hpp#L315) | Gated-adaptation policy: decide whether a parameter update is safe to apply |
| [`Chirp`](inc/damp/estimation/excitation/chirp.hpp#L132) | Linear or logarithmic chirp runtime generator |
| [`ChirpConfig`](inc/damp/estimation/excitation/chirp.hpp#L40) | Configuration for a sine chirp excitation |
| [`ChirpMode`](inc/damp/estimation/excitation/chirp.hpp#L26) | Chirp sweep law |
| [`ChirpResult`](inc/damp/estimation/excitation/chirp.hpp#L83) | Chirp design payload |
| [`ClassicalDOB`](inc/damp/estimation/dob.hpp#L370) | Classical Pn^-1·Q disturbance observer runtime (bolt-on compensator) |
| [`ClassicalDobResult`](inc/damp/estimation/dob.hpp#L295) | Design result for the classical Pn^-1·Q disturbance observer |
| [`CompensatorKind`](inc/damp/estimation/successive_compensator.hpp#L38) | Which biquad family to load for a feature |
| [`ComplementaryFilter`](inc/damp/estimation/sensor_fusion.hpp#L64) | Simple complementary filter for orientation estimation |
| [`DOB`](inc/damp/estimation/dob.hpp#L153) | Lightweight SISO disturbance observer runtime |
| [`DOBConfig`](inc/damp/estimation/dob.hpp#L31) | Configuration for a first-order disturbance observer |
| [`EKFMeasFn`](inc/damp/estimation/ekf.hpp#L82) | Concept for EKF measurement functions |
| [`EKFStateFn`](inc/damp/estimation/ekf.hpp#L54) | Concept for EKF state functions |
| [`ErrorStateJacobian`](inc/damp/estimation/eskf.hpp#L247) | ESKF prediction Jacobians F and G (nominal state updated outside) |
| [`ErrorStateKalmanFilter`](inc/damp/estimation/eskf.hpp#L294) | Runtime error-state KF: tracks δx and P, not the full nominal state |
| [`ESKFMeasFn`](inc/damp/estimation/eskf.hpp#L270) | Callable that returns a measurement linearization for one update |
| [`ESKFOrientationFilter`](inc/damp/estimation/sensor_fusion.hpp#L399) | Turnkey 6-state attitude ESKF: owns q, b_g, and the error filter |
| [`ESKFPredictFn`](inc/damp/estimation/eskf.hpp#L259) | Callable that returns ErrorStateJacobian for one predict step |
| [`ESKFResult`](inc/damp/estimation/eskf.hpp#L50) | Design payload for an ESKF: Q, R, P₀, and success |
| [`ExperimentSafety`](inc/damp/estimation/experiment_safety.hpp#L62) | Experiment lifecycle + command protection for on-target commissioning |
| [`ExperimentSafetyConfig`](inc/damp/estimation/experiment_safety.hpp#L39) | Configuration for ExperimentSafety |
| [`ExtendedKalmanFilter`](inc/damp/estimation/ekf.hpp#L112) | Extended Kalman Filter for nonlinear discrete-time systems |
| [`FirstOrderPlantEstimator`](inc/damp/estimation/parameter_estimation.hpp#L126) | Online estimator for a first-order plant's gain and time constant |
| [`FirstOrderPlantEstimatorConfig`](inc/damp/estimation/parameter_estimation.hpp#L85) | Configuration for the first-order grey-box estimator |
| [`FrequencyResponseEstimator`](inc/damp/estimation/frequency_response.hpp#L106) | Lock-in FRF estimator over a fixed table of NFreq bins |
| [`FrfFeatureKind`](inc/damp/estimation/frequency_response.hpp#L373) | Peak (resonance) or valley (anti-resonance) on an FRF magnitude curve |
| [`FrfFeatureResult`](inc/damp/estimation/frequency_response.hpp#L453) | Combined peak + valley reduction of an FRF table |
| [`FrfMargins`](inc/damp/estimation/frequency_response.hpp#L485) | Open-loop margins read from a discrete FRF table (on-target) |
| [`FrfPidAutotuneConfig`](inc/damp/estimation/commissioning_frontends.hpp#L109) | Configuration for FRF-based PID gain selection |
| [`FrfPidAutotuneResult`](inc/damp/estimation/commissioning_frontends.hpp#L125) | FRF PID autotune hand-off (gains + success) |
| [`FrfPoint`](inc/damp/estimation/frequency_response.hpp#L60) | One measured frequency-response point (on-target FRF table entry) |
| [`Impulse`](inc/damp/estimation/excitation/impulse.hpp#L86) | Rectangular impulse (finite-width Dirac stand-in) |
| [`ImpulseConfig`](inc/damp/estimation/excitation/impulse.hpp#L32) | Configuration for a one-shot rectangular impulse (Dirac stand-in) |
| [`ImpulseResult`](inc/damp/estimation/excitation/impulse.hpp#L52) | Impulse design payload |
| [`ImuMount`](inc/damp/estimation/serial_arm_pose.hpp#L88) | Fixed extrinsic of an IMU relative to its DH / link frame |
| [`ImuSample`](inc/damp/estimation/ins_mechanization.hpp#L125) | IMU sample in the body frame |
| [`InsEskfResult`](inc/damp/estimation/ins_eskf.hpp#L180) | Design result for the 15-state INS error-state filter |
| [`InsNavigator`](inc/damp/estimation/ins_eskf.hpp#L754) | Runtime strapdown navigator: nominal InsState + 15-state ESKF |
| [`InsState`](inc/damp/estimation/ins_mechanization.hpp#L139) | Nominal strapdown navigation state (external to the ESKF error state) |
| [`KalmanFilter`](inc/damp/estimation/kalman.hpp#L249) | Runtime Kalman filter for embedded systems |
| [`KalmanResult`](inc/damp/estimation/kalman.hpp#L45) | Steady-state Kalman filter design result |
| [`LinkAttitudeEskf`](inc/damp/estimation/serial_arm_pose.hpp#L143) | One-frame attitude for undercarriage / house (accel-only ESKF) |
| [`Luenberger`](inc/damp/estimation/luenberger.hpp#L400) | Luenberger state observer (runtime) |
| [`LuenbergerResult`](inc/damp/estimation/luenberger.hpp#L92) | Luenberger observer design result |
| [`MadgwickFilter`](inc/damp/estimation/sensor_fusion.hpp#L105) | Madgwick gradient-descent AHRS filter |
| [`MahonyFilter`](inc/damp/estimation/sensor_fusion.hpp#L163) | Mahony nonlinear complementary filter with PI correction |
| [`MeasJacobian`](inc/damp/estimation/ekf.hpp#L69) | Measurement prediction result from the user's observation function |
| [`MHE`](inc/damp/estimation/mhe.hpp#L283) | Runtime moving-horizon estimator (fixed per-tick iteration budget) (+1 more overload) |
| [`MHEArtifacts`](inc/damp/estimation/mhe.hpp#L110) | Window QP data produced by mhe(), consumed by damp::MHE |
| [`MHEConstraints`](inc/damp/estimation/mhe.hpp#L87) | Box constraints on the estimated states (applied to every window state) |
| [`ModeExtractorConfig`](inc/damp/estimation/frequency_response.hpp#L411) | Configuration for FRF peak / valley extraction |
| [`ModeExtractorResult`](inc/damp/estimation/frequency_response.hpp#L428) | Result of reducing an FRF table to at most MaxModes peaks (resonances) |
| [`MultiSine`](inc/damp/estimation/excitation/multi_sine.hpp#L200) | Sum-of-tones multi-sine runtime generator |
| [`MultiSineConfig`](inc/damp/estimation/excitation/multi_sine.hpp#L62) | Configuration for fixed-component multi-sine excitation |
| [`MultiSineResult`](inc/damp/estimation/excitation/multi_sine.hpp#L104) | Multi-sine design payload |
| [`NavFrame`](inc/damp/estimation/ins_mechanization.hpp#L77) | Local-level navigation frame for strapdown mechanization |
| [`ParameterDriftMonitor`](inc/damp/estimation/parameter_estimation.hpp#L252) | Residual-based drift detector: has the plant moved away from the model? |
| [`PositionAxes`](inc/damp/estimation/serial_arm_pose.hpp#L57) | Which world Cartesian components of a frame origin to fuse |
| [`PRBS`](inc/damp/estimation/excitation/prbs.hpp#L181) | Maximal-length PRBS runtime generator |
| [`PRBSConfig`](inc/damp/estimation/excitation/prbs.hpp#L81) | Configuration for maximal-length pseudo-random binary excitation |
| [`PRBSResult`](inc/damp/estimation/excitation/prbs.hpp#L125) | PRBS design payload |
| [`Ramp`](inc/damp/estimation/excitation/ramp.hpp#L110) | Rate-limited ramp runtime generator |
| [`RampConfig`](inc/damp/estimation/excitation/ramp.hpp#L32) | Configuration for a slew-rate-limited ramp excitation |
| [`RampResult`](inc/damp/estimation/excitation/ramp.hpp#L64) | Ramp design payload |
| [`ReducedLuenberger`](inc/damp/estimation/luenberger.hpp#L507) | Reduced-order (Gopinath) state observer (runtime) |
| [`ReducedLuenbergerResult`](inc/damp/estimation/luenberger.hpp#L217) | Reduced-order (Gopinath) observer design result |
| [`RelayAutotuneConfig`](inc/damp/estimation/relay_autotune.hpp#L120) | Configuration for the relay-feedback autotuning experiment |
| [`RelayAutotuneOutput`](inc/damp/estimation/relay_autotune.hpp#L238) | Per-tick output of RelayAutotuner::step |
| [`RelayAutotuner`](inc/damp/estimation/relay_autotune.hpp#L257) | Runtime relay-feedback autotuner |
| [`RelayAutotuneResult`](inc/damp/estimation/relay_autotune.hpp#L180) | Relay-autotuner design payload |
| [`RelayAutotuneStatus`](inc/damp/estimation/relay_autotune.hpp#L222) | Lifecycle state of a RelayAutotuner |
| [`ResonanceAutotuneResult`](inc/damp/estimation/commissioning_frontends.hpp#L44) | Result of resonance / anti-resonance compensator design from an FRF |
| [`ResonantMode`](inc/damp/estimation/frequency_response.hpp#L384) | One extracted FRF extremum (peak or valley) |
| [`Rls`](inc/damp/estimation/rls.hpp#L193) | Scalar runtime RLS estimator |
| [`RlsConfig`](inc/damp/estimation/rls.hpp#L40) | Common RLS configuration |
| [`RlsResult`](inc/damp/estimation/rls.hpp#L82) | Scalar RLS design payload |
| [`RlsState`](inc/damp/estimation/rls.hpp#L68) | Scalar RLS runtime state |
| [`RlsVector`](inc/damp/estimation/rls.hpp#L287) | Vector runtime RLS estimator (NP parameters) |
| [`RlsVectorResult`](inc/damp/estimation/rls.hpp#L125) | Vector RLS design payload for N parameters |
| [`RlsVectorState`](inc/damp/estimation/rls.hpp#L111) | Vector RLS runtime state for N parameters |
| [`SerialArmPoseFilter`](inc/damp/estimation/serial_arm_pose.hpp#L353) | Joint-centric arm pose: optional base attitude + peel gyros → q + FK |
| [`SerialArmPoseResult`](inc/damp/estimation/serial_arm_pose.hpp#L225) | Design payload for the joint-space arm filter (+ optional attitude) |
| [`StateJacobian`](inc/damp/estimation/ekf.hpp#L41) | State prediction result from the user's dynamics function |
| [`SteadyStateKalmanFilter`](inc/damp/estimation/kalman.hpp#L426) | Steady-state (fixed-gain) Kalman estimator for LQG-class designs |
| [`Step`](inc/damp/estimation/excitation/step.hpp#L108) | Single step with optional \|dy/dt\| settle detection |
| [`StepConfig`](inc/damp/estimation/excitation/step.hpp#L35) | Configuration for a single step with optional settle detection |
| [`SteppedSine`](inc/damp/estimation/excitation/stepped_sine.hpp#L125) | Stepped-sine excitation — one pure tone at a time across a frequency table |
| [`SteppedSineConfig`](inc/damp/estimation/excitation/stepped_sine.hpp#L33) | Configuration for single-tone stepped-sine excitation |
| [`SteppedSineResult`](inc/damp/estimation/excitation/stepped_sine.hpp#L73) | Stepped-sine design payload |
| [`StepResult`](inc/damp/estimation/excitation/step.hpp#L67) | Single-step design payload |
| [`StepTrain`](inc/damp/estimation/excitation/step_train.hpp#L110) | Alternating +/- step train runtime generator |
| [`StepTrainConfig`](inc/damp/estimation/excitation/step_train.hpp#L32) | Configuration for alternating +/- step excitation |
| [`StepTrainResult`](inc/damp/estimation/excitation/step_train.hpp#L64) | Step-train design payload |
| [`SuccessiveCompensatorCommissioner`](inc/damp/estimation/successive_compensator.hpp#L208) | Successive biquad compensator bank (notch and/or band-pass) |
| [`SuccessiveCompensatorConfig`](inc/damp/estimation/successive_compensator.hpp#L50) | Configuration for successive compensator commissioning |
| [`Tone`](inc/damp/estimation/excitation/multi_sine.hpp#L29) | One sinusoidal component in a multi-sine excitation |
| [`UKFMeasFn`](inc/damp/estimation/ukf.hpp#L65) | Concept for UKF measurement functions |
| [`UKFStateFn`](inc/damp/estimation/ukf.hpp#L54) | Concept for UKF state (process) functions |
| [`UnscentedKalmanFilter`](inc/damp/estimation/ukf.hpp#L126) | Unscented (sigma-point) Kalman Filter for nonlinear discrete-time systems |
| [`UnscentedParams`](inc/damp/estimation/ukf.hpp#L81) | Tuning parameters for the scaled unscented transform |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`arm_frames_world`](inc/damp/estimation/serial_arm_pose.hpp#L848) | World-frame DH poses: fixed chain origin, base body→world R0 |
| [`chirp`](inc/damp/estimation/excitation/chirp.hpp#L115) | Build a chirp design payload from a configuration |
| [`classical_dob`](inc/damp/estimation/dob.hpp#L338) | Synthesize a classical disturbance observer from a nominal plant and Q-filter |
| [`compensator_from_feature`](inc/damp/estimation/successive_compensator.hpp#L179) | Map an FRF feature to RBJ biquad coefficients |
| [`discrete_lqe`](inc/damp/estimation/kalman.hpp#L175) | Discrete linear-quadratic estimator (LQE) — dual-of-LQR spelling of kalman (+1 more overload) |
| [`dob`](inc/damp/estimation/dob.hpp#L128) | Validate and package DOB configuration into a runtime-ready design result |
| [`enu_from_ned`](inc/damp/estimation/ins_mechanization.hpp#L116) | Map a NED vector into ENU (axis permute) |
| [`eskf_imu`](inc/damp/estimation/eskf.hpp#L179) | Design Q, R, P₀ for a 6-state IMU attitude ESKF (gyro + accel) |
| [`eskf_marg`](inc/damp/estimation/eskf.hpp#L218) | Design Q, R, P₀ for a 6-state MARG attitude ESKF (gyro + accel + mag) |
| [`extract_frf_features`](inc/damp/estimation/frequency_response.hpp#L740) | Extract peaks and valleys in one pass over the FRF table |
| [`extract_modes`](inc/damp/estimation/frequency_response.hpp#L715) | Extract resonant peaks (local \|G\| maxima) from an FRF table |
| [`extract_valleys`](inc/damp/estimation/frequency_response.hpp#L729) | Extract anti-resonance valleys (local \|G\| minima) from an FRF table |
| [`finish_extremum`](inc/damp/estimation/frequency_response.hpp#L561) | Half-power (peaks) or double-power (valleys) bandwidth → ζ |
| [`fit_second_order_step`](inc/damp/estimation/successive_compensator.hpp#L85) | Fit one underdamped second-order mode from a step / ring-down capture |
| [`frf_pid_autotune`](inc/damp/estimation/commissioning_frontends.hpp#L152) | PID gains from an open-loop FRF table + margin / bandwidth targets |
| [`gravity_nav`](inc/damp/estimation/ins_mechanization.hpp#L93) | Gravity vector in the local-level frame [m/s²] |
| [`heading_from_baseline_nav`](inc/damp/estimation/ins_eskf.hpp#L536) | Heading [rad] from a nav-frame baseline vector (dual-antenna difference) |
| [`impulse`](inc/damp/estimation/excitation/impulse.hpp#L70) | Build an impulse design payload |
| [`ins_aid_position`](inc/damp/estimation/ins_eskf.hpp#L503) | Predicted antenna / marker position with body lever-arm |
| [`ins_apply_correction`](inc/damp/estimation/ins_eskf.hpp#L489) | Inject filter correction into x and reset the error state |
| [`ins_error_jacobian`](inc/damp/estimation/ins_eskf.hpp#L307) | Discrete error-state transition and noise input for one IMU step |
| [`ins_eskf_design`](inc/damp/estimation/ins_eskf.hpp#L219) | Build Q / R / P0 for a 15-state INS ESKF from IMU noise densities |
| [`ins_heading`](inc/damp/estimation/ins_eskf.hpp#L521) | Horizontal heading [rad] of a body axis expressed in the nav frame |
| [`ins_inject`](inc/damp/estimation/ins_eskf.hpp#L433) | Inject δx into the nominal INS state (right-multiplicative q) |
| [`ins_predict`](inc/damp/estimation/ins_eskf.hpp#L466) | Mechanize nominal state and propagate the error covariance |
| [`ins_propagate_covariance`](inc/damp/estimation/ins_eskf.hpp#L403) | P ← F P Fᵀ + Q using sparse F (exact for ins_error_jacobian structure) |
| [`ins_update_heading`](inc/damp/estimation/ins_eskf.hpp#L638) | Sparse dual-antenna / yaw heading update (scalar) |
| [`ins_update_orientation`](inc/damp/estimation/ins_eskf.hpp#L676) | Sparse orientation update from a measured body→nav quaternion |
| [`ins_update_pose`](inc/damp/estimation/ins_eskf.hpp#L701) | Generic pose alias: position then orientation (two sparse updates) |
| [`ins_update_pose_heading`](inc/damp/estimation/ins_eskf.hpp#L720) | Pose alias: position + dual-antenna heading (two sparse updates) |
| [`ins_update_position`](inc/damp/estimation/ins_eskf.hpp#L551) | Sparse position update: y = p + Rℓ (+1 more overload) |
| [`ins_update_velocity`](inc/damp/estimation/ins_eskf.hpp#L597) | Sparse velocity update: y = v at the IMU origin (no lever-arm) |
| [`ins_update_zupt`](inc/damp/estimation/ins_eskf.hpp#L616) | Zero-velocity update (ZUPT): velocity measurement of zero |
| [`kalman`](inc/damp/estimation/kalman.hpp#L103) | Steady-state Kalman filter design |
| [`luenberger`](inc/damp/estimation/luenberger.hpp#L138) | Design a Luenberger observer by robust pole placement (matrix form) (+3 more overloads) |
| [`margins_from_frf`](inc/damp/estimation/frequency_response.hpp#L771) | Gain / phase margins from a fixed on-target FRF table |
| [`mechanize_step_from_corrected`](inc/damp/estimation/ins_mechanization.hpp#L173) | One dead-reckoning step from an IMU sample (strapdown mechanization) |
| [`mhe`](inc/damp/estimation/mhe.hpp#L166) | Synthesize a constrained moving-horizon estimator |
| [`multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L133) | Build a multi-sine design payload from a configuration |
| [`ned_from_enu`](inc/damp/estimation/ins_mechanization.hpp#L106) | Map an ENU vector into NED (axis permute) |
| [`prbs`](inc/damp/estimation/excitation/prbs.hpp#L158) | Build a PRBS design payload from a configuration |
| [`ramp`](inc/damp/estimation/excitation/ramp.hpp#L94) | Build a ramp design payload from a configuration |
| [`reduced_luenberger`](inc/damp/estimation/luenberger.hpp#L271) | Design a reduced-order (Gopinath) observer by pole placement (matrix form) (+3 more overloads) |
| [`relay_autotune`](inc/damp/estimation/relay_autotune.hpp#L213) | Build a validated relay-autotune design payload |
| [`resonance_autotune_from_frf`](inc/damp/estimation/commissioning_frontends.hpp#L76) | Map FRF peaks and valleys to a fixed bank of biquad compensators |
| [`rls`](inc/damp/estimation/rls.hpp#L154) | Build scalar RLS design payload |
| [`rls_vector`](inc/damp/estimation/rls.hpp#L170) | Build vector RLS design payload |
| [`rows`](inc/damp/estimation/ins_eskf.hpp#L348) | Sparse y = F x for the 15-state INS first-order structure (G = I path) |
| [`schroeder_multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L180) | Build a multi-sine design payload with Schroeder phases applied |
| [`serial_arm_pose_design`](inc/damp/estimation/serial_arm_pose.hpp#L272) | Diagonal process / prior for joint-space pose filter |
| [`specific_force_at_rest`](inc/damp/estimation/ins_mechanization.hpp#L212) | Specific force [m/s²] a stationary IMU measures for the given orientation |
| [`step`](inc/damp/estimation/excitation/step.hpp#L91) | Build a single-step design payload |
| [`step_train`](inc/damp/estimation/excitation/step_train.hpp#L94) | Build a step-train design payload from a configuration |
| [`stepped_sine`](inc/damp/estimation/excitation/stepped_sine.hpp#L104) | Build a stepped-sine design payload from a configuration |
| [`synthetic_arm_imus`](inc/damp/estimation/serial_arm_pose.hpp#L879) | Build IMU samples for all frames (needs three configs); optional mounts |
| [`synthetic_link_imu`](inc/damp/estimation/serial_arm_pose.hpp#L813) | Synthetic IMU at a sense point (gyro + specific force) for SIL |
| [`with_attitude`](inc/damp/estimation/serial_arm_pose.hpp#L313) | Enable one-frame attitude ESKF on an existing joint design |
| [`with_schroeder_phases`](inc/damp/estimation/excitation/multi_sine.hpp#L160) | Assign Schroeder low-crest-factor phases to a multi-sine tone table |

## Filters & signal conditioning

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`Biquad`](inc/damp/filters/biquad.hpp#L34) | Second-order IIR (biquad) section runtime |
| [`BiquadCascade`](inc/damp/filters/biquad.hpp#L101) | Cascade of second-order sections (SOS) for higher-order IIR filters |
| [`Complementary`](inc/damp/filters/complementary.hpp#L32) | Scalar (1-D) complementary filter — fuse a fast rate with a slow absolute |
| [`Delay`](inc/damp/filters/delay.hpp#L27) | Discrete-time delay buffer |
| [`DsogiPll`](inc/damp/filters/pll.hpp#L282) | Dual-SOGI three-phase positive-sequence PLL (DSOGI-PLL) |
| [`FirDesignResult`](inc/damp/filters/fir.hpp#L67) | Window-method FIR design result |
| [`FirFilter`](inc/damp/filters/fir.hpp#L450) | Direct-form FIR filter runtime (tapped delay line + dot product) |
| [`FirstOrderCoeffs`](inc/damp/filters/iir_design.hpp#L30) | DSP coefficients for first-order IIR filter |
| [`FirType`](inc/damp/filters/fir.hpp#L38) | Ideal frequency-selective response for fir1 / fir_window |
| [`FirWindow`](inc/damp/filters/fir.hpp#L50) | Window applied to the ideal truncated impulse response |
| [`Goertzel`](inc/damp/filters/spectral.hpp#L46) | Generalized Goertzel single-bin DFT — amplitude/phase at one frequency |
| [`HarmonicAnalyzer`](inc/damp/filters/spectral.hpp#L173) | Harmonic analyzer — a Goertzel bank over a fundamental and K−1 harmonics |
| [`HighPass`](inc/damp/filters/highpass.hpp#L40) | First-order high-pass (washout) filter runtime |
| [`LowPass`](inc/damp/filters/lowpass.hpp#L49) | Nth-order low-pass filter |
| [`MedianFilter`](inc/damp/filters/median.hpp#L34) | Sliding-window median filter — nonlinear spike/outlier rejection |
| [`MovingAverage`](inc/damp/filters/moving_average.hpp#L62) | Moving-average (boxcar) filter — also a DC-preserving harmonic-notch comb |
| [`MSTOGI`](inc/damp/filters/sogi.hpp#L242) | Runtime MSTOGI with exact resonator and forward-Euler washout |
| [`NlmsFilter`](inc/damp/filters/fir.hpp#L546) | Normalized LMS adaptive FIR filter |
| [`RobustExactDifferentiator`](inc/damp/filters/differentiator.hpp#L63) | First-order robust exact differentiator (super-twisting differentiator) |
| [`SecondOrderCoeffs`](inc/damp/filters/iir_design.hpp#L65) | DSP coefficients for second-order IIR filter |
| [`SinglePhasePLL`](inc/damp/filters/pll.hpp#L43) | Single-Phase PLL |
| [`SOGI`](inc/damp/filters/sogi.hpp#L177) | Runtime SOGI wrapper around design::sogi(w0, alpha, Ts) |
| [`SogiFll`](inc/damp/filters/sogi.hpp#L330) | SOGI with a Frequency-Locked Loop — self-tuning single-tone tracker |
| [`ThreePhasePLL`](inc/damp/filters/pll.hpp#L141) | Synchronous-reference-frame (SRF) PLL for balanced three-phase input |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`bandpass`](inc/damp/filters/iir_design.hpp#L661) | Second-order band-pass filter (constant 0 dB peak gain) |
| [`butterworth_lowpass`](inc/damp/filters/iir_design.hpp#L314) | Butterworth low-pass filter design (+1 more overload) |
| [`comb_notch_window`](inc/damp/filters/moving_average.hpp#L34) | Window length for a moving-average comb that notches f_notch and all its harmonics: N = round(fs / f_notch) |
| [`continuous_lpf_exact_step`](inc/damp/filters/lowpass.hpp#L37) | Exact step of continuous first-order LPF ẏ = −ω_c (y − u) |
| [`fir1`](inc/damp/filters/fir.hpp#L218) | Window-method FIR design (normalized frequency, single cutoff) (+1 more overload) |
| [`fir1_hz`](inc/damp/filters/fir.hpp#L372) | Window-method FIR design from frequencies in Hz (+1 more overload) |
| [`fir_window`](inc/damp/filters/fir.hpp#L408) | Alias for fir1 (descriptive name) (+1 more overload) |
| [`highpass_1st`](inc/damp/filters/iir_design.hpp#L169) | First-order high-pass filter design (Tustin / bilinear) |
| [`highpass_2nd`](inc/damp/filters/iir_design.hpp#L684) | Second-order high-pass filter (RBJ) |
| [`highshelf`](inc/damp/filters/iir_design.hpp#L765) | High-shelf EQ filter: boost or cut everything above fc |
| [`lowpass_1st`](inc/damp/filters/iir_design.hpp#L136) | First-order low-pass filter design (+1 more overload) |
| [`lowpass_2nd`](inc/damp/filters/iir_design.hpp#L221) | Second-order low-pass filter design (+1 more overload) |
| [`lowpass_2nd_continuous`](inc/damp/filters/iir_design.hpp#L272) | Second-order low-pass filter design (continuous-time) |
| [`lowshelf`](inc/damp/filters/iir_design.hpp#L733) | Low-shelf EQ filter: boost or cut everything below fc |
| [`mstogi`](inc/damp/filters/sogi.hpp#L117) | Mixed Second/Third-Order Generalized Integrator (MSTOGI) (+1 more overload) |
| [`negative_sequence_ab`](inc/damp/filters/pll.hpp#L244) | Instantaneous negative-sequence αβ from a quadrature signal pair |
| [`notch`](inc/damp/filters/iir_design.hpp#L639) | Second-order band-reject (notch) filter |
| [`pade_delay_1st`](inc/damp/filters/iir_design.hpp#L378) | First-order Pade approximation of time delay (+1 more overload) |
| [`pade_delay_2nd`](inc/damp/filters/iir_design.hpp#L422) | Second-order Pade approximation of time delay (+1 more overload) |
| [`peaking`](inc/damp/filters/iir_design.hpp#L708) | Peaking (bell) EQ filter: boost or cut a band around f0 |
| [`positive_sequence_ab`](inc/damp/filters/pll.hpp#L223) | Instantaneous positive-sequence αβ from a quadrature signal pair |
| [`sogi`](inc/damp/filters/sogi.hpp#L39) | Second-Order Generalized Integrator (SOGI) design (+1 more overload) |
| [`to_coeffs`](inc/damp/filters/iir_design.hpp#L467) | Convert StateSpace system to first-order DSP coefficients (+3 more overloads) |

> **Pose / trajectory helpers (core)**

## Trajectory value types

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AxisMotionLimits`](inc/damp/trajectory/online_otg.hpp#L64) | Per-axis value / asymmetric accel–decel / jerk caps |
| [`AxisState`](inc/damp/trajectory/online_otg.hpp#L80) | Shared kinematic state for one axis: value and its derivative |
| [`CamFollower`](inc/damp/trajectory/cam.hpp#L481) | Runtime follower of an electronic cam table (design::CamProfile) |
| [`CamProfile`](inc/damp/trajectory/cam.hpp#L99) | Cam table: follower position as a function of leader parameter over one cycle |
| [`CartesianMove`](inc/damp/trajectory/cartesian_move.hpp#L114) | Path-preserving task-space move (Pipeline B / LIN) |
| [`DriveFeel`](inc/damp/trajectory/selector_governor.hpp#L69) | Stick / command feel knobs (deadband, expo, PIO pole frequency) |
| [`DriveSelectorGovernor`](inc/damp/trajectory/selector_governor.hpp#L87) | Two-axis drive governor: longitudinal v and yaw rate ω |
| [`ElectronicGear`](inc/damp/trajectory/cam.hpp#L155) | Electronic gearing law q = r·θ + φ (PLCopen MC_GearIn) |
| [`GearFollower`](inc/damp/trajectory/cam.hpp#L540) | Runtime electronic gear q = k_q(r θ + φ) + o_q with chain-rule rates |
| [`HolonomicSelectorGovernor`](inc/damp/trajectory/selector_governor.hpp#L189) | Three-axis holonomic governor: forward v_y, strafe v_x, yaw ω |
| [`InputShaper`](inc/damp/trajectory/input_shaper.hpp#L195) | Input-shaper runtime — convolves a command stream with the shaper impulses |
| [`InputShaperBank`](inc/damp/trajectory/input_shaper.hpp#L256) | Multi-axis input-shaper bank — one shaper per axis, shared buffer length |
| [`InputShaperResult`](inc/damp/trajectory/input_shaper.hpp#L70) | Input-shaper design result: impulse amplitudes and sample delays |
| [`JerkLimitedAxis`](inc/damp/trajectory/online_otg.hpp#L486) | One online OTG axis with optional override jerks (MaxOv stack array) |
| [`JerkRequest`](inc/damp/trajectory/online_otg.hpp#L89) | One requester's proposed jerk for this tick (no private integrator) |
| [`Jet`](inc/damp/trajectory/trajectory_types.hpp#L95) | Truncated jet of a scalar motion sample: d[k] = s⁽ᵏ⁾ |
| [`JointLimits`](inc/damp/trajectory/cartesian_move.hpp#L88) | Per-joint velocity and acceleration limits for a task-space move |
| [`LinearPath`](inc/damp/trajectory/cartesian_move.hpp#L62) | A straight-line path `p(s) = start + s·dir`, `s ∈ [0, length]` |
| [`PolynomialTrajectory`](inc/damp/trajectory/polynomial.hpp#L265) | Runtime evaluator for a precomputed polynomial trajectory |
| [`PolyTrajectory`](inc/damp/trajectory/polynomial.hpp#L73) | A synthesized polynomial trajectory: the coefficients of p(t) = Σ cᵢ·tⁱ over t ∈ [0, T], plus the duration |
| [`ScurveProfile`](inc/damp/trajectory/scurve.hpp#L72) | A synthesized jerk-limited (double-S) profile: a sequence of constant-jerk segments, evaluated exactly (cubic in t within a segment) |
| [`ScurveSegment`](inc/damp/trajectory/scurve.hpp#L56) | One constant-jerk segment, valid for t ∈ [t0, t0 + duration), with the position/velocity/acceleration cached at the segment start |
| [`ScurveTrajectory`](inc/damp/trajectory/scurve.hpp#L307) | Runtime evaluator for a precomputed jerk-limited (double-S) profile |
| [`ShaperType`](inc/damp/trajectory/input_shaper.hpp#L56) | Input-shaper family |
| [`SplineProfile`](inc/damp/trajectory/spline.hpp#L69) | A synthesized multi-waypoint spline: per-segment polynomial coefficients (ascending power, in segment-local time) plus the knot times |
| [`SplineTrajectory`](inc/damp/trajectory/spline.hpp#L391) | Runtime player for a multi-waypoint spline (design::SplineProfile) |
| [`ToppMove`](inc/damp/trajectory/topp.hpp#L123) | Time-optimal task-space move (path-preserving, pointwise minimum-time) |
| [`ToppProfile`](inc/damp/trajectory/topp.hpp#L66) | The scalar time-optimal path-timing produced by TOPP |
| [`TrajectoryBank`](inc/damp/trajectory/polynomial.hpp#L334) | Multi-axis coordination: time-scale each axis's profile to the slowest so a multi-DOF move starts and finishes synchronized ("linear" / coordinated joint moves — the feedforward reference for a manipulator) |
| [`TrajectoryBoundary`](inc/damp/trajectory/trajectory_types.hpp#L178) | Boundary conditions at one endpoint of a polynomial trajectory: a position and its time derivatives through jerk |
| [`TrajectoryLimits`](inc/damp/trajectory/trajectory_types.hpp#L45) | Asymmetric kinematic limits for a trapezoidal or S-curve motion profile |
| [`TrajectoryState`](inc/damp/trajectory/trajectory_types.hpp#L76) | A point on a motion profile: commanded position, velocity, acceleration |
| [`TrapezoidalProfile`](inc/damp/trajectory/trapezoidal.hpp#L56) | Planned trapezoidal profile: the segment durations, reached values, and boundary state needed to evaluate the trajectory at any time |
| [`TrapezoidalTrajectory`](inc/damp/trajectory/trapezoidal.hpp#L288) | Runtime evaluator for a precomputed trapezoidal profile |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`axis_snap_quantum`](inc/damp/trajectory/online_otg.hpp#L220) | Near-target snap quantum (~0.01% of scale; fixed-point-style residual floor) |
| [`cam_profile`](inc/damp/trajectory/cam.hpp#L136) | Build a CamProfile from a synthesized spline and cycle length |
| [`clamp_accel_asymmetric`](inc/damp/trajectory/online_otg.hpp#L101) | Asymmetric clamp of accel to +max_accel / −max_decel |
| [`clamp_jerk_accel_stop`](inc/damp/trajectory/online_otg.hpp#L179) | Clamp jerk so one step does not drive a past a_stop |
| [`clamp_rate_asymmetric`](inc/damp/trajectory/online_otg.hpp#L110) | Asymmetric clamp of accel (legacy name) |
| [`cubic_rise`](inc/damp/trajectory/cam.hpp#L237) | Cubic rest–rest rise (Hermite; C¹ ends — accel jumps into a dwell) |
| [`cubic_spline`](inc/damp/trajectory/spline.hpp#L239) | Cubic (C²) spline through points at times; clamped end velocities |
| [`cycloidal_rise`](inc/damp/trajectory/cam.hpp#L210) | Cycloidal rise of height h over unit leader fraction u ∈ [0, 1] |
| [`factorial`](inc/damp/trajectory/trajectory_types.hpp#L206) | k! (exact for modest k in float/double; used up through nonic BVP / jets) |
| [`falling_factorial`](inc/damp/trajectory/trajectory_types.hpp#L217) | Falling factorial i·(i−1)···(i−k+1) = i! / (i−k)! — the k-th derivative coefficient of tⁱ. Zero when i < k |
| [`in_jerk_terminal`](inc/damp/trajectory/online_otg.hpp#L153) | True when residual a is closing on rem and inside the jmax stop budget |
| [`input_shaper`](inc/damp/trajectory/input_shaper.hpp#L113) | Synthesize an input shaper for a second-order mode |
| [`integrate_jerk_axis`](inc/damp/trajectory/online_otg.hpp#L372) | Integrate jerk → accel → value (semi-implicit), land/snap on value_target (+1 more overload) |
| [`jerk_stop_budget`](inc/damp/trajectory/online_otg.hpp#L125) | Velocity change while ramping residual accel to zero at jmax |
| [`jerk_toward_accel`](inc/damp/trajectory/online_otg.hpp#L196) | Jerk that slews a → a_des in one step, sat to ±jmax, no reverse past a_des |
| [`make_cartesian_move`](inc/damp/trajectory/cartesian_move.hpp#L266) | Deduction-friendly factory (deduces NJoints from joint_limits, callables from args) |
| [`make_topp_move`](inc/damp/trajectory/topp.hpp#L361) | Deduction-friendly factory (deduces NJoints from joint_limits; pick NGrid) |
| [`max_accel_for_velocity_error`](inc/damp/trajectory/online_otg.hpp#L139) | Max \|a\| still unloadable under jmax for remaining velocity error |
| [`min_accel`](inc/damp/trajectory/polynomial.hpp#L213) | Minimum-acceleration (cubic) rest-to-rest move p0 → pT over duration T |
| [`min_jerk`](inc/damp/trajectory/polynomial.hpp#L207) | Minimum-jerk (quintic) rest-to-rest move p0 → pT over duration T (Flash–Hogan) |
| [`min_snap`](inc/damp/trajectory/polynomial.hpp#L219) | Minimum-snap (septic) rest-to-rest move p0 → pT over duration T (Mellinger–Kumar) |
| [`modified_sine_rise`](inc/damp/trajectory/cam.hpp#L327) | Modified-sine rise (Norton): low peak accel, common industrial law |
| [`modified_trapezoid_rise`](inc/damp/trajectory/cam.hpp#L368) | Modified-trapezoid rise (constant-accel flanks + harmonic transitions) |
| [`nest_a_des`](inc/damp/trajectory/online_otg.hpp#L249) | Nest a_des under hard a limits and the jmax unload budget for rem |
| [`operator_a_des`](inc/damp/trajectory/online_otg.hpp#L300) | PIO a_des (legacy name used by drive governor docs) |
| [`otg_jerk`](inc/damp/trajectory/online_otg.hpp#L315) | OTG jerk from any a_des generator, with terminal unload toward value_target |
| [`periodic_cubic_spline`](inc/damp/trajectory/spline.hpp#L363) | Periodic cubic (C² across the cycle wrap) |
| [`periodic_quintic_spline`](inc/damp/trajectory/spline.hpp#L371) | Periodic quintic (C⁴ across the cycle wrap) |
| [`periodic_spline`](inc/damp/trajectory/spline.hpp#L279) | Periodic multi-waypoint spline (cyclic wrap of derivatives 1…Order−1) |
| [`pio_a_des`](inc/damp/trajectory/online_otg.hpp#L269) | PIO a_des generator: a_des = α rem + nest + rest min |
| [`poly_rise`](inc/damp/trajectory/polynomial.hpp#L242) | Rest–rest polynomial rise of height h over unit parameter u ∈ [0,1] |
| [`poly_trajectory`](inc/damp/trajectory/polynomial.hpp#L149) | Synthesize a fixed-duration polynomial matching endpoint derivative jets (+1 more overload) |
| [`predicted_abs_value`](inc/damp/trajectory/online_otg.hpp#L440) | Predict \|value\| after one step of jerk j (selector ranking) |
| [`quintic_rise`](inc/damp/trajectory/cam.hpp#L266) | Quintic rest–rest rise (degree 5; C² ends) |
| [`quintic_spline`](inc/damp/trajectory/spline.hpp#L251) | Quintic (C⁴ — jerk- and snap-continuous) spline; clamped end velocity + accel |
| [`sample_rise_law`](inc/damp/trajectory/cam.hpp#L445) | Sample a unit-parameter rise law onto NPts knots on [0, rise_length] |
| [`scurve`](inc/damp/trajectory/scurve.hpp#L196) | Synthesize a minimum-time jerk-limited (7-segment double-S) profile from (Xi, Vi) to (Xf, Vf) under asymmetric kinematic limits (+1 more overload) |
| [`select_jerk`](inc/damp/trajectory/online_otg.hpp#L452) | Pick the jerk that most reduces \|value\| if any override beats the operator |
| [`septic_rise`](inc/damp/trajectory/cam.hpp#L299) | Septic rest–rest rise (degree 7; zero end jerk; C³ into dwells) |
| [`shape_stick`](inc/damp/trajectory/selector_governor.hpp#L56) | Joystick / RC stick feel: scaled center dead zone then RC expo |
| [`sliding_jerk_toward_value`](inc/damp/trajectory/online_otg.hpp#L347) | PIO operator path: pio_a_des then otg_jerk |
| [`spline`](inc/damp/trajectory/spline.hpp#L155) | Synthesize a multi-waypoint spline through points at times |
| [`terminal_unload_jerk`](inc/damp/trajectory/online_otg.hpp#L167) | Continuous land jerk toward (v , a → 0): j = −a²/(2 rem) |
| [`trapezoidal`](inc/damp/trajectory/trapezoidal.hpp#L206) | Synthesize the minimum-time asymmetric trapezoidal profile from (Xi, Vi) to (Xf, Vf) under the given limits (+1 more overload) |
| [`truncate_residual`](inc/damp/trajectory/online_otg.hpp#L231) | Kill residuals smaller than q (float dust → 0); leave larger values alone |
| [`unload_accel_jerk`](inc/damp/trajectory/online_otg.hpp#L207) | Pure unload of residual a toward 0 under jmax (no reverse through 0) |
| [`wrap_leader`](inc/damp/trajectory/cam.hpp#L72) | Wrap leader position into [0, L) (true modulo; supports reverse) |

## Kinematics / pose

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`ArmIkResult`](inc/damp/kinematics/serial_arm.hpp#L164) | Result of a numerical inverse-kinematics solve |
| [`CartesianMap`](inc/damp/kinematics/motion_maps.hpp#L44) | Cartesian gantry: independent per-axis affine map `task = scale·act + offset` (the "kinematics" is the identity, exposed for a uniform forward/inverse interface) |
| [`CoreXY`](inc/damp/kinematics/motion_maps.hpp#L74) | CoreXY belt mapping (2 motors A/B → Cartesian X/Y) |
| [`DeltaForward`](inc/damp/kinematics/motion_maps.hpp#L129) | Result of a delta forward solve: the end-effector pose (orientation fixed to identity — deltas are 3-DOF translational) + validity |
| [`DeltaInverse`](inc/damp/kinematics/motion_maps.hpp#L121) | Result of a delta inverse solve: the three actuator values + reachability |
| [`DhChain`](inc/damp/kinematics/serial_arm.hpp#L108) | An N-joint DH chain (the arm geometry) |
| [`DhJoint`](inc/damp/kinematics/serial_arm.hpp#L90) | One joint's standard (distal) DH parameters and motion limits |
| [`FiveBar`](inc/damp/kinematics/scara.hpp#L92) | Planar five-bar parallel manipulator (parallel SCARA) |
| [`FiveBarForward`](inc/damp/kinematics/scara.hpp#L76) | Forward-kinematics result: the end-effector point (x, y) + validity |
| [`FiveBarGeometry`](inc/damp/kinematics/scara.hpp#L57) | Symmetric five-bar geometry (two base motors, equal proximal/distal links) |
| [`FiveBarInverse`](inc/damp/kinematics/scara.hpp#L69) | Inverse-kinematics result: the two motor angles [rad] + reachability |
| [`JointType`](inc/damp/kinematics/serial_arm.hpp#L74) | Joint actuation type for a DH joint |
| [`LinearDelta`](inc/damp/kinematics/motion_maps.hpp#L285) | Linear delta robot — per-carriage closed-form inverse, sphere- trilateration forward. Towers at 90°, 210°, 330° |
| [`LinearDeltaGeometry`](inc/damp/kinematics/motion_maps.hpp#L269) | Linear delta geometry (three vertical carriages, fixed-length rods) |
| [`PolarMap`](inc/damp/kinematics/motion_maps.hpp#L95) | Polar / R-θ mapping (radius + angle ↔ Cartesian X/Y) |
| [`Pose`](inc/damp/kinematics/pose.hpp#L60) | Rigid-body pose: a translation and an orientation (unit quaternion) |
| [`RotaryDelta`](inc/damp/kinematics/motion_maps.hpp#L161) | Rotary delta robot — closed-form inverse, quadratic-intersection forward |
| [`RotaryDeltaGeometry`](inc/damp/kinematics/motion_maps.hpp#L140) | Rotary delta geometry (three base servos, parallelogram arms) |
| [`SerialArm`](inc/damp/kinematics/serial_arm.hpp#L209) | Serial N-DOF manipulator runtime (joint-space N free; task space SE(3)) |
| [`SerialArmConfig`](inc/damp/kinematics/serial_arm.hpp#L176) | Validated serial-arm configuration (the design payload) |
| [`StewartConfig`](inc/damp/kinematics/stewart.hpp#L113) | Validated Stewart configuration (the design payload) |
| [`StewartForward`](inc/damp/kinematics/stewart.hpp#L102) | Result of a forward (Newton–Raphson) solve |
| [`StewartGeometry`](inc/damp/kinematics/stewart.hpp#L65) | Rig geometry: the six fixed base anchors `bᵢ`, the six moving-platform anchors `pᵢ`, the actuator stroke limits, and the nominal home height |
| [`StewartInverse`](inc/damp/kinematics/stewart.hpp#L95) | Result of an inverse solve: the six leg lengths + a stroke-window flag |
| [`StewartPlatform`](inc/damp/kinematics/stewart.hpp#L165) | Gough–Stewart platform runtime — closed-form inverse, Newton forward |
| [`Translation3`](inc/damp/kinematics/pose.hpp#L37) | A 3-D translation — a thin Vec3 with domain-named conveniences |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`arm_spherical_wrist`](inc/damp/kinematics/serial_arm.hpp#L601) | Tier-2 builder for a standard 6R elbow arm with a spherical wrist |
| [`dh_joint_axis_in_parent`](inc/damp/kinematics/serial_arm.hpp#L505) | Standard-DH joint axis in the parent frame: $`z_{j-1} = (0,0,1)`$ |
| [`five_bar_symmetric`](inc/damp/kinematics/scara.hpp#L241) | Build a symmetric five-bar parallel SCARA |
| [`frame_orientations`](inc/damp/kinematics/serial_arm.hpp#L546) | Orientations of DH frames {0,…,N} at configuration q |
| [`is_spherical_wrist`](inc/damp/kinematics/serial_arm.hpp#L564) | Spherical-wrist (Pieper) criterion for a 6R chain: axes 4-5-6 intersect, i.e. `a₄ = a₅ = a₆ = 0` and `d₅ = 0` (within eps). `a₄` (joints[3].a) is required so the wrist centre sits at the joint-4/5/6 concurrency point |
| [`joint_angles_from_frame_orientations`](inc/damp/kinematics/serial_arm.hpp#L525) | Recover revolute joint angles from world orientations of DH frames {0,…,N} |
| [`revolute_angle_about_axis`](inc/damp/kinematics/serial_arm.hpp#L491) | Signed revolute angle about a known axis from two body→world orientations |
| [`scara_arm`](inc/damp/kinematics/scara.hpp#L262) | Build a series SCARA (RRPR) as a 4-joint DH chain |
| [`select_nearest`](inc/damp/kinematics/serial_arm.hpp#L455) | Pick the solution branch nearest a reference configuration |
| [`serial_arm`](inc/damp/kinematics/serial_arm.hpp#L578) | Validate a serial-arm DH chain and flag a spherical wrist |
| [`stewart`](inc/damp/kinematics/stewart.hpp#L317) | Validate a hand-entered Stewart geometry and confirm the home pose is reachable |
| [`stewart_symmetric`](inc/damp/kinematics/stewart.hpp#L352) | Tier-2 builder for the common symmetric hexagonal layout |

## Motor control pack (if present)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`Acim`](inc/damp/motor/acim.hpp#L282) | AC induction machine (ACIM) electrical stage — indirect FOC |
| [`AcimConfig`](inc/damp/motor/acim.hpp#L254) | Configuration for an Acim electrical stage |
| [`AcimServoConfig`](inc/damp/motor/acim.hpp#L440) | Flat configuration for an ACIM cascade servo (make_acim_servo) |
| [`AlphaBeta`](inc/damp/transforms.hpp#L152) | Alpha-beta (stationary-frame) component pair |
| [`AlphaBetaZero`](inc/damp/transforms.hpp#L218) | Alpha-beta-zero (stationary-frame) component triple |
| [`AngleGenerator`](inc/damp/motor/scalar_control.hpp#L54) | Fixed-frequency electrical-angle generator with a slew-limited frequency |
| [`BrushedDc`](inc/damp/motor/brushed_dc.hpp#L66) | Brushed-DC electrical stage: single armature current loop, one H-bridge duty |
| [`BrushedDcConfig`](inc/damp/motor/brushed_dc.hpp#L44) | Configuration for a BrushedDc electrical stage |
| [`BrushedDcServoConfig`](inc/damp/motor/brushed_dc.hpp#L144) | Flat configuration for a brushed-DC cascade servo (make_brushed_servo) |
| [`Cascade`](inc/damp/motor/cascade.hpp#L94) | Position→velocity→torque cascade: the default Servo outer controller |
| [`CascadeBandwidths`](inc/damp/motor/cascade.hpp#L57) | The three bandwidth knobs of the position/velocity/current cascade |
| [`CascadeConfig`](inc/damp/motor/cascade.hpp#L73) | Configuration for a Cascade controller (the mechanical / outer-loop half) |
| [`ControlMode`](inc/damp/motor/cascade.hpp#L42) | Which cascade outer loops are active |
| [`Convention`](inc/damp/transforms.hpp#L69) | Scaling convention for the Clarke/Park family |
| [`DcBusLimiter`](inc/damp/motor/limits.hpp#L66) | Holds the inverter's torque current within DC-bus current/power limits |
| [`DcBusLimits`](inc/damp/motor/limits.hpp#L31) | DC-bus current and voltage limits for an inverter |
| [`DcBusState`](inc/damp/motor/limits.hpp#L42) | DC-bus state and the torque-current derate it implies |
| [`DeadTimeBlanking`](inc/damp/motor/dead_time.hpp#L78) | Blanking geometry assumed by the average error model |
| [`DeadTimeCompensator`](inc/damp/motor/dead_time.hpp#L279) | Runtime dead-time compensator (holds blanking parameters) |
| [`DirectQuadrature`](inc/damp/transforms.hpp#L84) | Direct-quadrature (rotor-frame) component pair |
| [`DirectQuadratureZero`](inc/damp/transforms.hpp#L238) | Direct-quadrature-zero (rotor-frame) component triple |
| [`DqCommand`](inc/damp/motor/foc.hpp#L331) | Result of FOController::current_controller(): the dq voltage command plus its saturation signals |
| [`Drive`](inc/damp/motor/drive.hpp#L88) | The motor-agnostic drive core: encoder tracker + a Motor stage + torque limiting |
| [`DriveNameplateLimits`](inc/damp/motor/drive_design.hpp#L41) | Nameplate-derived drive limits at one current operating point for base speed |
| [`DriveResult`](inc/damp/motor/drive.hpp#L46) | Output of one electrical step: the switch duties plus realized torque and flags |
| [`EemfObserver`](inc/damp/motor/sensorless_eemf.hpp#L93) | Extended-EMF sensorless observer (salient mid-speed) |
| [`EemfResult`](inc/damp/motor/sensorless_eemf.hpp#L32) | EEMF design payload |
| [`EncoderTracker`](inc/damp/motor/encoder_tracker.hpp#L61) | Kinematic angle/speed tracking observer (PLL) over a wrapped absolute position |
| [`EncoderTrackerConfig`](inc/damp/motor/encoder_tracker.hpp#L23) | Configuration for EncoderTracker |
| [`FetLossModel`](inc/damp/motor/thermal.hpp#L161) | First-order inverter FET loss model (conduction + switching) |
| [`FieldExcitation`](inc/damp/motor/sepex.hpp#L240) | How the field winding is powered |
| [`FieldWeakening`](inc/damp/motor/field_weakening.hpp#L141) | Field-weakening current-reference regulator (voltage-feedback or feedforward) |
| [`FieldWeakeningConfig`](inc/damp/motor/field_weakening.hpp#L89) | Configuration for FieldWeakening |
| [`FieldWeakeningPolicy`](inc/damp/motor/field_weakening.hpp#L205) | Concept for a pluggable field-weakening / current-reference policy |
| [`FluxSensorlessAdapter`](inc/damp/motor/sensorless.hpp#L104) | Adapter: SensorlessEstimator → SensorlessObserver |
| [`FocCascadeExperimentOutput`](inc/damp/motor/foc_cascade_experiment.hpp#L70) | Per-tick output of a FOC cascade velocity experiment |
| [`FocCascadeExperimentStatus`](inc/damp/motor/foc_cascade_experiment.hpp#L56) | Lifecycle of a FOC cascade velocity experiment |
| [`FocCascadeTuneResult`](inc/damp/motor/foc_autotune.hpp#L78) | FOC cascade design result (current + velocity + position) |
| [`FocElectricalParams`](inc/damp/motor/foc_cascade_experiment.hpp#L89) | Electrical nameplate inputs shared by velocity experiment drivers |
| [`FOController`](inc/damp/motor/foc.hpp#L346) | Field-oriented current controller (dq PI + decoupling + SVPWM duties) |
| [`FocResult`](inc/damp/motor/foc.hpp#L312) | Result of one FOController::step(), carrying the actuator command plus the saturation/measurement signals an outer (velocity/position) loop needs to propagate anti-windup back up a cascade |
| [`FocVelocityAmigoConfig`](inc/damp/motor/foc_cascade_experiment.hpp#L179) | Configuration for FocVelocityAmigoExperiment |
| [`FocVelocityAmigoExperiment`](inc/damp/motor/foc_cascade_experiment.hpp#L275) | Runtime step + relay experiment → FOC cascade velocity AMIGO gains |
| [`FocVelocityAmigoResult`](inc/damp/motor/foc_cascade_experiment.hpp#L221) | Design payload for the velocity-AMIGO experiment (validated config) |
| [`FocVelocityFrfConfig`](inc/damp/motor/foc_cascade_experiment.hpp#L506) | Configuration for FocVelocityFrfExperiment |
| [`FocVelocityFrfExperiment`](inc/damp/motor/foc_cascade_experiment.hpp#L578) | Runtime stepped-sine FRF experiment → FOC cascade velocity PI gains |
| [`FocVelocityFrfResult`](inc/damp/motor/foc_cascade_experiment.hpp#L541) | Validated design payload for the velocity FRF experiment |
| [`FwMethod`](inc/damp/motor/field_weakening.hpp#L79) | Field-weakening law selection |
| [`HallDecoder`](inc/damp/motor/commutation.hpp#L69) | 120° Hall-effect decoder: sector, coarse angle, direction, and speed |
| [`HallState`](inc/damp/motor/commutation.hpp#L40) | One three-phase Hall reading decoded to an electrical sector and its kinematics |
| [`HfInjector`](inc/damp/motor/sensorless_hfi.hpp#L142) | High-frequency voltage injector (estimated dq) |
| [`HfiObserver`](inc/damp/motor/sensorless_hfi.hpp#L193) | Pulsating-HFI saliency tracker (pre-polarity: locks θ or θ+π) |
| [`HfiPolarityIdent`](inc/damp/motor/sensorless_hfi.hpp#L318) | Polarity ID via short +d voltage pulse after saliency lock |
| [`HfiResult`](inc/damp/motor/sensorless_hfi.hpp#L39) | Pulsating-HFI design payload |
| [`HybridSensorlessFrontEnd`](inc/damp/motor/sensorless_hfi.hpp#L399) | Hybrid HFI + mid-speed front-end |
| [`IfCommand`](inc/damp/motor/scalar_control.hpp#L156) | Result of one IfController step: an open-loop angle and a dq current reference |
| [`IfController`](inc/damp/motor/scalar_control.hpp#L182) | Forced-current (I-F) open-loop startup controller |
| [`InstantaneousPower`](inc/damp/transforms.hpp#L566) | Instantaneous active and reactive power |
| [`JunctionEstimator`](inc/damp/motor/thermal.hpp#L220) | FET junction-temperature estimator: case temperature plus a thermal model |
| [`MechanicalEstimator`](inc/damp/motor/mechanical_estimator.hpp#L128) | Cheap-predict mechanical estimator for position, speed, and load torque |
| [`MechanicalEstimatorConfig`](inc/damp/motor/mechanical_estimator.hpp#L92) | Configuration for MechanicalEstimator |
| [`Modulator`](inc/damp/motor/modulation.hpp#L496) | Thin scheme-holding wrapper for FOC / deploy paths |
| [`Motor`](inc/damp/motor/drive.hpp#L64) | The pluggable electrical stage of a drive for one machine type |
| [`MtpaReference`](inc/damp/motor/mtpa.hpp#L127) | Maximum-torque-per-ampere current-reference generator (PMSM / IPMSM / SynRM) |
| [`NameplateDesignResult`](inc/damp/motor/acim.hpp#L527) | Nameplate → cascade tune + filled motor::AcimServoConfig (Design Is Deploy entry) (+2 more overloads) |
| [`NameplateLimits`](inc/damp/motor/acim.hpp#L192) | ACIM nameplate limits (base speed, peak torque, rated power) (+2 more overloads) |
| [`NoFieldWeakening`](inc/damp/motor/field_weakening.hpp#L220) | Null field-weakening policy — passes the base reference through unchanged |
| [`PhaseCalibrationCommand`](inc/damp/motor/calibration.hpp#L54) | One step's output from PhaseParameterCalibrator |
| [`PhaseCalibrationConfig`](inc/damp/motor/calibration.hpp#L34) | Configuration for online phase resistance/inductance commissioning |
| [`PhaseDrive`](inc/damp/motor/commutation.hpp#L143) | Per-phase drive state for one half-bridge under six-step commutation |
| [`PhaseParameterCalibrator`](inc/damp/motor/calibration.hpp#L111) | Online phase R/L identification by recursive least squares (PRBS injected) |
| [`Pmsm`](inc/damp/motor/pmsm.hpp#L71) | Field-oriented PMSM/BLDC electrical stage |
| [`PmsmConfig`](inc/damp/motor/pmsm.hpp#L45) | Configuration for a Pmsm electrical stage |
| [`PmsmNameplateDesignResult`](inc/damp/motor/drive_design.hpp#L99) | Full nameplate design result: limits + cascade tune + deployable config |
| [`PmsmServoConfig`](inc/damp/motor/servo.hpp#L105) | Flat configuration for a PMSM cascade servo (make_pmsm_servo) |
| [`PwmScheme`](inc/damp/motor/modulation.hpp#L82) | Carrier-based three-phase VSI modulation scheme |
| [`ReferredRL`](inc/damp/motor/sepex.hpp#L164) | Field plant as seen from the rotary-transformer primary (referred) |
| [`ResistiveLossModel`](inc/damp/motor/thermal.hpp#L196) | Minimal conduction-only loss model for a weak datasheet |
| [`RotaryTransformerField`](inc/damp/motor/sepex.hpp#L250) | Rotary-transformer field path parameters (primary-referred control) |
| [`SensorlessEstimator`](inc/damp/motor/sensorless_flux.hpp#L82) | Sensorless rotor flux/position estimator for a PMSM, with optional sensor fusion |
| [`SensorlessObserver`](inc/damp/motor/sensorless.hpp#L62) | Sensorless observer concept (electrical frame) |
| [`Sepex`](inc/damp/motor/sepex.hpp#L290) | Separately-excited DC electrical stage (armature + field bridges) |
| [`SepexConfig`](inc/damp/motor/sepex.hpp#L262) | Configuration for a Sepex electrical stage |
| [`SepexServoConfig`](inc/damp/motor/sepex.hpp#L468) | Flat configuration for a SEPEX cascade servo (make_sepex_servo) |
| [`SequenceComponents`](inc/damp/transforms.hpp#L623) | Symmetrical (sequence) components of a three-phase phasor set |
| [`Servo`](inc/damp/motor/servo.hpp#L52) | A drive plus an outer torque law — convenience for the cascade daily path |
| [`SinCosCalibration`](inc/damp/motor/resolver.hpp#L48) | Per-channel gain/offset trims for a SinCosDecoder (hardware calibration) |
| [`SinCosDecoder`](inc/damp/motor/resolver.hpp#L78) | Resolver / sin-cos-encoder decoder: (sin, cos) → angle, speed, multi-turn position |
| [`SinCosState`](inc/damp/motor/resolver.hpp#L57) | Decoded state from one SinCosDecoder step |
| [`SixStepCommutator`](inc/damp/motor/commutation.hpp#L171) | Six-step (trapezoidal) commutator: electrical sector → phase energization |
| [`SixStepDrive`](inc/damp/motor/commutation.hpp#L150) | One commutation step: the drive state of each of the three phases {a, b, c} |
| [`SmoObserver`](inc/damp/motor/sensorless_smo.hpp#L97) | αβ sliding-mode back-EMF sensorless observer |
| [`SmoResult`](inc/damp/motor/sensorless_smo.hpp#L33) | SMO design payload |
| [`Stepper`](inc/damp/motor/stepper.hpp#L71) | Two-phase field-oriented stepper electrical stage |
| [`StepperConfig`](inc/damp/motor/stepper.hpp#L49) | Configuration for a Stepper electrical stage |
| [`StepperServoConfig`](inc/damp/motor/stepper.hpp#L150) | Flat configuration for a stepper cascade servo (make_stepper_servo) |
| [`SvmDuties`](inc/damp/motor/modulation.hpp#L106) | Result of a duty-map: half-bridge duties plus an over-modulation flag |
| [`ThermalLimiter`](inc/damp/motor/thermal.hpp#L308) | Derates the current command from a temperature (Tj for FETs, winding for the motor) |
| [`ThermalLimits`](inc/damp/motor/thermal.hpp#L272) | A derating curve plus a hard fault threshold |
| [`ThermalLossModel`](inc/damp/motor/thermal.hpp#L142) | A loss model usable by JunctionEstimator |
| [`ThermalState`](inc/damp/motor/thermal.hpp#L292) | State from a ThermalLimiter evaluation |
| [`TorqueSpeedPoint`](inc/damp/motor/spm.hpp#L94) | One sample of the SPM torque–speed envelope (MTPA i_d=0 branch) |
| [`TwoMassEstimator`](inc/damp/motor/two_mass_estimator.hpp#L209) | Dual-encoder two-mass mechanical estimator (shaft compliance) |
| [`TwoMassEstimatorConfig`](inc/damp/motor/two_mass_estimator.hpp#L166) | Configuration for TwoMassEstimator |
| [`VfController`](inc/damp/motor/scalar_control.hpp#L113) | Constant Volts-per-Hertz open-loop scalar controller |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`angle_delay_compensation`](inc/damp/motor/scalar_control.hpp#L224) | Compensate an electrical angle for the control/PWM transport delay |
| [`apply_voltage_current_limits`](inc/damp/motor/ipm.hpp#L94) | Apply analytic voltage FW then current-circle clamp (open-loop, no integrator) |
| [`back_emf`](inc/damp/motor/sepex.hpp#L81) | Armature back-EMF e_a = K_af i_f ω [V] (omega in rad/s) |
| [`bake_torque_limit_lut`](inc/damp/motor/ipm.hpp#L235) | Bake a 1-D torque ceiling vs mechanical speed into a Lut1D (float-friendly) |
| [`base_speed`](inc/damp/motor/foc.hpp#L294) | Base (corner) electrical speed where the voltage circle is first hit (+1 more overload) |
| [`base_speed_elec`](inc/damp/motor/acim.hpp#L179) | Electrical base speed [elec rad/s] at rated flux and voltage ceiling |
| [`cauer_thermal_ss`](inc/damp/motor/thermal.hpp#L109) | Continuous state-space model of a physical Cauer RC thermal ladder |
| [`clamp_current_circle`](inc/damp/motor/ipm.hpp#L78) | Clamp a dq reference to the current circle of radius Imax (preserve angle) |
| [`clarke_park_transform`](inc/damp/transforms.hpp#L444) | Fused Clarke-Park transform (abc → dq) |
| [`clarke_park_zero_transform`](inc/damp/transforms.hpp#L527) | Fused Clarke-Park transform with zero (abc → dq0) |
| [`clarke_transform`](inc/damp/transforms.hpp#L343) | Clarke transform (abc → αβ) |
| [`clarke_zero_transform`](inc/damp/transforms.hpp#L273) | Zero-retaining Clarke transform (abc → αβ0) |
| [`compensate_duties`](inc/damp/motor/dead_time.hpp#L227) | Apply duty compensation to three half-bridge duties (clamped to [0, 1]) |
| [`compensate_phase_voltages`](inc/damp/motor/dead_time.hpp#L247) | Apply phase-voltage compensation before modulation |
| [`constant_power_torque_limit`](inc/damp/motor/foc.hpp#L80) | Constant-power torque ceiling Tₘₐₓ = P_rated / \|ω\| |
| [`current_loop_pi`](inc/damp/motor/foc.hpp#L63) | Current-loop PI gains by closed-loop pole placement on the R–L plant |
| [`current_magnitude`](inc/damp/motor/ipm.hpp#L70) | Stator current magnitude on the current circle |
| [`dead_time_current_sign`](inc/damp/motor/dead_time.hpp#L111) | Soft sign used by complementary blanking: ∈ {−1, 0, +1} |
| [`dead_time_duty_compensation`](inc/damp/motor/dead_time.hpp#L182) | Duty correction δd to add to a commanded half-bridge duty |
| [`dead_time_phase_voltage_compensation`](inc/damp/motor/dead_time.hpp#L206) | Phase-voltage correction [V] dual of dead_time_duty_compensation |
| [`dead_time_pole_voltage_error`](inc/damp/motor/dead_time.hpp#L145) | Average pole-voltage error Δv = v_avg − d V_dc |
| [`dead_time_ratio`](inc/damp/motor/dead_time.hpp#L97) | Dead-time ratio r = T_d / T_s (fraction of the PWM period) |
| [`derate_window`](inc/damp/motor/thermal.hpp#L40) | A two-breakpoint derating curve: 1 below derate_start, 0 at cutoff |
| [`dpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L268) | DPWM zero-sequence for a chosen discontinuous scheme |
| [`drive_nameplate_limits`](inc/damp/motor/drive_design.hpp#L66) | Build voltage/current/base-speed limits from bus, current ceiling, and machine |
| [`duties_from_phase_voltages`](inc/damp/motor/modulation.hpp#L360) | Map phase voltages + zero-sequence to clamped half-bridge duties |
| [`electromagnetic_torque`](inc/damp/motor/foc.hpp#L248) | Electromagnetic torque produced by a dq current (salient PMSM) |
| [`field_for_speed`](inc/damp/motor/sepex.hpp#L134) | Maximum field current that keeps e_a ≤ V_{a,max} at speed |
| [`field_weakening_id`](inc/damp/motor/field_weakening.hpp#L58) | Feedforward field-weakening d-axis current from the voltage ellipse |
| [`flux_for_speed`](inc/damp/motor/acim.hpp#L163) | Rotor-flux command that keeps referred flux inside the voltage circle |
| [`flux_from_id`](inc/damp/motor/acim.hpp#L99) | Steady-state rotor flux from magnetizing current: λ_r = Lₘ i_d [Wb] |
| [`flux_from_Kv`](inc/damp/motor/foc.hpp#L184) | PM flux linkage from the datasheet velocity constant Kᵥ |
| [`flux_from_torque_constant`](inc/damp/motor/foc.hpp#L127) | PM flux linkage from a motor's torque constant (amplitude-invariant) |
| [`foc_autotune`](inc/damp/motor/foc_autotune.hpp#L225) | Short alias for the nameplate FOC cascade designer |
| [`foc_cascade_tune`](inc/damp/motor/foc_autotune.hpp#L132) | Nameplate / post-calibration FOC cascade tune (no experiment) |
| [`foc_cascade_tune_velocity_amigo`](inc/damp/motor/foc_autotune.hpp#L194) | Velocity-loop PI from a biased-relay ultimate point (current loop already closed) |
| [`foc_velocity_amigo`](inc/damp/motor/foc_cascade_experiment.hpp#L236) | Build a validated velocity-AMIGO experiment payload |
| [`foc_velocity_frf`](inc/damp/motor/foc_cascade_experiment.hpp#L558) | Build a validated velocity-FRF experiment payload |
| [`foster_thermal_ss`](inc/damp/motor/thermal.hpp#L72) | Continuous state-space model of a Foster RC thermal network |
| [`hfi`](inc/damp/motor/sensorless_hfi.hpp#L79) | Design pulsating HFI from nameplate / PWM headroom |
| [`ia_from_torque`](inc/damp/motor/sepex.hpp#L101) | Armature current for a torque command at a field current |
| [`instantaneous_power`](inc/damp/transforms.hpp#L605) | Instantaneous active and reactive power from dq quantities (+1 more overload) |
| [`inverse_clarke_transform`](inc/damp/transforms.hpp#L360) | Inverse Clarke transform (αβ → abc) |
| [`inverse_clarke_zero_transform`](inc/damp/transforms.hpp#L303) | Inverse zero-retaining Clarke transform (αβ0 → abc) |
| [`inverse_park_clarke_transform`](inc/damp/transforms.hpp#L473) | Fused inverse Park-Clarke transform (dq → abc) |
| [`inverse_park_clarke_zero_transform`](inc/damp/transforms.hpp#L540) | Fused inverse Park-Clarke transform with zero (dq0 → abc) |
| [`inverse_park_transform`](inc/damp/transforms.hpp#L417) | Inverse Park transform (dq → αβ) |
| [`inverse_park_zero_transform`](inc/damp/transforms.hpp#L508) | Inverse Park transform with zero passthrough (dq0 → αβ0) |
| [`inverse_symmetrical_components`](inc/damp/transforms.hpp#L684) | Inverse symmetrical-component transform (012 → abc) |
| [`iq_from_torque`](inc/damp/motor/acim.hpp#L131) | q-axis current for a torque command at a given rotor flux (+1 more overload) |
| [`is_salient`](inc/damp/motor/ipm.hpp#L62) | True when \|Lq − Ld\| is meaningful saliency (same test as MTPA) |
| [`leakage_factor`](inc/damp/motor/acim.hpp#L72) | Leakage factor σ = 1 − Lₘ²/(L_s L_r) |
| [`linear_modulation_voltage`](inc/damp/motor/modulation.hpp#L130) | Peak phase voltage at the linear SVPWM hexagon limit |
| [`make_acim_servo`](inc/damp/motor/acim.hpp#L491) | Build the common `Servo<Cascade\<T\>, Acim\<T\>>` from one flat config |
| [`make_brushed_servo`](inc/damp/motor/brushed_dc.hpp#L164) | Build the common `Servo<Cascade\<T\>, BrushedDc\<T\>>` from one flat config |
| [`make_ipm_servo`](inc/damp/motor/servo.hpp#L214) | Build `Servo<Cascade, Pmsm<T, FieldWeakening>>` for the IPM daily path |
| [`make_pmsm_servo`](inc/damp/motor/servo.hpp#L190) | Build the common `Servo<Cascade\<T\>, Pmsm<T, Policy>>` from one flat config |
| [`make_stepper_servo`](inc/damp/motor/stepper.hpp#L170) | Build the common `Servo<Cascade\<T\>, Stepper\<T\>>` from one flat config |
| [`max_iq`](inc/damp/motor/spm.hpp#L61) | Max positive i_q for an SPM (i_d = 0) at electrical speed |
| [`max_torque_at_speed`](inc/damp/motor/ipm.hpp#L124) | Max \|Te\| at a mechanical speed under V and I limits (MTPA + open-loop FW) (+1 more overload) |
| [`mech_turns_from_elec`](inc/damp/motor/sensorless.hpp#L78) | Mechanical single-turn angle [turns] from electrical angle [rad] |
| [`mech_turns_per_s_from_elec`](inc/damp/motor/sensorless.hpp#L88) | Mechanical speed [turns/s] from electrical velocity [elec rad/s] |
| [`modulation_duty_cycles`](inc/damp/motor/modulation.hpp#L399) | Carrier-based VSI duty cycles from an αβ voltage command |
| [`modulation_index`](inc/damp/motor/modulation.hpp#L168) | Modulation index relative to the linear SVPWM circle |
| [`modulation_zero_sequence`](inc/damp/motor/modulation.hpp#L326) | Zero-sequence for any PwmScheme |
| [`motor_constant`](inc/damp/motor/foc.hpp#L206) | Motor constant Kₘ (torque per √copper-loss) — a figure of merit |
| [`mtpa_id_from_iq`](inc/damp/motor/mtpa.hpp#L49) | MTPA d-axis current on the trajectory for a given q-axis current |
| [`mtpa_reference`](inc/damp/motor/mtpa.hpp#L85) | MTPA dq current reference for a commanded torque |
| [`nameplate_design`](inc/damp/motor/drive_design.hpp#L232) | Nameplate → MTPA limits + foc_cascade_tune + servo + field-weakening configs |
| [`nameplate_limits`](inc/damp/motor/ipm.hpp#L301) | Voltage/current/base-speed/MTPA peak limits from bus and machine nameplate |
| [`park_transform`](inc/damp/transforms.hpp#L390) | Park transform (αβ → dq) |
| [`park_zero_transform`](inc/damp/transforms.hpp#L494) | Park transform with zero passthrough (αβ0 → dq0) |
| [`pmsm_config`](inc/damp/motor/servo.hpp#L169) | Split a PmsmServoConfig into a PmsmConfig for the motor stage |
| [`pmsm_nameplate_design`](inc/damp/motor/drive_design.hpp#L134) | Nameplate → limits + foc_cascade_tune + filled motor::PmsmServoConfig |
| [`pwm_center_aligned_ccr`](inc/damp/motor/inverter.hpp#L94) | Single CCR for up-down centre-aligned PWM (ARR = peak count) |
| [`pwm_center_aligned_ccrs`](inc/damp/motor/inverter.hpp#L133) | Three single-CCR levels for up-down centre-aligned PWM |
| [`pwm_center_aligned_compare`](inc/damp/motor/inverter.hpp#L69) | Dual-edge centre-aligned compare levels (full carrier period in ticks) |
| [`pwm_center_aligned_compare_deadtime`](inc/damp/motor/inverter.hpp#L151) | Dual-edge compares with delayed turn-on (software dead-time on rise) |
| [`pwm_center_aligned_compares`](inc/damp/motor/inverter.hpp#L113) | Dual-edge compares for three legs (full-period tick count) |
| [`pwm_center_aligned_edges`](inc/damp/motor/inverter.hpp#L42) | Centre-aligned ideal high window for duty d in period Ts |
| [`rotary_transformer_primary_voltage`](inc/damp/motor/sepex.hpp#L182) | Primary voltage that produces secondary field voltage vf_sec under ideal RT |
| [`rotational_load_ss`](inc/damp/motor/mechanical_estimator.hpp#L54) | Continuous state-space model of a 1-DOF rotational drivetrain with an augmented load-torque state |
| [`rotor_time_constant`](inc/damp/motor/acim.hpp#L90) | Rotor time constant T_r = L_r / R_r [s] |
| [`six_step_duty_cycles`](inc/damp/motor/modulation.hpp#L465) | Classical six-step (full-wave) duty pattern from the αβ angle |
| [`six_step_fundamental_voltage`](inc/damp/motor/modulation.hpp#L149) | Six-step (square-wave) fundamental peak phase voltage |
| [`slip_frequency`](inc/damp/motor/acim.hpp#L147) | Slip frequency [elec rad/s] for IFOC |
| [`smo`](inc/damp/motor/sensorless_smo.hpp#L67) | Design SMO gains from SPM electrical parameters |
| [`steady_state_vdq`](inc/damp/motor/pmsm_equations.hpp#L47) | Steady-state dq voltage at an operating point (rotor-synchronous frame) |
| [`steady_state_voltage_magnitude`](inc/damp/motor/pmsm_equations.hpp#L64) | \|Vdq\| magnitude at a steady-state operating point |
| [`suggest_cascade_bandwidths`](inc/damp/motor/foc_autotune.hpp#L60) | Nested cascade bandwidths from a current-loop target (rad/s) |
| [`svm_duty_cycles`](inc/damp/motor/modulation.hpp#L442) | Space-vector PWM duty cycles from an αβ voltage command |
| [`svpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L203) | Min-max zero-sequence injection for space-vector PWM |
| [`symmetrical_components`](inc/damp/transforms.hpp#L653) | Forward symmetrical-component (Fortescue) transform (abc → 012) |
| [`thipwm_zero_sequence`](inc/damp/motor/modulation.hpp#L229) | Third-harmonic injection (THIPWM) zero-sequence — 1/6 of fundamental |
| [`torque`](inc/damp/motor/acim.hpp#L120) | Electromagnetic torque from rotor flux and q-axis current (+3 more overloads) |
| [`torque_constant`](inc/damp/motor/acim.hpp#L111) | Torque constant at a rotor flux (amplitude-invariant Clarke/Park) (+1 more overload) |
| [`torque_constant_from_flux`](inc/damp/motor/foc.hpp#L103) | Torque constant Kₜ of a PMSM (amplitude-invariant convention) |
| [`torque_constant_from_Kv`](inc/damp/motor/foc.hpp#L165) | Torque constant from the datasheet velocity constant Kᵥ |
| [`torque_speed_envelope`](inc/damp/motor/ipm.hpp#L208) | IPM torque–speed envelope table (host / Lut1D bake) (+1 more overload) |
| [`transient_inductance`](inc/damp/motor/acim.hpp#L81) | Stator transient inductance L_σ = σ L_s [H] |
| [`two_mass_rotational_ss`](inc/damp/motor/two_mass_estimator.hpp#L98) | Continuous two-mass rotational drivetrain (compliant shaft) |
| [`two_mass_ss`](inc/damp/motor/two_mass_estimator.hpp#L135) | Alias of two_mass_rotational_ss (roadmap short name) |
| [`two_mass_torsional_wn`](inc/damp/motor/two_mass_estimator.hpp#L145) | Undamped two-mass torsional natural frequency [rad/s] |
| [`voltage_circle_radius`](inc/damp/motor/foc.hpp#L268) | Radius of the SVPWM voltage circle (max synthesizable \|V_dq\|) |

## Power electronics pack (if present)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AcPolarity`](inc/damp/power/totem_pole.hpp#L45) | AC line polarity for totem-pole slow (LF) leg |
| [`AnpcMap`](inc/damp/power/anpc.hpp#L37) | Single-phase ANPC leg map — voltage command → 6 duties (layer 3) |
| [`BoostMap`](inc/damp/power/boost.hpp#L42) | Boost CCM plant map — voltage command → duty (layer 3) |
| [`BuckBoostMap`](inc/damp/power/buck_boost.hpp#L31) | Inverting buck-boost CCM plant map — HardSwitchingConverter (layer 3) |
| [`BuckMap`](inc/damp/power/buck.hpp#L33) | Buck CCM plant map — voltage command → duty (layer 3) |
| [`CukMap`](inc/damp/power/cuk.hpp#L30) | Ćuk CCM plant map — HardSwitchingConverter (layer 3) |
| [`CurrentRegulatedStage`](inc/damp/power/concepts.hpp#L124) | Current-regulated PE stage: control(r, y) → switch duties (layer 2) |
| [`CurrentStage`](inc/damp/power/current_stage.hpp#L58) | Current-regulated stage: PIController + map |
| [`FourSwitchBuckBoostMap`](inc/damp/power/four_switch_buck_boost.hpp#L33) | Four-switch buck–boost region map — voltage → Q1/Q3 duties (layer 3) |
| [`FullBridge1phMap`](inc/damp/power/full_bridge.hpp#L56) | 1φ full-bridge map — voltage command → two high-side duties (layer 3) |
| [`FullBridgePwm`](inc/damp/power/full_bridge.hpp#L33) | 1φ full-bridge carrier style |
| [`HalfBridge1phMap`](inc/damp/power/half_bridge.hpp#L40) | Half-bridge pole map — voltage command → one high-side duty (layer 3) |
| [`HardSwitchingConverter`](inc/damp/power/concepts.hpp#L94) | Voltage-shaped command → duties (topology / plant map plug) |
| [`HasSwitchDuties`](inc/damp/power/concepts.hpp#L76) | Duck-typed switch command: duty vector + clip flag |
| [`NpcMap`](inc/damp/power/npc.hpp#L38) | Single-phase diode NPC leg map — voltage command → 4 duties (layer 3) |
| [`SepicMap`](inc/damp/power/sepic.hpp#L29) | SEPIC CCM plant map — HardSwitchingConverter (layer 3) |
| [`ThreeLevelDwells`](inc/damp/power/multilevel.hpp#L34) | P/O/N dwell fractions for one pole (sum to 1 in the linear region) |
| [`TotemPoleMap`](inc/damp/power/totem_pole.hpp#L56) | Totem-pole bridgeless map — boost plant + slow-leg polarity fan-out |
| [`TtypeMap`](inc/damp/power/ttype.hpp#L42) | Single-phase T-type (TNPC) leg map — voltage command → 4 duties (layer 3) |
| [`ZetaMap`](inc/damp/power/zeta.hpp#L29) | Zeta CCM plant map — HardSwitchingConverter (layer 3) |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`three_level_modulator_1ph`](inc/damp/power/multilevel.hpp#L102) | 1φ unit modulator from pole-to-mid voltage command + NP balance |
| [`three_level_np_zero_sequence`](inc/damp/power/multilevel.hpp#L73) | Neutral-point zero-sequence offset on unit modulator m |
| [`three_level_pd_dwells`](inc/damp/power/multilevel.hpp#L44) | Phase-disposition (PD) dwell map for one pole |
| [`three_level_pole_to_mid`](inc/damp/power/multilevel.hpp#L128) | Pole-to-mid: balanced bus → P=+Vdc/2, O=0, N=−Vdc/2 |
| [`three_level_pole_to_neg`](inc/damp/power/multilevel.hpp#L120) | Average pole-to-negative voltage: N=0, O=v_mid, P=Vdc |

> **Embedded helpers — controls-adjacent utilities (not ETL)**

## Embedded helpers (controls-adjacent utilities)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AdaptiveLut1D`](inc/damp/toolbox/lookup.hpp#L427) | Adaptive 1-D lookup: fixed breakpoints, online-updated cell values |
| [`AdaptiveLut2D`](inc/damp/toolbox/lookup.hpp#L502) | Adaptive 2-D lookup: fixed grid, online-updated cells (nearest grid node) |
| [`AdaptOutOfRange`](inc/damp/toolbox/lookup.hpp#L48) | How AdaptiveLut1D / AdaptiveLut2D treat samples outside the breakpoint span |
| [`AffineCal`](inc/damp/toolbox/scaling.hpp#L83) | Affine sensor calibration `y = gain·x + offset` |
| [`AnalogCrossMode`](inc/damp/toolbox/conditioning.hpp#L498) | How two redundant analog channels should relate in engineering units |
| [`AnalogCrossResult`](inc/damp/toolbox/conditioning.hpp#L520) | Output of cross_check_analog |
| [`AnalogCrossStatus`](inc/damp/toolbox/conditioning.hpp#L506) | Result of comparing two redundant analog channels |
| [`AnalogInput`](inc/damp/toolbox/io.hpp#L46) | A single analog input: range/fault check on the raw reading, then affine calibration to engineering units |
| [`AngleUnwrapper`](inc/damp/toolbox/encoder.hpp#L190) | Wrap-safe accumulator: successive wrapped absolute readings → continuous position |
| [`AxisInput`](inc/damp/toolbox/io.hpp#L191) | Operator-axis conditioning chain (joystick / RC stick → command) |
| [`BLINK`](inc/damp/toolbox/iec61131.hpp#L550) | BLINK (free-running square-wave / flasher) |
| [`Bounds`](inc/damp/toolbox/bounds.hpp#L43) | A per-channel closed-interval box constraint |
| [`Button`](inc/damp/toolbox/io.hpp#L238) | Debounced momentary push-button with edge and long-press detection |
| [`ConstantInertiaFeedforward`](inc/damp/toolbox/actuator.hpp#L197) | Per-axis decoupled torque feedforward: `τ = J·a + b·v + τ_c·sign(v) + g` |
| [`Counter`](inc/damp/toolbox/logic.hpp#L273) | Edge-counting up/down counter: increments on each rising edge of up, decrements on each rising edge of down. Returns the running count |
| [`CTD`](inc/damp/toolbox/iec61131.hpp#L364) | CTD Counter (Count Down) |
| [`CTU`](inc/damp/toolbox/iec61131.hpp#L324) | CTU Counter (Count Up) |
| [`CTUD`](inc/damp/toolbox/iec61131.hpp#L404) | CTUD Counter (Count Up Down) |
| [`Debounce`](inc/damp/toolbox/logic.hpp#L213) | Debounce: the output adopts in only after in differs from the current output continuously for stable_time. Rejects contact bounce and brief glitches. (Not an IEC block — the one everyone hand-rolls.) |
| [`DFF`](inc/damp/toolbox/iec61131.hpp#L461) | D Flip-Flop (edge-triggered data latch) |
| [`DLATCH`](inc/damp/toolbox/iec61131.hpp#L492) | D Latch (level-sensitive / transparent latch) |
| [`DrivePair`](inc/damp/toolbox/io.hpp#L315) | Left/right actuator pair (differential/arcade drive output) |
| [`Extrapolation`](inc/damp/toolbox/lookup.hpp#L36) | Out-of-range behaviour for a Lut1D / Lut2D / Lut3D query beyond its breakpoints |
| [`F_TRIG`](inc/damp/toolbox/iec61131.hpp#L153) | F_TRIG (Falling Edge Trigger) |
| [`FallingEdge`](inc/damp/toolbox/logic.hpp#L44) | Falling-edge detector: true on the tick x goes true → false |
| [`HolonomicOutput`](inc/damp/toolbox/io.hpp#L463) | Four-wheel actuator set (mecanum/holonomic drive output), X-config |
| [`Hysteresis`](inc/damp/toolbox/conditioning.hpp#L332) | Hysteresis comparator (Schmitt trigger): bool output with separate on/off thresholds to reject chatter |
| [`Interpolation`](inc/damp/toolbox/lookup.hpp#L42) | In-range interpolant for Lut1D (2-D stays bilinear; use SplineSurface for C¹ grids) |
| [`Latch`](inc/damp/toolbox/logic.hpp#L62) | Set/reset latch. @tparam SetDominant which input wins when both are asserted (default: set-dominant, e.g. a trip overriding a clear) |
| [`Lut1D`](inc/damp/toolbox/lookup.hpp#L108) | 1-D interpolating lookup table over monotonic breakpoints |
| [`Lut2D`](inc/damp/toolbox/lookup.hpp#L200) | 2-D bilinear interpolating lookup table over a regular grid |
| [`Lut3D`](inc/damp/toolbox/lookup.hpp#L274) | 3-D trilinear interpolating lookup table over a regular grid |
| [`OffDelayTimer`](inc/damp/toolbox/logic.hpp#L133) | Off-delay timer: output goes true immediately when in is true and stays true until in has been false continuously for delay |
| [`OnDelayTimer`](inc/damp/toolbox/logic.hpp#L96) | On-delay timer: output goes true once in has been held true continuously for delay; drops immediately when in goes false |
| [`Periodic`](inc/damp/toolbox/timing.hpp#L125) | Periodic trigger — fires once per elapsed period |
| [`ProportionalCurrentMap`](inc/damp/toolbox/io.hpp#L139) | Affine map from a normalized command to a proportional current (mA) |
| [`PulseTimer`](inc/damp/toolbox/logic.hpp#L170) | Pulse timer (non-retriggerable): a rising edge of in emits a fixed |
| [`QuadMode`](inc/damp/toolbox/encoder.hpp#L47) | Quadrature decode resolution (edges counted per A/B cycle) |
| [`QuadratureDecoder`](inc/damp/toolbox/encoder.hpp#L69) | Software A/B quadrature decoder with optional index |
| [`R_TRIG`](inc/damp/toolbox/iec61131.hpp#L118) | R_TRIG (Rising Edge Trigger) |
| [`RangeMonitor`](inc/damp/toolbox/conditioning.hpp#L430) | Analog-input range/fault monitor (NAMUR NE43 pattern) |
| [`RedundantAnalogPair`](inc/damp/toolbox/io.hpp#L90) | Dual-channel analog input with same-slope or opposite-slope cross-check |
| [`RisingEdge`](inc/damp/toolbox/logic.hpp#L30) | Rising-edge detector: true on the tick x goes false → true |
| [`RS`](inc/damp/toolbox/iec61131.hpp#L91) | RS Latch (Reset-Set Latch) |
| [`ServoAxis`](inc/damp/toolbox/actuator.hpp#L107) | One servoactuator transmission: SI joint unit ⟷ drive (motor) units |
| [`ServoBank`](inc/damp/toolbox/actuator.hpp#L229) | A bank of ServoAxis transmissions: maps a synchronized multi-axis TrajectoryState array straight to drive commands in one call |
| [`ServoCommand`](inc/damp/toolbox/actuator.hpp#L84) | A drive-native servoactuator setpoint: position, velocity, torque |
| [`SignalStatus`](inc/damp/toolbox/conditioning.hpp#L368) | Classification of an analog input against its valid/fault bands |
| [`SlewLimiter`](inc/damp/toolbox/conditioning.hpp#L268) | Slew-rate limiter: bound how fast the output may follow the target |
| [`SplineSurface`](inc/damp/toolbox/lookup.hpp#L393) | 2-D interpolating surface with smooth (Catmull-Rom spline) blending |
| [`SR`](inc/damp/toolbox/iec61131.hpp#L63) | SR Latch (Set-dominant Set-Reset Latch) |
| [`Stopwatch`](inc/damp/toolbox/timing.hpp#L43) | Free-running elapsed-time accumulator |
| [`Switch`](inc/damp/toolbox/io.hpp#L282) | Debounced maintained switch (toggle/selector contact) with change flag |
| [`Tachometer`](inc/damp/toolbox/encoder.hpp#L144) | Pulse-based speed (tachometer) with frequency/period crossover |
| [`TFF`](inc/damp/toolbox/iec61131.hpp#L517) | T Flip-Flop (toggle on rising edge) |
| [`ThermalKalmanObserver`](inc/damp/toolbox/thermal.hpp#L248) | Runtime thermal-network Kalman observer |
| [`Thermistor`](inc/damp/toolbox/thermistor.hpp#L166) | NTC thermistor linearization (resistance → temperature) |
| [`ThermistorCoeffs`](inc/damp/toolbox/thermistor.hpp#L41) | Fitted NTC coefficients in Steinhart-Hart form |
| [`Timeout`](inc/damp/toolbox/timing.hpp#L78) | One-shot timeout |
| [`TOF`](inc/damp/toolbox/iec61131.hpp#L226) | TOF Timer (Timer Off Delay) |
| [`Toggle`](inc/damp/toolbox/logic.hpp#L245) | Toggle (T flip-flop): output flips on each rising edge of in |
| [`TON`](inc/damp/toolbox/iec61131.hpp#L183) | TON Timer (Timer On Delay) |
| [`TP`](inc/damp/toolbox/iec61131.hpp#L271) | TP Timer (Timer Pulse) |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`arcade_drive`](inc/damp/toolbox/io.hpp#L361) | Arcade (single-stick) drive mixer: x (steer) + y (throttle) → left/right |
| [`beta`](inc/damp/toolbox/thermistor.hpp#L77) | Fit NTC coefficients from the Beta-parameter model |
| [`catmull_1d`](inc/damp/toolbox/lookup.hpp#L339) | 1-D non-uniform Catmull-Rom (cubic Hermite) evaluation |
| [`cauer_thermal_ss_ambient`](inc/damp/toolbox/thermal.hpp#L135) | Continuous Cauer RC ladder in absolute temperature with ambient input |
| [`cauer_thermal_ss_ambient_mimo`](inc/damp/toolbox/thermal.hpp#L162) | Continuous Cauer RC ladder, absolute temps, all nodes as outputs |
| [`classify_range`](inc/damp/toolbox/conditioning.hpp#L396) | Classify x against the four band edges `[fault_lo (valid_lo, valid_hi) fault_hi]` (assumed ordered, non-decreasing) |
| [`compensate_coulomb_viscous`](inc/damp/toolbox/conditioning.hpp#L178) | Coulomb + viscous friction compensation (static inverse, massless) |
| [`compensate_stribeck`](inc/damp/toolbox/conditioning.hpp#L205) | Stribeck-style friction compensation (static + Coulomb + viscous) |
| [`cross_check_analog`](inc/damp/toolbox/conditioning.hpp#L573) | Cross-check two calibrated analog readings for redundant sensors |
| [`deadband`](inc/damp/toolbox/conditioning.hpp#L47) | Dead zone over `[lower, upper]`, matching Simulink®'s Dead Zone block (+1 more overload) |
| [`differential_drive`](inc/damp/toolbox/io.hpp#L340) | Differential (tank) drive mixer: throttle + turn → left/right |
| [`expo`](inc/damp/toolbox/conditioning.hpp#L252) | Exponential response curve `y = (1−k)·x + k·x³` (RC "expo") |
| [`h_pattern`](inc/damp/toolbox/io.hpp#L377) | H-pattern (dual-lever / tank) drive: one lever per track, no mixing |
| [`holonomic_drive`](inc/damp/toolbox/io.hpp#L495) | Four-wheel holonomic drive mixer (mecanum / omni) |
| [`inverse_deadband`](inc/damp/toolbox/conditioning.hpp#L79) | Inverse dead zone: add an offset to overcome a physical dead zone (valve overlap, static friction, motor stiction), with independent negative/positive offsets (+1 more overload) |
| [`inverse_deadband_sine`](inc/damp/toolbox/conditioning.hpp#L152) | Sine-shaped inverse dead zone (full boost once \|x\|≥w) |
| [`inverse_deadband_soft`](inc/damp/toolbox/conditioning.hpp#L128) | Smooth inverse dead zone via soft-sign (tanh-like Coulomb boost) |
| [`inverse_lerp`](inc/damp/toolbox/scaling.hpp#L48) | Inverse of lerp: the fraction t such that `lerp(a, b, t) == x` |
| [`is_fault`](inc/damp/toolbox/conditioning.hpp#L382) | True for a wire fault (FaultLow/FaultHigh) — i.e. not a real reading at all |
| [`is_usable`](inc/damp/toolbox/conditioning.hpp#L539) | True when a channel reading may still be used (in-span or saturated) |
| [`is_valid`](inc/damp/toolbox/conditioning.hpp#L377) | True only for the in-span status |
| [`iso_c_drive`](inc/damp/toolbox/io.hpp#L412) | ISO-C (single-stick, coordinated / curvature) drive |
| [`iso_s_drive`](inc/damp/toolbox/io.hpp#L393) | ISO-S (single-stick, speed-summed) drive — arcade travel stick |
| [`lerp`](inc/damp/toolbox/scaling.hpp#L35) | Linear interpolation between a and b by fraction t |
| [`LIMIT`](inc/damp/toolbox/iec61131.hpp#L598) | LIMIT (IEC 61131-3 selection function): clamp in to [mn, mx] |
| [`linear_screw`](inc/damp/toolbox/actuator.hpp#L166) | Build a ServoAxis for a linear axis driven by a leadscrew/belt |
| [`lut_segment`](inc/damp/toolbox/lookup.hpp#L63) | Index of the interpolation segment containing x |
| [`mecanum_drive`](inc/damp/toolbox/io.hpp#L514) | Mecanum (H-layout, 45° rollers) mixer — alias for holonomic_drive |
| [`MUX`](inc/damp/toolbox/iec61131.hpp#L617) | MUX (IEC 61131-3 multiplexer): select input k of N (0-based) |
| [`omniwheel_drive`](inc/damp/toolbox/io.hpp#L520) | Omni-wheel (X-layout) mixer — alias for holonomic_drive |
| [`poly_horner`](inc/damp/toolbox/scaling.hpp#L134) | Evaluate a polynomial at x by Horner's method |
| [`rescale`](inc/damp/toolbox/scaling.hpp#L69) | Affine map of x from the input range to the output range |
| [`rotary_gearbox`](inc/damp/toolbox/actuator.hpp#L149) | Build a ServoAxis for a rotary joint behind a gearbox |
| [`scaled_deadband`](inc/damp/toolbox/conditioning.hpp#L233) | Center dead zone that rescales the surviving range back to full span |
| [`SEL`](inc/damp/toolbox/iec61131.hpp#L606) | SEL (IEC 61131-3 binary selection): g ? in1 : in0 |
| [`soft_sign`](inc/damp/toolbox/conditioning.hpp#L103) | Smooth unit direction x / √(x²+ε²) ∈ (−1,1) |
| [`steer_map`](inc/damp/toolbox/io.hpp#L452) | Curve-driven steering map — arbitrary `(throttle, turn)` geometry |
| [`steinhart_hart`](inc/damp/toolbox/thermistor.hpp#L118) | Fit the Steinhart-Hart coefficients from three calibration points |
| [`thermal_kalman_plant`](inc/damp/toolbox/thermal.hpp#L195) | Prepare a discretized thermal plant for design::kalman |
| [`two_point_cal`](inc/damp/toolbox/scaling.hpp#L111) | Fit an AffineCal through two `(raw, engineering)` points |
| [`wrapped_delta`](inc/damp/toolbox/encoder.hpp#L41) | Signed difference between two unsigned counter readings, wrap-safe |

> **Host — analysis, sim, MATLAB® aliases**

## Frequency-domain analysis (host)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AdaptationGateDecision`](inc/damp/analysis/identification.hpp#L226) | Gate decision for whether adaptive updates are currently allowed |
| [`ARXModel`](inc/damp/analysis/identification.hpp#L103) | ARX candidate model with fixed numerator/denominator orders |
| [`BodeResult`](inc/damp/analysis/frequency.hpp#L60) | Bode plot data for a SISO system |
| [`DriftMonitorResult`](inc/damp/analysis/identification.hpp#L215) | Result of model-drift monitoring used for adaptation triggers |
| [`FitMetrics`](inc/damp/analysis/identification.hpp#L138) | Scalar fit metrics used for model comparison and gating |
| [`FOPDTModel`](inc/damp/analysis/identification.hpp#L76) | First-order plus dead-time candidate model |
| [`FrequencyPoint`](inc/damp/analysis/frequency.hpp#L38) | Single-point frequency response result |
| [`FrequencyResponsePoint`](inc/damp/analysis/identification.hpp#L52) | One frequency-response point used by sweep/chirp-based identification |
| [`FrequencyResponseSummary`](inc/damp/analysis/identification.hpp#L64) | Summary statistics for a frequency-response data set |
| [`FRFEstimateResult`](inc/damp/analysis/identification.hpp#L116) | Frequency-response estimate quality summary |
| [`GreyBoxCandidateModels`](inc/damp/analysis/identification.hpp#L172) | Holds all candidate models produced during an identification pass |
| [`GreyBoxIdentificationResult`](inc/damp/analysis/identification.hpp#L188) | Selected model used by downstream model-based controller design |
| [`IdentificationHandoff`](inc/damp/analysis/identification.hpp#L202) | Aggregated handoff payload from identification to model-based control design |
| [`IdentificationModelKind`](inc/damp/analysis/identification.hpp#L161) | Enumerates which reduced model family is selected for downstream design |
| [`IdentifiedModelLike`](inc/damp/analysis/identification.hpp#L237) | Concept for identified models consumed by downstream design modules |
| [`ImpedanceResult`](inc/damp/analysis/frequency.hpp#L872) | Result of impedance frequency response evaluation |
| [`InputOutputSample`](inc/damp/analysis/identification.hpp#L24) | One time-aligned sample used for closed-loop or open-loop identification |
| [`LoopResponseResult`](inc/damp/analysis/frequency.hpp#L415) | Open-loop and closed-loop frequency response package |
| [`LoopSummary`](inc/damp/analysis/frequency.hpp#L448) | Compact loop summary metrics for quick stability/robustness checks |
| [`LsimInfo`](inc/damp/analysis/time_response.hpp#L484) | Transient characteristics of an arbitrary response signal |
| [`LsimResult`](inc/damp/analysis/time_response.hpp#L66) | Result of a single-trajectory simulation: time, output, and state history |
| [`MiddlebrookResult`](inc/damp/analysis/frequency.hpp#L896) | Result of Middlebrook minor loop gain analysis |
| [`NicholsPoint`](inc/damp/analysis/frequency.hpp#L654) | Single-point Nichols chart sample (open-loop phase vs magnitude) |
| [`NicholsResult`](inc/damp/analysis/frequency.hpp#L666) | Nichols chart locus across a frequency sweep |
| [`NyquistPoint`](inc/damp/analysis/frequency.hpp#L369) | Single-point Nyquist response data |
| [`NyquistResult`](inc/damp/analysis/frequency.hpp#L381) | Nyquist response data across a frequency sweep |
| [`PoleInfo`](inc/damp/analysis/poles.hpp#L75) | Natural frequency and damping ratio for each pole |
| [`PoleZeroMap`](inc/damp/analysis/poles.hpp#L115) | Poles and zeros of a system, for pole-zero plotting |
| [`RootLocusResult`](inc/damp/analysis/poles.hpp#L203) | Root-locus data: closed-loop poles along a gain grid |
| [`SigmaPoint`](inc/damp/analysis/frequency.hpp#L739) | Singular values of G(jω) at one frequency |
| [`SigmaResult`](inc/damp/analysis/frequency.hpp#L753) | Singular-value frequency response over a grid |
| [`SOPDTModel`](inc/damp/analysis/identification.hpp#L89) | Second-order plus dead-time candidate model |
| [`StepInfo`](inc/damp/analysis/time_response.hpp#L354) | Step-response characteristics of a single output signal |
| [`StepResponseSummary`](inc/damp/analysis/identification.hpp#L36) | Compact statistics extracted from a step experiment |
| [`SweepIdentificationResult`](inc/damp/analysis/identification.hpp#L127) | End-to-end output of a sweep-based identification pass |
| [`TimeResponse`](inc/damp/analysis/time_response.hpp#L51) | Multi-channel time-domain response sampled on a time grid |
| [`ValidationResult`](inc/damp/analysis/identification.hpp#L150) | Validation status for a model against held-out or replayed data |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`arange`](inc/damp/analysis/linspace.hpp#L134) | Half-open arithmetic range [start, stop) with step step (+2 more overloads) |
| [`bode`](inc/damp/analysis/frequency.hpp#L250) | Compute Bode plot data for a SISO state-space system (+2 more overloads) |
| [`bode_discrete`](inc/damp/analysis/frequency.hpp#L358) | Compute Bode plot data for a discrete-time SISO state-space system |
| [`canonical_phase_margin`](inc/damp/analysis/frequency.hpp#L148) | Normalize phase margin to (-180, 180] |
| [`damp`](inc/damp/analysis/poles.hpp#L90) | Compute natural frequency and damping for each pole |
| [`dcgain`](inc/damp/analysis/norms.hpp#L40) | Compute DC gain of a continuous-time system |
| [`frf_s_or_z`](inc/damp/analysis/frequency.hpp#L974) | FRF evaluation point: jω (continuous) or e^{jωTs} (discrete) |
| [`gain_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L202) | Find gain margin using unwrapped phase trajectory |
| [`geomspace`](inc/damp/analysis/linspace.hpp#L88) | Geometric sequence from start to end (endpoint values, not exponents) |
| [`impedance`](inc/damp/analysis/frequency.hpp#L999) | Compute impedance frequency response from a SISO admittance system |
| [`impedance_direct`](inc/damp/analysis/frequency.hpp#L1039) | Compute impedance frequency response from a SISO impedance transfer function |
| [`impulse`](inc/damp/analysis/time_response.hpp#L185) | Impulse response of a (MIMO) state-space system (+1 more overload) |
| [`initial`](inc/damp/analysis/time_response.hpp#L219) | Initial-condition (free) response of a (MIMO) state-space system |
| [`is_stable_continuous`](inc/damp/analysis/poles.hpp#L58) | Check continuous-time stability |
| [`loop_metrics`](inc/damp/analysis/frequency.hpp#L480) | One-call loop analysis: compute L/S/T response and return compact metrics (+1 more overload) |
| [`loop_response`](inc/damp/analysis/frequency.hpp#L558) | Compute open-loop L, sensitivity S, complementary sensitivity T, and Nyquist data (+1 more overload) |
| [`lsim`](inc/damp/analysis/time_response.hpp#L257) | Forced time response of a (MIMO) state-space system to an input signal (+1 more overload) |
| [`lsiminfo`](inc/damp/analysis/time_response.hpp#L503) | Compute transient characteristics from an output/time signal |
| [`middlebrook`](inc/damp/analysis/frequency.hpp#L1083) | Middlebrook stability analysis for cascaded source-load systems (+1 more overload) |
| [`nichols`](inc/damp/analysis/frequency.hpp#L683) | Build Nichols points from existing Bode data (+2 more overloads) |
| [`norm_h2`](inc/damp/analysis/norms.hpp#L76) | H2 norm of a state-space system |
| [`norm_hinf`](inc/damp/analysis/norms.hpp#L134) | H∞ norm of a state-space system: sup_ω σ̄(G(jω)) |
| [`nyquist`](inc/damp/analysis/frequency.hpp#L497) | Compute Nyquist data for a SISO state-space system (+1 more overload) |
| [`phase_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L166) | Find phase margin using unwrapped phase trajectory |
| [`poles`](inc/damp/analysis/poles.hpp#L41) | Compute open-loop poles (eigenvalues of A matrix) |
| [`poly_roots`](inc/damp/analysis/poles.hpp#L128) | Roots of a polynomial given in ascending powers (MATLAB® `roots`, reversed order) |
| [`pzmap`](inc/damp/analysis/poles.hpp#L154) | Pole-zero map of a SISO transfer function (MATLAB® `pzmap(tf)`) (+2 more overloads) |
| [`rlocus`](inc/damp/analysis/poles.hpp#L262) | Root locus of a SISO state-space plant over an explicit gain grid (+1 more overload) |
| [`sigma`](inc/damp/analysis/frequency.hpp#L796) | Singular-value frequency response of a (possibly MIMO) state-space system (+1 more overload) |
| [`step`](inc/damp/analysis/time_response.hpp#L153) | Step response of a (MIMO) state-space system (+1 more overload) |
| [`stepinfo`](inc/damp/analysis/time_response.hpp#L374) | Compute step-response characteristics from an output/time signal (+1 more overload) |
| [`summarize_loop_response`](inc/damp/analysis/frequency.hpp#L460) | Summarize loop_response() results into one compact metrics struct |
| [`unwrap_phase_deg`](inc/damp/analysis/frequency.hpp#L125) | Unwrap phase data in degrees to avoid +/-180 discontinuities |

## Simulation / SIL harness (host)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AdaptiveOptions`](inc/damp/simulation/solver.hpp#L306) | Options for adaptive-step ODE integration |
| [`AdaptiveStepIntegrator`](inc/damp/simulation/integrator.hpp#L65) | Integrator that reports a genuine embedded local-error estimate |
| [`AdaptiveStepSolver`](inc/damp/simulation/solver.hpp#L358) | Adaptive-step ODE solver (+1 more overload) |
| [`BackwardEuler`](inc/damp/simulation/integrator.hpp#L295) | Backward Euler integrator |
| [`BDF2`](inc/damp/simulation/integrator.hpp#L366) | Backward Differentiation Formula 2 (BDF2) integrator |
| [`CachedZoh`](inc/damp/simulation/cached_zoh.hpp#L62) | Cached discrete ZOH maps Ad, Bd for step size h |
| [`Discrete`](inc/damp/simulation/integrator.hpp#L83) | Discrete-time integrator (no integration, just one step) |
| [`DP45`](inc/damp/simulation/integrator.hpp#L757) | Dormand-Prince 5(4) adaptive integrator — the ode45 pair |
| [`ErrorTolerance`](inc/damp/simulation/solver.hpp#L259) | Scalar-or-vector local-error tolerance for adaptive ODE control |
| [`EventSimConfig`](inc/damp/simulation/hybrid.hpp#L67) | Configuration for event-driven hybrid / multi-rate style runs |
| [`Exact`](inc/damp/simulation/integrator.hpp#L117) | Exact integrator for LTI systems |
| [`FixedStepSolver`](inc/damp/simulation/solver.hpp#L118) | Fixed-step ODE solver (+1 more overload) |
| [`ForwardEuler`](inc/damp/simulation/integrator.hpp#L169) | Forward Euler integrator |
| [`Heun`](inc/damp/simulation/integrator.hpp#L840) | Heun's method (Improved Euler, RK2) integrator |
| [`HybridEventAction`](inc/damp/simulation/hybrid.hpp#L102) | Action applied when an event fires (after Exact advance to the event time) |
| [`HybridSimulationResult`](inc/damp/simulation/hybrid.hpp#L80) | Result of a piecewise-LTI hybrid simulation |
| [`IntegrationResult`](inc/damp/simulation/integrator.hpp#L31) | Result of an integration step |
| [`LogPolicy`](inc/damp/simulation/hybrid.hpp#L58) | How often to record samples along a hybrid trajectory |
| [`MultiRateConfig`](inc/damp/simulation/multirate.hpp#L81) | Configuration for multi-rate closed-loop runs |
| [`PiecewiseLTI`](inc/damp/simulation/hybrid.hpp#L94) | Plant: NModes continuous LTI pieces (same state/input/output dimensions) |
| [`RK23`](inc/damp/simulation/integrator.hpp#L902) | Bogacki-Shampine 2(3) adaptive integrator |
| [`RK3`](inc/damp/simulation/integrator.hpp#L969) | Classical 3rd-order Runge-Kutta (RK3) integrator |
| [`RK4`](inc/damp/simulation/integrator.hpp#L683) | Classical 4th-order Runge-Kutta (RK4) integrator |
| [`SimulationResult`](inc/damp/simulation/simulate.hpp#L103) | Result of a closed-loop simulation |
| [`SisoReferenceAdapter`](inc/damp/simulation/simulate.hpp#L71) | Adapt a SISO `control(r, y)` controller for callables that want `u = f(y)` or `u = f(t, y)` (the older simulate contracts) |
| [`SolveResult`](inc/damp/simulation/solver.hpp#L59) | Result of an ODE solve operation |
| [`SymplecticEuler`](inc/damp/simulation/integrator.hpp#L229) | Semi-implicit (symplectic) Euler for mechanical systems |
| [`Trapezoidal`](inc/damp/simulation/integrator.hpp#L623) | Trapezoidal (Tustin) integrator |
| [`TRBDF2`](inc/damp/simulation/integrator.hpp#L497) | TR-BDF2 composite integrator — the stiff adaptive pair (ode23tb) |
| [`TwoRateSimulationResult`](inc/damp/simulation/multirate.hpp#L94) | Result of a two-rate cascade simulation |
| [`XyPlotOpts`](inc/damp/simulation/plot_plotly.hpp#L535) | Optional axis framing for plot_xy |
| [`XySeries`](inc/damp/simulation/plot_plotly.hpp#L517) | One series for a planar (x, y) scatter / path plot |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`adaptive_solve`](inc/damp/simulation/solver.hpp#L774) | Adaptive-step solve in one call |
| [`apply_wet_plotly_shell`](inc/damp/simulation/plot_plotly.hpp#L64) | Apply damp's full-window HTML shell to plotlypp-generated markup |
| [`bodemag`](inc/damp/simulation/plot_plotly.hpp#L796) | Plot a magnitude-only Bode diagram (log frequency, dB magnitude) |
| [`bodeplot`](inc/damp/simulation/plot_plotly.hpp#L782) | Plot magnitude and phase Bode subplots |
| [`exact_lti_step`](inc/damp/simulation/hybrid.hpp#L315) | Exact open-loop step of a single continuous LTI system |
| [`extract_channel`](inc/damp/simulation/plot_plotly.hpp#L223) | Extract the i-th element from each vector entry into a std::vector<double> |
| [`fires_at`](inc/damp/simulation/multirate.hpp#L73) | True when base tick tick is a firing instant for divisor |
| [`fixed_solve`](inc/damp/simulation/solver.hpp#L243) | Fixed-step solve in one call |
| [`impulseplot`](inc/damp/simulation/plot_plotly.hpp#L733) | Plot an impulse response, one trace per input/output pair |
| [`locate_zero_crossing_exact`](inc/damp/simulation/hybrid.hpp#L335) | Locate a state zero-crossing of g(x) on an Exact LTI segment |
| [`lsimplot`](inc/damp/simulation/plot_plotly.hpp#L749) | Plot a forced (lsim) simulation, one trace per output |
| [`nicholsplot`](inc/damp/simulation/plot_plotly.hpp#L936) | Plot a Nichols chart (open-loop phase vs magnitude) with M-circle grid |
| [`nyquistplot`](inc/damp/simulation/plot_plotly.hpp#L821) | Plot a Nyquist locus with the -1 critical point marked |
| [`panel_color`](inc/damp/simulation/plot_plotly.hpp#L314) | Plotly default colorway entry for index i within a single subplot |
| [`panel_legend`](inc/damp/simulation/plot_plotly.hpp#L297) | Per-panel legend box just right of a stacked subplot band |
| [`panel_legend_id`](inc/damp/simulation/plot_plotly.hpp#L283) | Legend id for subplot row (1-based): "legend", "legend2", … |
| [`panel_x_domain`](inc/damp/simulation/plot_plotly.hpp#L273) | Horizontal domain for stacked multi-panel x-axes |
| [`plot_bode`](inc/damp/simulation/plot_plotly.hpp#L422) | Plot Bode magnitude and phase as subplots |
| [`plot_line`](inc/damp/simulation/plot_plotly.hpp#L491) | Simple line plot of time vs value |
| [`plot_simulation`](inc/damp/simulation/plot_plotly.hpp#L343) | Plot simulation results with subplots for states, outputs, and inputs |
| [`plot_step`](inc/damp/simulation/plot_plotly.hpp#L621) | Plot step response data |
| [`plot_xy`](inc/damp/simulation/plot_plotly.hpp#L556) | Multi-series planar scatter (true 2-D path / locus) (+1 more overload) |
| [`plots_root_for`](inc/damp/simulation/plot_plotly.hpp#L118) | Walk up from dir until a directory named "plots" is found |
| [`pzplot`](inc/damp/simulation/plot_plotly.hpp#L848) | Plot a pole-zero map on the complex plane (poles as ×, zeros as ○) |
| [`rlocusplot`](inc/damp/simulation/plot_plotly.hpp#L1033) | Plot a root locus (closed-loop poles vs gain) on the complex plane |
| [`sigmaplot`](inc/damp/simulation/plot_plotly.hpp#L990) | Plot singular-value frequency response (log frequency, dB) |
| [`simulate`](inc/damp/simulation/simulate.hpp#L132) | Simulate a nonlinear plant with a controller in closed loop |
| [`simulate_cached_zoh`](inc/damp/simulation/cached_zoh.hpp#L148) | Simulate a continuous LTI plant with fixed ZOH step and held input (+1 more overload) |
| [`simulate_cached_zoh_multirate`](inc/damp/simulation/cached_zoh.hpp#L227) | Multi-rate ZOH: held control samples, dense plant samples for plotting |
| [`simulate_discrete`](inc/damp/simulation/simulate.hpp#L478) | Simulate a discrete-time system with a controller |
| [`simulate_discrete_nonlinear`](inc/damp/simulation/simulate.hpp#L429) | Simulate a discrete-time nonlinear plant with a controller |
| [`simulate_lti`](inc/damp/simulation/simulate.hpp#L319) | Simulate a continuous LTI system with a controller |
| [`simulate_lti_siso`](inc/damp/simulation/simulate.hpp#L377) | SIL: continuous SISO LTI plant under a sampled SISO controller with constant reference |
| [`simulate_multirate`](inc/damp/simulation/multirate.hpp#L128) | Multi-rate closed-loop simulation with MultiRateConfig (+1 more overload) |
| [`simulate_piecewise_lti`](inc/damp/simulation/hybrid.hpp#L182) | Simulate a piecewise-LTI plant with Exact integration between events |
| [`simulate_sampled`](inc/damp/simulation/simulate.hpp#L205) | Simulate a continuous plant under a discrete (sampled) controller — multi-rate |
| [`simulate_state_feedback`](inc/damp/simulation/simulate.hpp#L264) | Simulate a nonlinear plant with state-feedback controller |
| [`simulate_two_rate`](inc/damp/simulation/multirate.hpp#L244) | Two-rate cascade with MultiRateConfig (+1 more overload) |
| [`simulate_two_rate_lti`](inc/damp/simulation/multirate.hpp#L388) | Two-rate cascade on a continuous LTI plant (A, B, C) |
| [`siso_ref`](inc/damp/simulation/simulate.hpp#L90) | Build a SisoReferenceAdapter for a SISO `control(r,y)` controller |
| [`stepplot`](inc/damp/simulation/plot_plotly.hpp#L717) | Plot a step response, one trace per input/output pair |
| [`to_double_vector`](inc/damp/simulation/plot_plotly.hpp#L236) | Convert std::vector\<T\> to std::vector<double> |
| [`to_std_vector`](inc/damp/simulation/plot_plotly.hpp#L211) | Convert a ColVec<N,T> to std::vector<double> for plotlypp |
| [`write_html`](inc/damp/simulation/plot_plotly.hpp#L145) | Write a plotlypp figure to HTML with damp's shell (prefer over writeHtml) |
| [`write_hybrid_trace_npy`](inc/damp/simulation/npy_export.hpp#L100) | Pack HybridSimulationResult into (N, 2+NX+NU+NY) float64: t, mode, x..., u..., y |
| [`write_npy_f64`](inc/damp/simulation/npy_export.hpp#L40) | Write a C-order float64 array as a .npy file |

## MATLAB®-style aliases (host)

**Blocks (structs, classes, enums, concepts)**

| Name | Description |
| ---- | ----------- |
| [`AllMarginResult`](inc/damp/matlab.hpp#L916) | All classical margins including delay margin (superset of margin) |
| [`MarginResult`](inc/damp/matlab.hpp#L868) | Gain/phase margins and their crossover frequencies |

**Functions**

| Name | Description |
| ---- | ----------- |
| [`acker`](inc/damp/matlab.hpp#L554) | Pole placement for state-feedback control |
| [`allmargin`](inc/damp/matlab.hpp#L937) | Gain, phase, and delay margins of a SISO loop over a frequency grid |
| [`bandwidth`](inc/damp/matlab.hpp#L1042) | -3 dB bandwidth of a SISO system over a frequency grid |
| [`blkdiag`](inc/damp/matlab.hpp#L349) | Block diagonal matrix construction |
| [`c2d`](inc/damp/matlab.hpp#L296) | MATLAB® interface function c2d to discretize a continuous-time state-space system (+1 more overload) |
| [`ctrb`](inc/damp/matlab.hpp#L258) | MATLAB® short alias for controllability_matrix (+1 more overload) |
| [`diag`](inc/damp/matlab.hpp#L372) | Returns a square diagonal matrix from the given array (+1 more overload) |
| [`dlqr`](inc/damp/matlab.hpp#L618) | Discrete-time Linear-Quadratic Regulator design |
| [`dlyap`](inc/damp/matlab.hpp#L1062) | MATLAB® alias for the discrete Lyapunov solve AXAᵀ−X+Q=0 |
| [`eig`](inc/damp/matlab.hpp#L473) | MATLAB® short alias for the eigenvalues of a square matrix |
| [`estim`](inc/damp/matlab.hpp#L490) | Form state estimator from system and estimator gain |
| [`eye`](inc/damp/matlab.hpp#L413) | Create an identity matrix of size n x n |
| [`gram`](inc/damp/matlab.hpp#L1109) | MATLAB® alias for the controllability/observability Gramian of a system |
| [`hinfnorm`](inc/damp/matlab.hpp#L1131) | MATLAB® alias for the H∞ system norm norm(sys,Inf) / hinfnorm(sys) |
| [`isstable`](inc/damp/matlab.hpp#L956) | Continuous-time stability predicate on a state matrix (+2 more overloads) |
| [`linmod`](inc/damp/matlab.hpp#L329) | MATLAB®-style nonlinear linearization about an operating point |
| [`lqg`](inc/damp/matlab.hpp#L737) | Linear-Quadratic-Gaussian regulator design |
| [`lqgreg`](inc/damp/matlab.hpp#L753) | Combine separate Kalman filter and LQR designs into an LQG controller |
| [`lqgtrack`](inc/damp/matlab.hpp#L765) | Linear-Quadratic-Gaussian design with integral action for tracking |
| [`lqi`](inc/damp/matlab.hpp#L724) | Linear-Quadratic Integral design for tracking |
| [`lqr`](inc/damp/matlab.hpp#L603) | Continuous-time LQR design (MATLAB®'s lqr) |
| [`lqrd`](inc/damp/matlab.hpp#L633) | Design discrete LQR from continuous-time system via discretization (+1 more overload) |
| [`lqry`](inc/damp/matlab.hpp#L680) | Output-weighted continuous LQR (state cost Q = Cᵀ Q_y C) (+1 more overload) |
| [`lyap`](inc/damp/matlab.hpp#L1052) | MATLAB® alias for the continuous Lyapunov solve AX+XAᵀ+Q=0 |
| [`make1DOF`](inc/damp/matlab.hpp#L220) | Force 1-DOF setpoint weights on a PID design result (b=c=1) |
| [`make2DOF`](inc/damp/matlab.hpp#L236) | Apply 2-DOF setpoint weights on a PID design result |
| [`margin`](inc/damp/matlab.hpp#L889) | Gain and phase margins of a SISO loop over a frequency grid |
| [`minreal`](inc/damp/matlab.hpp#L250) | MATLAB® short alias for design::minreal |
| [`norm`](inc/damp/matlab.hpp#L1122) | MATLAB® alias for the H2 system norm norm(sys,2) |
| [`null`](inc/damp/matlab.hpp#L459) | MATLAB® short alias for an orthonormal null-space basis |
| [`obsv`](inc/damp/matlab.hpp#L274) | MATLAB® short alias for observability_matrix (+1 more overload) |
| [`pade`](inc/damp/matlab.hpp#L998) | First-order Padé approximation of pure delay e^{−sT} (+1 more overload) |
| [`pade2`](inc/damp/matlab.hpp#L1020) | Second-order Padé approximation of pure delay e^{−sT} (+1 more overload) |
| [`pid`](inc/damp/matlab.hpp#L161) | MATLAB®-style parallel-form continuous PID constructor |
| [`pidstd`](inc/damp/matlab.hpp#L189) | Standard-form continuous PID constructor (1-DOF) |
| [`pidstd2`](inc/damp/matlab.hpp#L207) | Standard-form continuous 2-DOF PID constructor |
| [`pidtune`](inc/damp/matlab.hpp#L787) | PID controller tuning using frequency domain method |
| [`pinv`](inc/damp/matlab.hpp#L446) | MATLAB® short alias for the Moore–Penrose pseudoinverse |
| [`place`](inc/damp/matlab.hpp#L582) | Robust multi-input pole placement (MATLAB®'s place) |
| [`pole`](inc/damp/matlab.hpp#L853) | MATLAB® short alias for the open-loop poles of a system |
| [`quadprog`](inc/damp/matlab.hpp#L1073) | MATLAB® alias for the dense inequality-constrained QP solve (+1 more overload) |
| [`reg`](inc/damp/matlab.hpp#L519) | Form dynamic regulator from system, state-feedback gain, and estimator gain |
| [`ss`](inc/damp/matlab.hpp#L97) | MATLAB®-style state-space model constructor |
| [`svd`](inc/damp/matlab.hpp#L437) | MATLAB® short alias for the singular value decomposition |
| [`tf`](inc/damp/matlab.hpp#L61) | MATLAB®-style transfer function constructor (+1 more overload) |
| [`zpk`](inc/damp/matlab.hpp#L120) | MATLAB®-style zero-pole-gain model constructor |

## Math backends

Internal, compile-time-selected implementations of the `damp::` scalar-math surface (`sin`/`cos`/`sqrt`/`exp`/…), chosen via `damp/config.hpp`. Every backend exposes the same functions, so they're listed once here as files rather than repeated in the tables above. The public dispatcher is [`damp/math/math.hpp`](inc/damp/math/math.hpp).

| File | Role |
| ---- | ---- |
| [`damp/math/math_backend.hpp`](inc/damp/math/math_backend.hpp) | Pluggable runtime math backend selection (freestanding-safe) |
| [`damp/math/std_fallback.hpp`](inc/damp/math/std_fallback.hpp) | Composable std:: (\<cmath\>) base for the hosted math backends |
| [`damp/math/damp_backend.hpp`](inc/damp/math/damp_backend.hpp) | Fast float/double math backend (trig.hpp) over the std:: fallback |
| [`damp/math/series_backend.hpp`](inc/damp/math/series_backend.hpp) | Freestanding runtime math backend: routes scalar math to the constexpr series in constexpr_math.hpp, pulling no hosted headers (no \<cmath\>) |
| [`damp/math/constexpr_math.hpp`](inc/damp/math/constexpr_math.hpp) | Compile-time scalar math for the damp:: dispatch layer |
| [`damp/math/trig.hpp`](inc/damp/math/trig.hpp) | Fast float/double sine, cosine, and inverse trig with Cody–Waite reduction |

## Examples

Runnable programs in `examples/` (104 total). Build with `make` (or `tup --quiet examples`); outputs go to `examples/build/`.

| Example | Description |
| ------- | ----------- |
| [`acim_sil.cpp`](examples/motor/acim/acim_sil.cpp) | ACIM IFOC cascade — host smoke (calls acim_controller.hpp) |
| [`acim_sketch.cpp`](examples/motor/acim/acim_sketch.cpp) | ACIM IFOC cascade servo — flashable sketch |
| [`anpc_sil.cpp`](examples/power/anpc/anpc_sil.cpp) | anpc — host hybrid SIL (calls anpc_controller.hpp) |
| [`anpc_sketch.cpp`](examples/power/anpc/anpc_sketch.cpp) | anpc — flashable sketch (Design Is Deploy) |
| [`boost_pfc_sil.cpp`](examples/power/boost_pfc/boost_pfc_sil.cpp) | boost_pfc — host hybrid SIL (calls boost_pfc_controller.hpp) |
| [`boost_pfc_sketch.cpp`](examples/power/boost_pfc/boost_pfc_sketch.cpp) | boost_pfc — flashable sketch (Design Is Deploy) |
| [`boost_sil.cpp`](examples/power/boost/boost_sil.cpp) | boost — host hybrid SIL (calls boost_controller.hpp) |
| [`boost_sketch.cpp`](examples/power/boost/boost_sketch.cpp) | boost — flashable sketch (Design Is Deploy) |
| [`buck_boost_sil.cpp`](examples/power/buck_boost/buck_boost_sil.cpp) | buck_boost — host hybrid SIL (calls buck_boost_controller.hpp) |
| [`buck_boost_sketch.cpp`](examples/power/buck_boost/buck_boost_sketch.cpp) | buck_boost — flashable sketch (Design Is Deploy) |
| [`buck_sil.cpp`](examples/power/buck/buck_sil.cpp) | buck — host hybrid SIL (calls buck_controller.hpp) |
| [`buck_sketch.cpp`](examples/power/buck/buck_sketch.cpp) | buck — flashable sketch (Design Is Deploy) |
| [`camming_gallery_sil.cpp`](examples/motion/camming_gallery/camming_gallery_sil.cpp) | Electronic camming gallery host SIL (calls camming_gallery_deploy.hpp) |
| [`camming_gallery_sketch.cpp`](examples/motion/camming_gallery/camming_gallery_sketch.cpp) | Camming gallery — flashable rise-law smoke |
| [`cart_pole_sil.cpp`](examples/control/cart_pole/cart_pole_sil.cpp) | Cart-pole LQR — host nonlinear SIL (calls cart_pole_controller.hpp) |
| [`cart_pole_sketch.cpp`](examples/control/cart_pole/cart_pole_sketch.cpp) | Cart-pole LQR — flashable sketch (Design Is Deploy) |
| [`coordinated_joints_sil.cpp`](examples/motion/coordinated_joints/coordinated_joints_sil.cpp) | Coordinated joints host SIL (calls coordinated_joints_deploy.hpp) |
| [`coordinated_joints_sketch.cpp`](examples/motion/coordinated_joints/coordinated_joints_sketch.cpp) | Coordinated joint move — flashable TrajectoryBank smoke |
| [`corexy_sil.cpp`](examples/motion/corexy/corexy_sil.cpp) | CoreXY time-optimal move host SIL (calls corexy_deploy.hpp) |
| [`corexy_sketch.cpp`](examples/motion/corexy/corexy_sketch.cpp) | CoreXY move — flashable command stream smoke |
| [`cuk_sil.cpp`](examples/power/cuk/cuk_sil.cpp) | cuk — host hybrid SIL (calls cuk_controller.hpp) |
| [`cuk_sketch.cpp`](examples/power/cuk/cuk_sketch.cpp) | cuk — flashable sketch (Design Is Deploy) |
| [`drive_governor_sil.cpp`](examples/motor/drive_governor/drive_governor_sil.cpp) | Drive selector governor — SIL plots + optional flight-stick HID / live UI |
| [`drive_governor_sketch.cpp`](examples/motor/drive_governor/drive_governor_sketch.cpp) | Drive selector governor — thin sketch (mock sticks) |
| [`encoder_velocity_sil.cpp`](examples/estimation/encoder_velocity/encoder_velocity_sil.cpp) | Encoder / tach velocity estimation — host comparison + plots |
| [`encoder_velocity_sketch.cpp`](examples/estimation/encoder_velocity/encoder_velocity_sketch.cpp) | Encoder PLL velocity — thin deploy-shaped smoke |
| [`eskf_sil.cpp`](examples/estimation/eskf/eskf_sil.cpp) | MARG attitude ESKF — finite-tick host smoke (calls eskf_estimator.hpp) |
| [`eskf_sketch.cpp`](examples/estimation/eskf/eskf_sketch.cpp) | MARG attitude ESKF — flashable sketch (Design Is Deploy) |
| [`example_math_backend.cpp`](examples/example_math_backend.cpp) | Pluggable math backend example |
| [`excavator_sil.cpp`](examples/motion/excavator/excavator_sil.cpp) | Excavator dig SIL + optional live cylinder teleop (calls excavator_estimator.hpp) |
| [`excavator_sketch.cpp`](examples/motion/excavator/excavator_sketch.cpp) | Excavator pose filter — flashable sketch (Design Is Deploy) |
| [`fo_plant_inertia_sil.cpp`](examples/estimation/fo_plant_inertia/fo_plant_inertia_sil.cpp) | Online J,b from FirstOrderPlantEstimator — host identification demo |
| [`fo_plant_inertia_sketch.cpp`](examples/estimation/fo_plant_inertia/fo_plant_inertia_sketch.cpp) | FO plant inertia ID — thin deploy-shaped smoke |
| [`foc_sil.cpp`](examples/motor/foc/foc_sil.cpp) | Host teaching plot: FOC current-loop PI vs I-P (calls foc_controller.hpp) |
| [`foc_sketch.cpp`](examples/motor/foc/foc_sketch.cpp) | FOC current-loop — thin smoke sketch (host teaching demo) |
| [`foc_switching_sil.cpp`](examples/motor/foc_switching/foc_switching_sil.cpp) | Switching FOC PMSM — host Exact αβ plant SIL (calls foc_switching_controller.hpp) |
| [`foc_switching_sketch.cpp`](examples/motor/foc_switching/foc_switching_sketch.cpp) | Switching FOC PMSM — flashable sketch (Design Is Deploy) |
| [`four_switch_buck_boost_sil.cpp`](examples/power/four_switch_buck_boost/four_switch_buck_boost_sil.cpp) | four_switch_buck_boost — host hybrid SIL (calls four_switch_buck_boost_controller.hpp) |
| [`four_switch_buck_boost_sketch.cpp`](examples/power/four_switch_buck_boost/four_switch_buck_boost_sketch.cpp) | four_switch_buck_boost — flashable sketch (Design Is Deploy) |
| [`grid_3ph_sil.cpp`](examples/power/grid_3ph/grid_3ph_sil.cpp) | grid_3ph — host hybrid SIL (calls grid_3ph_controller.hpp) |
| [`grid_3ph_sketch.cpp`](examples/power/grid_3ph/grid_3ph_sketch.cpp) | grid_3ph — flashable sketch (Design Is Deploy) |
| [`imu_pose_sil.cpp`](examples/estimation/imu_pose/imu_pose_sil.cpp) | Animated body-frame triad from gyro integration (IMU attitude only) |
| [`imu_pose_sketch.cpp`](examples/estimation/imu_pose/imu_pose_sketch.cpp) | Open-loop gyro attitude — thin smoke (teaching) |
| [`input_shaper_sil.cpp`](examples/control/input_shaper/input_shaper_sil.cpp) | Input-shaper calibration from accelerometer (Klipper-style host demo) |
| [`input_shaper_sketch.cpp`](examples/control/input_shaper/input_shaper_sketch.cpp) | Input shaper — flashable sketch (fixed ZVD at nameplate fn, ζ) |
| [`ins_eskf_sil.cpp`](examples/estimation/ins_eskf/ins_eskf_sil.cpp) | Animated INS: free-run bias vs InsNavigator (position + dual-antenna heading) |
| [`ins_eskf_sketch.cpp`](examples/estimation/ins_eskf/ins_eskf_sketch.cpp) | InsNavigator thin smoke — same design as ins_eskf_sil (float deploy) |
| [`ins_mechanization_sil.cpp`](examples/estimation/ins_mechanization/ins_mechanization_sil.cpp) | Animated strapdown INS path: IMU samples → p, v, q (mechanization) |
| [`ins_mechanization_sketch.cpp`](examples/estimation/ins_mechanization/ins_mechanization_sketch.cpp) | Strapdown mechanization — thin rest-frame smoke |
| [`ins_navigator_sil.cpp`](examples/estimation/ins_navigator/ins_navigator_sil.cpp) | 15-state INS navigator — finite-tick host smoke (calls ins_navigator_estimator.hpp) |
| [`ins_navigator_sketch.cpp`](examples/estimation/ins_navigator/ins_navigator_sketch.cpp) | 15-state INS navigator — flashable sketch (Design Is Deploy) |
| [`ipm_sil.cpp`](examples/motor/ipm/ipm_sil.cpp) | IPM cascade servo — host smoke (calls ipm_controller.hpp) |
| [`ipm_sketch.cpp`](examples/motor/ipm/ipm_sketch.cpp) | IPM cascade servo — flashable sketch (Design Is Deploy) |
| [`kinematic_maps_sil.cpp`](examples/motion/kinematic_maps/kinematic_maps_sil.cpp) | Kinematic maps host SIL (calls kinematic_maps_deploy.hpp) |
| [`kinematic_maps_sketch.cpp`](examples/motion/kinematic_maps/kinematic_maps_sketch.cpp) | Kinematic maps — flashable inverse-map smoke |
| [`lpf_sil.cpp`](examples/control/lpf/lpf_sil.cpp) | Low-pass filter — host step-response demo (calls lpf_filter.hpp) |
| [`lpf_sketch.cpp`](examples/control/lpf/lpf_sketch.cpp) | Low-pass filter — flashable sketch |
| [`multirate_rig_sil.cpp`](examples/control/multirate_rig/multirate_rig_sil.cpp) | Multi-rate pose → joint cascade — host SIL (calls multirate_rig_controller.hpp) |
| [`multirate_rig_sketch.cpp`](examples/control/multirate_rig/multirate_rig_sketch.cpp) | Multi-rate pose → joint cascade — flashable sketch |
| [`npc_sil.cpp`](examples/power/npc/npc_sil.cpp) | npc — host hybrid SIL (calls npc_controller.hpp) |
| [`npc_sketch.cpp`](examples/power/npc/npc_sketch.cpp) | npc — flashable sketch (Design Is Deploy) |
| [`pendulum_sil.cpp`](examples/control/pendulum/pendulum_sil.cpp) | Upright pendulum LQR — host nonlinear SIL (calls pendulum_controller.hpp) |
| [`pendulum_sketch.cpp`](examples/control/pendulum/pendulum_sketch.cpp) | Upright pendulum LQR — flashable sketch (Design Is Deploy) |
| [`pid_sil.cpp`](examples/control/pid/pid_sil.cpp) | PI — thin host smoke (calls pid_controller.hpp) |
| [`pid_sketch.cpp`](examples/control/pid/pid_sketch.cpp) | PI — flashable sketch (Design Is Deploy) |
| [`pmac_bus_limit_sil.cpp`](examples/PMAC/pmac_bus_limit/pmac_bus_limit_sil.cpp) | pmac_bus_limit host SIL (calls pmac_bus_limit_deploy.hpp) |
| [`pmac_bus_limit_sketch.cpp`](examples/PMAC/pmac_bus_limit/pmac_bus_limit_sketch.cpp) | DC-bus limiter — flashable evaluate smoke |
| [`pmac_calibration_sil.cpp`](examples/PMAC/pmac_calibration/pmac_calibration_sil.cpp) | pmac_calibration host SIL (calls pmac_calibration_deploy.hpp) |
| [`pmac_calibration_sketch.cpp`](examples/PMAC/pmac_calibration/pmac_calibration_sketch.cpp) | Phase calibrator — flashable PRBS step smoke |
| [`pmac_estimator_sil.cpp`](examples/PMAC/pmac_estimator/pmac_estimator_sil.cpp) | pmac_estimator host SIL (calls pmac_estimator_estimator.hpp) |
| [`pmac_estimator_sketch.cpp`](examples/PMAC/pmac_estimator/pmac_estimator_sketch.cpp) | Mechanical estimator — flashable predict/update smoke |
| [`pmac_servo_sil.cpp`](examples/PMAC/pmac_servo/pmac_servo_sil.cpp) | pmac_servo host SIL (calls pmac_servo_deploy.hpp) |
| [`pmac_servo_sketch.cpp`](examples/PMAC/pmac_servo/pmac_servo_sketch.cpp) | PMAC servo — flashable make_pmsm_servo smoke |
| [`pmac_thermal_sil.cpp`](examples/PMAC/pmac_thermal/pmac_thermal_sil.cpp) | pmac_thermal host SIL (calls pmac_thermal_deploy.hpp) |
| [`pmac_thermal_sketch.cpp`](examples/PMAC/pmac_thermal/pmac_thermal_sketch.cpp) | Junction thermal — flashable Foster step smoke |
| [`pmsm_sil.cpp`](examples/motor/pmsm/pmsm_sil.cpp) | PMSM cascade servo — host smoke (calls pmsm_controller.hpp) |
| [`pmsm_sketch.cpp`](examples/motor/pmsm/pmsm_sketch.cpp) | PMSM cascade servo — flashable sketch (Design Is Deploy) |
| [`reaction_wheel_sil.cpp`](examples/control/reaction_wheel/reaction_wheel_sil.cpp) | Reaction-wheel gallery — print gains + one tick each (host) |
| [`reaction_wheel_sketch.cpp`](examples/control/reaction_wheel/reaction_wheel_sketch.cpp) | Reaction-wheel cascade PI — flashable sketch (Design Is Deploy) |
| [`sensorless_ipm_sil.cpp`](examples/motor/sensorless_ipm/sensorless_ipm_sil.cpp) | IPM sensorless HFI hybrid — host SIL + plots (calls sensorless_ipm_controller.hpp) |
| [`sensorless_ipm_sketch.cpp`](examples/motor/sensorless_ipm/sensorless_ipm_sketch.cpp) | IPM sensorless hybrid — flashable sketch |
| [`sepex_sil.cpp`](examples/motor/sepex/sepex_sil.cpp) | SEPEX cascade — host smoke (calls sepex_controller.hpp) |
| [`sepex_sketch.cpp`](examples/motor/sepex/sepex_sketch.cpp) | SEPEX cascade servo — flashable sketch |
| [`sepic_sil.cpp`](examples/power/sepic/sepic_sil.cpp) | sepic — host hybrid SIL (calls sepic_controller.hpp) |
| [`sepic_sketch.cpp`](examples/power/sepic/sepic_sketch.cpp) | sepic — flashable sketch (Design Is Deploy) |
| [`serial_arm_sil.cpp`](examples/motion/serial_arm/serial_arm_sil.cpp) | 3R serial-arm host FK scrub (calls serial_arm_deploy.hpp) |
| [`serial_arm_sketch.cpp`](examples/motion/serial_arm/serial_arm_sketch.cpp) | 3R serial arm — flashable FK smoke (Design Is Deploy) |
| [`servo_commands_sil.cpp`](examples/motor/servo_commands/servo_commands_sil.cpp) | Trajectory → servo drive command stream (calls servo_commands_controller.hpp) |
| [`servo_commands_sketch.cpp`](examples/motor/servo_commands/servo_commands_sketch.cpp) | Trajectory → servo commands — thin stream sketch |
| [`servo_drive_sil.cpp`](examples/motor/servo_drive/servo_drive_sil.cpp) | Two-mass PMSM SIL with Cascade + FOController (calls servo_drive_controller.hpp) |
| [`servo_drive_sketch.cpp`](examples/motor/servo_drive/servo_drive_sketch.cpp) | Two-mass PMSM cascade + FOC — thin sketch smoke |
| [`servo_sim.cpp`](examples/servo_drive/servo_sim.cpp) | Servo sim |
| [`swashplate_sil.cpp`](examples/motion/swashplate/swashplate_sil.cpp) | Swash-plate STSMC host comparison SIL (calls swashplate_controller.hpp) |
| [`swashplate_sketch.cpp`](examples/motion/swashplate/swashplate_sketch.cpp) | Swash-plate STSMC — flashable sketch |
| [`totem_pole_pfc_sil.cpp`](examples/power/totem_pole_pfc/totem_pole_pfc_sil.cpp) | totem_pole_pfc — host hybrid SIL (calls totem_pole_pfc_controller.hpp) |
| [`totem_pole_pfc_sketch.cpp`](examples/power/totem_pole_pfc/totem_pole_pfc_sketch.cpp) | totem_pole_pfc — flashable sketch (Design Is Deploy) |
| [`trajectory_gallery_sil.cpp`](examples/motion/trajectory_gallery/trajectory_gallery_sil.cpp) | Motion-profile gallery host SIL (calls trajectory_gallery_deploy.hpp) |
| [`trajectory_gallery_sketch.cpp`](examples/motion/trajectory_gallery/trajectory_gallery_sketch.cpp) | Trajectory gallery — flashable S-curve smoke on first scenario |
| [`ttype_sil.cpp`](examples/power/ttype/ttype_sil.cpp) | ttype — host hybrid SIL (calls ttype_controller.hpp) |
| [`ttype_sketch.cpp`](examples/power/ttype/ttype_sketch.cpp) | ttype — flashable sketch (Design Is Deploy) |
| [`workflow_end_to_end_sil.cpp`](examples/control/workflow_end_to_end/workflow_end_to_end_sil.cpp) | Workflow end-to-end host SIL (calls workflow_end_to_end_controller.hpp) |
| [`workflow_end_to_end_sketch.cpp`](examples/control/workflow_end_to_end/workflow_end_to_end_sketch.cpp) | Workflow e2e — flashable float runtime smoke (after host design) |
| [`zeta_sil.cpp`](examples/power/zeta/zeta_sil.cpp) | zeta — host hybrid SIL (calls zeta_controller.hpp) |
| [`zeta_sketch.cpp`](examples/power/zeta/zeta_sketch.cpp) | zeta — flashable sketch (Design Is Deploy) |
