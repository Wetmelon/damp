# API Reference — Alphabetical Index

Auto-generated from `@brief` doc comments in `inc/damp/`. Regenerate with `python tools/gen_reference.py`. Grouped-by-domain view: [REFERENCE.md](REFERENCE.md).

| Name | Kind | Domain | Description |
| ---- | ---- | ------ | ----------- |
| [`abs`](inc/damp/math/complex.hpp#L358) | function | Scalar math & complex | Compute magnitude (absolute value) of a complex number (+1 more overload) |
| [`Acim`](inc/damp/motor/acim.hpp#L282) | block | Motor control pack (if present) | AC induction machine (ACIM) electrical stage — indirect FOC |
| [`AcimConfig`](inc/damp/motor/acim.hpp#L254) | block | Motor control pack (if present) | Configuration for an Acim electrical stage |
| [`AcimServoConfig`](inc/damp/motor/acim.hpp#L440) | block | Motor control pack (if present) | Flat configuration for an ACIM cascade servo (make_acim_servo) |
| [`acker`](inc/damp/matlab.hpp#L554) | function | MATLAB®-style aliases (host) | Pole placement for state-feedback control |
| [`ackermann`](inc/damp/design/pole_placement.hpp#L1582) | function | Design-time synthesis (not PWM-rate) | Single-input pole placement via Ackermann's formula |
| [`acos`](inc/damp/math/math.hpp#L150) | function | Scalar math & complex | Arccosine ∈ [0, π]. Input is clamped to [−1, 1] in both paths |
| [`AcPolarity`](inc/damp/power/totem_pole.hpp#L45) | enum | Power electronics pack (if present) | AC line polarity for totem-pole slow (LF) leg |
| [`ActionGovernor`](inc/damp/controllers/action_governor.hpp#L331) | block | Runtime controllers | Runtime action governor — QP projection of u_des onto A u ≤ b |
| [`ActiveSetSolver`](inc/damp/design/qp.hpp#L618) | block | Design-time synthesis (not PWM-rate) | Default QP solver policy: the Goldfarb–Idnani active-set solve_qp() |
| [`AdaptationGate`](inc/damp/estimation/parameter_estimation.hpp#L315) | block | Observers & estimators | Gated-adaptation policy: decide whether a parameter update is safe to apply |
| [`AdaptationGateDecision`](inc/damp/analysis/identification.hpp#L226) | block | Frequency-domain analysis (host) | Gate decision for whether adaptive updates are currently allowed |
| [`adaptive_solve`](inc/damp/simulation/solver.hpp#L774) | function | Simulation / SIL harness (host) | Adaptive-step solve in one call |
| [`AdaptiveLut1D`](inc/damp/toolbox/lookup.hpp#L427) | block | Embedded helpers (controls-adjacent utilities) | Adaptive 1-D lookup: fixed breakpoints, online-updated cell values |
| [`AdaptiveLut2D`](inc/damp/toolbox/lookup.hpp#L502) | block | Embedded helpers (controls-adjacent utilities) | Adaptive 2-D lookup: fixed grid, online-updated cells (nearest grid node) |
| [`AdaptiveOptions`](inc/damp/simulation/solver.hpp#L306) | block | Simulation / SIL harness (host) | Options for adaptive-step ODE integration |
| [`AdaptiveStepIntegrator`](inc/damp/simulation/integrator.hpp#L65) | concept | Simulation / SIL harness (host) | Integrator that reports a genuine embedded local-error estimate |
| [`AdaptiveStepSolver`](inc/damp/simulation/solver.hpp#L358) | block | Simulation / SIL harness (host) | Adaptive-step ODE solver (+1 more overload) |
| [`AdaptOutOfRange`](inc/damp/toolbox/lookup.hpp#L48) | enum | Embedded helpers (controls-adjacent utilities) | How AdaptiveLut1D / AdaptiveLut2D treat samples outside the breakpoint span |
| [`AdmmSettings`](inc/damp/design/qp.hpp#L713) | block | Design-time synthesis (not PWM-rate) | Tuning parameters for the ADMM QP solver |
| [`AdmmSolver`](inc/damp/design/qp.hpp#L752) | block | Design-time synthesis (not PWM-rate) | ADMM (OSQP-style) QP solver policy — warm-started, fixed-cost iterations |
| [`adrc`](inc/damp/controllers/adrc.hpp#L86) | function | Design-time synthesis (not PWM-rate) | Active Disturbance Rejection Control design |
| [`ADRCController`](inc/damp/controllers/adrc.hpp#L148) | block | Runtime controllers | Active Disturbance Rejection Control (ADRC) |
| [`ADRCResult`](inc/damp/controllers/adrc.hpp#L33) | block | Design-time synthesis (not PWM-rate) | Active Disturbance Rejection Control design result |
| [`AffineCal`](inc/damp/toolbox/scaling.hpp#L83) | block | Embedded helpers (controls-adjacent utilities) | Affine sensor calibration `y = gain·x + offset` |
| [`allmargin`](inc/damp/matlab.hpp#L937) | function | MATLAB®-style aliases (host) | Gain, phase, and delay margins of a SISO loop over a frequency grid |
| [`AllMarginResult`](inc/damp/matlab.hpp#L916) | block | MATLAB®-style aliases (host) | All classical margins including delay margin (superset of margin) |
| [`AlphaBeta`](inc/damp/transforms.hpp#L152) | block | Motor control pack (if present) | Alpha-beta (stationary-frame) component pair |
| [`AlphaBetaZero`](inc/damp/transforms.hpp#L218) | block | Motor control pack (if present) | Alpha-beta-zero (stationary-frame) component triple |
| [`amigo_kappa_tau`](inc/damp/design/pid_design.hpp#L237) | function | Design-time synthesis (not PWM-rate) | AMIGO PI from ultimate gain/period and static gain Kₛ |
| [`AnalogCrossMode`](inc/damp/toolbox/conditioning.hpp#L498) | enum | Embedded helpers (controls-adjacent utilities) | How two redundant analog channels should relate in engineering units |
| [`AnalogCrossResult`](inc/damp/toolbox/conditioning.hpp#L520) | block | Embedded helpers (controls-adjacent utilities) | Output of cross_check_analog |
| [`AnalogCrossStatus`](inc/damp/toolbox/conditioning.hpp#L506) | enum | Embedded helpers (controls-adjacent utilities) | Result of comparing two redundant analog channels |
| [`AnalogInput`](inc/damp/toolbox/io.hpp#L46) | block | Embedded helpers (controls-adjacent utilities) | A single analog input: range/fault check on the raw reading, then affine calibration to engineering units |
| [`angle_delay_compensation`](inc/damp/motor/scalar_control.hpp#L224) | function | Motor control pack (if present) | Compensate an electrical angle for the control/PWM transport delay |
| [`AngleGenerator`](inc/damp/motor/scalar_control.hpp#L54) | block | Motor control pack (if present) | Fixed-frequency electrical-angle generator with a slew-limited frequency |
| [`AngleUnwrapper`](inc/damp/toolbox/encoder.hpp#L190) | block | Embedded helpers (controls-adjacent utilities) | Wrap-safe accumulator: successive wrapped absolute readings → continuous position |
| [`AnpcMap`](inc/damp/power/anpc.hpp#L37) | block | Power electronics pack (if present) | Single-phase ANPC leg map — voltage command → 6 duties (layer 3) |
| [`apply_voltage_current_limits`](inc/damp/motor/ipm.hpp#L94) | function | Motor control pack (if present) | Apply analytic voltage FW then current-circle clamp (open-loop, no integrator) |
| [`apply_wet_plotly_shell`](inc/damp/simulation/plot_plotly.hpp#L64) | function | Simulation / SIL harness (host) | Apply damp's full-window HTML shell to plotlypp-generated markup |
| [`arange`](inc/damp/analysis/linspace.hpp#L134) | function | Frequency-domain analysis (host) | Half-open arithmetic range [start, stop) with step step (+2 more overloads) |
| [`arcade_drive`](inc/damp/toolbox/io.hpp#L361) | function | Embedded helpers (controls-adjacent utilities) | Arcade (single-stick) drive mixer: x (steer) + y (throttle) → left/right |
| [`arg`](inc/damp/math/complex.hpp#L369) | function | Scalar math & complex | Compute argument (phase angle) of a complex number |
| [`arm_frames_world`](inc/damp/estimation/serial_arm_pose.hpp#L848) | function | Observers & estimators | World-frame DH poses: fixed chain origin, base body→world R0 |
| [`arm_spherical_wrist`](inc/damp/kinematics/serial_arm.hpp#L601) | function | Kinematics / pose | Tier-2 builder for a standard 6R elbow arm with a spherical wrist |
| [`ArmIkResult`](inc/damp/kinematics/serial_arm.hpp#L164) | block | Kinematics / pose | Result of a numerical inverse-kinematics solve |
| [`ARXModel`](inc/damp/analysis/identification.hpp#L103) | block | Frequency-domain analysis (host) | ARX candidate model with fixed numerator/denominator orders |
| [`asin`](inc/damp/math/math.hpp#L131) | function | Scalar math & complex | Arcsine ∈ [−π/2, π/2]. Input is clamped to [−1, 1] in both paths so behavior matches at compile and run time (std::asin would return NaN for \|x\| > 1) |
| [`atan`](inc/damp/math/math.hpp#L117) | function | Scalar math & complex | Single-argument arctangent ∈ (−π/2, π/2) |
| [`atan2`](inc/damp/math/math.hpp#L105) | function | Scalar math & complex | Two-argument arctangent, atan2(y, x) ∈ [−π, π] |
| [`axis_snap_quantum`](inc/damp/trajectory/online_otg.hpp#L220) | function | Trajectory value types | Near-target snap quantum (~0.01% of scale; fixed-point-style residual floor) |
| [`AxisInput`](inc/damp/toolbox/io.hpp#L191) | block | Embedded helpers (controls-adjacent utilities) | Operator-axis conditioning chain (joystick / RC stick → command) |
| [`AxisMotionLimits`](inc/damp/trajectory/online_otg.hpp#L64) | block | Trajectory value types | Per-axis value / asymmetric accel–decel / jerk caps |
| [`AxisState`](inc/damp/trajectory/online_otg.hpp#L80) | block | Trajectory value types | Shared kinematic state for one axis: value and its derivative |
| [`back_emf`](inc/damp/motor/sepex.hpp#L81) | function | Motor control pack (if present) | Armature back-EMF e_a = K_af i_f ω [V] (omega in rad/s) |
| [`backward_substitute_transpose`](inc/damp/matrix/solve.hpp#L59) | function | Linear algebra | Backward substitution for Lᵀx = b (real) / Lᴴx = b layout |
| [`BackwardEuler`](inc/damp/simulation/integrator.hpp#L295) | block | Simulation / SIL harness (host) | Backward Euler integrator |
| [`bake_torque_limit_lut`](inc/damp/motor/ipm.hpp#L235) | function | Motor control pack (if present) | Bake a 1-D torque ceiling vs mechanical speed into a Lut1D (float-friendly) |
| [`balanced_realization`](inc/damp/design/model_reduction.hpp#L569) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for balreal |
| [`balanced_truncation`](inc/damp/design/model_reduction.hpp#L586) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for balred |
| [`BalancedRealizationResult`](inc/damp/design/model_reduction.hpp#L74) | block | Design-time synthesis (not PWM-rate) | Balanced realization result (Moore square-root method) |
| [`BalancedReductionResult`](inc/damp/design/model_reduction.hpp#L100) | block | Design-time synthesis (not PWM-rate) | Reduced-order model from balanced truncation or residualization |
| [`balreal`](inc/damp/design/model_reduction.hpp#L383) | function | Design-time synthesis (not PWM-rate) | Balanced realization via the square-root (Moore / Laub) method |
| [`balred`](inc/damp/design/model_reduction.hpp#L471) | function | Design-time synthesis (not PWM-rate) | Balanced truncation to NR states |
| [`bandpass`](inc/damp/filters/iir_design.hpp#L661) | function | Filters & signal conditioning | Second-order band-pass filter (constant 0 dB peak gain) |
| [`bandwidth`](inc/damp/matlab.hpp#L1042) | function | MATLAB®-style aliases (host) | -3 dB bandwidth of a SISO system over a frequency grid |
| [`bandwidth_from_settling_time`](inc/damp/design/pid_design.hpp#L571) | function | Design-time synthesis (not PWM-rate) | Map settling-time and damping-ratio targets to a bandwidth estimate |
| [`base_speed`](inc/damp/motor/foc.hpp#L294) | function | Motor control pack (if present) | Base (corner) electrical speed where the voltage circle is first hit (+1 more overload) |
| [`base_speed_elec`](inc/damp/motor/acim.hpp#L179) | function | Motor control pack (if present) | Electrical base speed [elec rad/s] at rated flux and voltage ceiling |
| [`BDF2`](inc/damp/simulation/integrator.hpp#L366) | block | Simulation / SIL harness (host) | Backward Differentiation Formula 2 (BDF2) integrator |
| [`beta`](inc/damp/toolbox/thermistor.hpp#L77) | function | Embedded helpers (controls-adjacent utilities) | Fit NTC coefficients from the Beta-parameter model |
| [`Biquad`](inc/damp/filters/biquad.hpp#L34) | block | Filters & signal conditioning | Second-order IIR (biquad) section runtime |
| [`BiquadCascade`](inc/damp/filters/biquad.hpp#L101) | block | Filters & signal conditioning | Cascade of second-order sections (SOS) for higher-order IIR filters |
| [`BLINK`](inc/damp/toolbox/iec61131.hpp#L550) | block | Embedded helpers (controls-adjacent utilities) | BLINK (free-running square-wave / flasher) |
| [`blkdiag`](inc/damp/matlab.hpp#L349) | function | MATLAB®-style aliases (host) | Block diagonal matrix construction |
| [`Block`](inc/damp/matrix/block.hpp#L38) | block | Linear algebra | Block view (non-owning) into a parent matrix |
| [`bode`](inc/damp/analysis/frequency.hpp#L250) | function | Frequency-domain analysis (host) | Compute Bode plot data for a SISO state-space system (+2 more overloads) |
| [`bode_discrete`](inc/damp/analysis/frequency.hpp#L358) | function | Frequency-domain analysis (host) | Compute Bode plot data for a discrete-time SISO state-space system |
| [`bodemag`](inc/damp/simulation/plot_plotly.hpp#L796) | function | Simulation / SIL harness (host) | Plot a magnitude-only Bode diagram (log frequency, dB magnitude) |
| [`bodeplot`](inc/damp/simulation/plot_plotly.hpp#L782) | function | Simulation / SIL harness (host) | Plot magnitude and phase Bode subplots |
| [`BodeResult`](inc/damp/analysis/frequency.hpp#L60) | block | Frequency-domain analysis (host) | Bode plot data for a SISO system |
| [`BoostMap`](inc/damp/power/boost.hpp#L42) | block | Power electronics pack (if present) | Boost CCM plant map — voltage command → duty (layer 3) |
| [`Bounds`](inc/damp/toolbox/bounds.hpp#L43) | block | Embedded helpers (controls-adjacent utilities) | A per-channel closed-interval box constraint |
| [`BoxCommandFilter`](inc/damp/controllers/action_governor.hpp#L119) | block | Runtime controllers | Runtime box command filter — closed-form clamp per tick |
| [`BrushedDc`](inc/damp/motor/brushed_dc.hpp#L66) | block | Motor control pack (if present) | Brushed-DC electrical stage: single armature current loop, one H-bridge duty |
| [`BrushedDcConfig`](inc/damp/motor/brushed_dc.hpp#L44) | block | Motor control pack (if present) | Configuration for a BrushedDc electrical stage |
| [`BrushedDcServoConfig`](inc/damp/motor/brushed_dc.hpp#L144) | block | Motor control pack (if present) | Flat configuration for a brushed-DC cascade servo (make_brushed_servo) |
| [`BuckBoostMap`](inc/damp/power/buck_boost.hpp#L31) | block | Power electronics pack (if present) | Inverting buck-boost CCM plant map — HardSwitchingConverter (layer 3) |
| [`BuckMap`](inc/damp/power/buck.hpp#L33) | block | Power electronics pack (if present) | Buck CCM plant map — voltage command → duty (layer 3) |
| [`build_lqg_analysis_models`](inc/damp/design/synthesis.hpp#L211) | function | Design-time synthesis (not PWM-rate) | Build analysis models from an LQG design |
| [`build_lqgi_analysis_models`](inc/damp/design/synthesis.hpp#L290) | function | Design-time synthesis (not PWM-rate) | Build analysis models from an LQGI servo design |
| [`build_lqi_analysis_models`](inc/damp/design/synthesis.hpp#L259) | function | Design-time synthesis (not PWM-rate) | Build analysis models from an LQI servo design |
| [`build_mpc_analysis_models`](inc/damp/controllers/mpc.hpp#L690) | function | Design-time synthesis (not PWM-rate) | Build the unconstrained-MPC LTI analysis models |
| [`butterworth_lowpass`](inc/damp/filters/iir_design.hpp#L314) | function | Filters & signal conditioning | Butterworth low-pass filter design (+1 more overload) |
| [`Button`](inc/damp/toolbox/io.hpp#L238) | block | Embedded helpers (controls-adjacent utilities) | Debounced momentary push-button with edge and long-press detection |
| [`c2d`](inc/damp/matlab.hpp#L296) | function | MATLAB®-style aliases (host) | MATLAB® interface function c2d to discretize a continuous-time state-space system (+1 more overload) |
| [`CachedZoh`](inc/damp/simulation/cached_zoh.hpp#L62) | block | Simulation / SIL harness (host) | Cached discrete ZOH maps Ad, Bd for step size h |
| [`cam_profile`](inc/damp/trajectory/cam.hpp#L136) | function | Trajectory value types | Build a CamProfile from a synthesized spline and cycle length |
| [`CamFollower`](inc/damp/trajectory/cam.hpp#L481) | block | Trajectory value types | Runtime follower of an electronic cam table (design::CamProfile) |
| [`CamProfile`](inc/damp/trajectory/cam.hpp#L99) | block | Trajectory value types | Cam table: follower position as a function of leader parameter over one cycle |
| [`canonical_phase_margin`](inc/damp/analysis/frequency.hpp#L148) | function | Frequency-domain analysis (host) | Normalize phase margin to (-180, 180] |
| [`care`](inc/damp/design/riccati.hpp#L763) | function | Design-time synthesis (not PWM-rate) | Solve the Continuous-time Algebraic Riccati Equation (CARE) |
| [`CartesianMap`](inc/damp/kinematics/motion_maps.hpp#L44) | block | Kinematics / pose | Cartesian gantry: independent per-axis affine map `task = scale·act + offset` (the "kinematics" is the identity, exposed for a uniform forward/inverse interface) |
| [`CartesianMove`](inc/damp/trajectory/cartesian_move.hpp#L114) | block | Trajectory value types | Path-preserving task-space move (Pipeline B / LIN) |
| [`Cascade`](inc/damp/motor/cascade.hpp#L94) | block | Motor control pack (if present) | Position→velocity→torque cascade: the default Servo outer controller |
| [`Cascade`](inc/damp/controllers/composition.hpp#L49) | block | Runtime controllers | Series cascade of two SISO controllers: outer → inner reference |
| [`CascadeBandwidths`](inc/damp/motor/cascade.hpp#L57) | block | Motor control pack (if present) | The three bandwidth knobs of the position/velocity/current cascade |
| [`CascadeConfig`](inc/damp/motor/cascade.hpp#L73) | block | Motor control pack (if present) | Configuration for a Cascade controller (the mechanical / outer-loop half) |
| [`catmull_1d`](inc/damp/toolbox/lookup.hpp#L339) | function | Embedded helpers (controls-adjacent utilities) | 1-D non-uniform Catmull-Rom (cubic Hermite) evaluation |
| [`cauer_thermal_ss`](inc/damp/motor/thermal.hpp#L109) | function | Motor control pack (if present) | Continuous state-space model of a physical Cauer RC thermal ladder |
| [`cauer_thermal_ss_ambient`](inc/damp/toolbox/thermal.hpp#L135) | function | Embedded helpers (controls-adjacent utilities) | Continuous Cauer RC ladder in absolute temperature with ambient input |
| [`cauer_thermal_ss_ambient_mimo`](inc/damp/toolbox/thermal.hpp#L162) | function | Embedded helpers (controls-adjacent utilities) | Continuous Cauer RC ladder, absolute temps, all nodes as outputs |
| [`cbf_relative_degree_1`](inc/damp/controllers/action_governor.hpp#L299) | function | Design-time synthesis (not PWM-rate) | Build a relative-degree-1 CBF inequality row |
| [`CBFConstraint`](inc/damp/controllers/action_governor.hpp#L283) | block | Runtime controllers | One affine row from a relative-degree-1 CBF condition |
| [`cbrt`](inc/damp/math/math.hpp#L92) | function | Scalar math & complex | Cube root (preserves sign for negative x) |
| [`ceil`](inc/damp/math/math.hpp#L329) | function | Scalar math & complex | Ceiling — smallest integer ≥ x |
| [`chirp`](inc/damp/estimation/excitation/chirp.hpp#L115) | function | Observers & estimators | Build a chirp design payload from a configuration |
| [`Chirp`](inc/damp/estimation/excitation/chirp.hpp#L132) | block | Observers & estimators | Linear or logarithmic chirp runtime generator |
| [`ChirpConfig`](inc/damp/estimation/excitation/chirp.hpp#L40) | block | Observers & estimators | Configuration for a sine chirp excitation |
| [`ChirpMode`](inc/damp/estimation/excitation/chirp.hpp#L26) | enum | Observers & estimators | Chirp sweep law |
| [`ChirpResult`](inc/damp/estimation/excitation/chirp.hpp#L83) | block | Observers & estimators | Chirp design payload |
| [`cholesky`](inc/damp/matrix/decomposition.hpp#L82) | function | Linear algebra | Cholesky decomposition for positive-definite matrices |
| [`cholesky_solve`](inc/damp/matrix/solve.hpp#L159) | function | Linear algebra | Solve AX = B via Cholesky (A = LLᴴ) |
| [`clamp_accel_asymmetric`](inc/damp/trajectory/online_otg.hpp#L101) | function | Trajectory value types | Asymmetric clamp of accel to +max_accel / −max_decel |
| [`clamp_current_circle`](inc/damp/motor/ipm.hpp#L78) | function | Motor control pack (if present) | Clamp a dq reference to the current circle of radius Imax (preserve angle) |
| [`clamp_jerk_accel_stop`](inc/damp/trajectory/online_otg.hpp#L179) | function | Trajectory value types | Clamp jerk so one step does not drive a past a_stop |
| [`clamp_rate_asymmetric`](inc/damp/trajectory/online_otg.hpp#L110) | function | Trajectory value types | Asymmetric clamp of accel (legacy name) |
| [`clarke_park_transform`](inc/damp/transforms.hpp#L444) | function | Motor control pack (if present) | Fused Clarke-Park transform (abc → dq) |
| [`clarke_park_zero_transform`](inc/damp/transforms.hpp#L527) | function | Motor control pack (if present) | Fused Clarke-Park transform with zero (abc → dq0) |
| [`clarke_transform`](inc/damp/transforms.hpp#L343) | function | Motor control pack (if present) | Clarke transform (abc → αβ) |
| [`clarke_zero_transform`](inc/damp/transforms.hpp#L273) | function | Motor control pack (if present) | Zero-retaining Clarke transform (abc → αβ0) |
| [`classical_dob`](inc/damp/estimation/dob.hpp#L338) | function | Observers & estimators | Synthesize a classical disturbance observer from a nominal plant and Q-filter |
| [`ClassicalDOB`](inc/damp/estimation/dob.hpp#L370) | block | Observers & estimators | Classical Pn^-1·Q disturbance observer runtime (bolt-on compensator) |
| [`ClassicalDobResult`](inc/damp/estimation/dob.hpp#L295) | block | Observers & estimators | Design result for the classical Pn^-1·Q disturbance observer |
| [`classify_range`](inc/damp/toolbox/conditioning.hpp#L396) | function | Embedded helpers (controls-adjacent utilities) | Classify x against the four band edges `[fault_lo (valid_lo, valid_hi) fault_hi]` (assumed ordered, non-decreasing) |
| [`closed_loop_poles`](inc/damp/design/stability.hpp#L335) | function | Design-time synthesis (not PWM-rate) | Compute closed-loop poles (eigenvalues) with state feedback |
| [`cohen_coon`](inc/damp/design/pid_design.hpp#L280) | function | Design-time synthesis (not PWM-rate) | Cohen-Coon tuning from first-order-plus-dead-time model |
| [`ColVec`](inc/damp/matrix/colvec.hpp#L27) | block | Linear algebra | Concrete Column vector specialization of Matrix<N, 1, T> |
| [`ColView`](inc/damp/matrix/views.hpp#L257) | block | Linear algebra | Non-owning column view of a matrix |
| [`comb_notch_window`](inc/damp/filters/moving_average.hpp#L34) | function | Filters & signal conditioning | Window length for a moving-average comb that notches f_notch and all its harmonics: N = round(fs / f_notch) |
| [`CommandProjectionResult`](inc/damp/controllers/action_governor.hpp#L185) | block | Design-time synthesis (not PWM-rate) | Result of an affine command projection |
| [`compensate_coulomb_viscous`](inc/damp/toolbox/conditioning.hpp#L178) | function | Embedded helpers (controls-adjacent utilities) | Coulomb + viscous friction compensation (static inverse, massless) |
| [`compensate_duties`](inc/damp/motor/dead_time.hpp#L227) | function | Motor control pack (if present) | Apply duty compensation to three half-bridge duties (clamped to [0, 1]) |
| [`compensate_phase_voltages`](inc/damp/motor/dead_time.hpp#L247) | function | Motor control pack (if present) | Apply phase-voltage compensation before modulation |
| [`compensate_stribeck`](inc/damp/toolbox/conditioning.hpp#L205) | function | Embedded helpers (controls-adjacent utilities) | Stribeck-style friction compensation (static + Coulomb + viscous) |
| [`compensator_from_feature`](inc/damp/estimation/successive_compensator.hpp#L179) | function | Observers & estimators | Map an FRF feature to RBJ biquad coefficients |
| [`CompensatorKind`](inc/damp/estimation/successive_compensator.hpp#L38) | enum | Observers & estimators | Which biquad family to load for a feature |
| [`Complementary`](inc/damp/filters/complementary.hpp#L32) | block | Filters & signal conditioning | Scalar (1-D) complementary filter — fuse a fast rate with a slow absolute |
| [`ComplementaryFilter`](inc/damp/estimation/sensor_fusion.hpp#L64) | block | Observers & estimators | Simple complementary filter for orientation estimation |
| [`complex`](inc/damp/math/complex.hpp#L38) | block | Scalar math & complex | Constexpr complex number class for compile-time computations |
| [`compute_eigenvalues`](inc/damp/matrix/eigen.hpp#L382) | function | Linear algebra | Compute the eigenvalues (and Schur vectors) of a real square matrix |
| [`constant_power_torque_limit`](inc/damp/motor/foc.hpp#L80) | function | Motor control pack (if present) | Constant-power torque ceiling Tₘₐₓ = P_rated / \|ω\| |
| [`ConstantInertiaFeedforward`](inc/damp/toolbox/actuator.hpp#L197) | block | Embedded helpers (controls-adjacent utilities) | Per-axis decoupled torque feedforward: `τ = J·a + b·v + τ_c·sign(v) + g` |
| [`continuous_lpf_exact_step`](inc/damp/filters/lowpass.hpp#L37) | function | Filters & signal conditioning | Exact step of continuous first-order LPF ẏ = −ω_c (y − u) |
| [`continuous_lqr`](inc/damp/controllers/lqr.hpp#L377) | function | Design-time synthesis (not PWM-rate) | Continuous-time Linear-Quadratic Regulator design (+1 more overload) |
| [`ContinuousPID`](inc/damp/controllers/pid.hpp#L633) | block | Runtime controllers | Continuous-gain PID with per-tick sample time (variable-rate secondary form) |
| [`controllability_gramian`](inc/damp/design/stability.hpp#L118) | function | Design-time synthesis (not PWM-rate) | Continuous/discrete controllability Gramian W_c |
| [`controllability_matrix`](inc/damp/design/stability.hpp#L52) | function | Design-time synthesis (not PWM-rate) | Compute the controllability matrix [B, AB, A²B, ..., A^(N-1)B] |
| [`ControlMode`](inc/damp/motor/cascade.hpp#L42) | enum | Motor control pack (if present) | Which cascade outer loops are active |
| [`Convention`](inc/damp/transforms.hpp#L69) | enum | Motor control pack (if present) | Scaling convention for the Clarke/Park family |
| [`copysign`](inc/damp/math/math.hpp#L406) | function | Scalar math & complex | Copy sign — magnitude of mag with the sign of sgn_src |
| [`CoreXY`](inc/damp/kinematics/motion_maps.hpp#L74) | block | Kinematics / pose | CoreXY belt mapping (2 motors A/B → Cartesian X/Y) |
| [`cos`](inc/damp/matrix/functions.hpp#L738) | function | Linear algebra | Matrix cosine via scaling and double-angle reconstruction |
| [`cos`](inc/damp/math/math.hpp#L169) | function | Scalar math & complex | Cosine |
| [`cosh`](inc/damp/matrix/functions.hpp#L809) | function | Linear algebra | Matrix hyperbolic cosine cosh(A) = (exp(A) + exp(−A))/2 |
| [`Counter`](inc/damp/toolbox/logic.hpp#L273) | block | Embedded helpers (controls-adjacent utilities) | Edge-counting up/down counter: increments on each rising edge of up, decrements on each rising edge of down. Returns the running count |
| [`cross_check_analog`](inc/damp/toolbox/conditioning.hpp#L573) | function | Embedded helpers (controls-adjacent utilities) | Cross-check two calibrated analog readings for redundant sensors |
| [`CTD`](inc/damp/toolbox/iec61131.hpp#L364) | block | Embedded helpers (controls-adjacent utilities) | CTD Counter (Count Down) |
| [`ctrb`](inc/damp/matlab.hpp#L258) | function | MATLAB®-style aliases (host) | MATLAB® short alias for controllability_matrix (+1 more overload) |
| [`CTU`](inc/damp/toolbox/iec61131.hpp#L324) | block | Embedded helpers (controls-adjacent utilities) | CTU Counter (Count Up) |
| [`CTUD`](inc/damp/toolbox/iec61131.hpp#L404) | block | Embedded helpers (controls-adjacent utilities) | CTUD Counter (Count Up Down) |
| [`cubic_rise`](inc/damp/trajectory/cam.hpp#L237) | function | Trajectory value types | Cubic rest–rest rise (Hermite; C¹ ends — accel jumps into a dwell) |
| [`cubic_spline`](inc/damp/trajectory/spline.hpp#L239) | function | Trajectory value types | Cubic (C²) spline through points at times; clamped end velocities |
| [`CukMap`](inc/damp/power/cuk.hpp#L30) | block | Power electronics pack (if present) | Ćuk CCM plant map — HardSwitchingConverter (layer 3) |
| [`current_loop_pi`](inc/damp/motor/foc.hpp#L63) | function | Motor control pack (if present) | Current-loop PI gains by closed-loop pole placement on the R–L plant |
| [`current_magnitude`](inc/damp/motor/ipm.hpp#L70) | function | Motor control pack (if present) | Stator current magnitude on the current circle |
| [`CurrentRegulatedStage`](inc/damp/power/concepts.hpp#L124) | concept | Power electronics pack (if present) | Current-regulated PE stage: control(r, y) → switch duties (layer 2) |
| [`CurrentStage`](inc/damp/power/current_stage.hpp#L58) | block | Power electronics pack (if present) | Current-regulated stage: PIController + map |
| [`cycloidal_rise`](inc/damp/trajectory/cam.hpp#L210) | function | Trajectory value types | Cycloidal rise of height h over unit leader fraction u ∈ [0, 1] |
| [`damp`](inc/damp/analysis/poles.hpp#L90) | function | Frequency-domain analysis (host) | Compute natural frequency and damping for each pole |
| [`damping_ratio_from_overshoot_percent`](inc/damp/design/pid_design.hpp#L510) | function | Design-time synthesis (not PWM-rate) | Map percent overshoot target to equivalent damping ratio |
| [`dare`](inc/damp/design/riccati.hpp#L634) | function | Design-time synthesis (not PWM-rate) | Solve the Discrete Algebraic Riccati Equation (DARE) |
| [`db2mag`](inc/damp/math/math.hpp#L470) | function | Scalar math & complex | Decibels to magnitude, 10^(db/20) |
| [`DcBusLimiter`](inc/damp/motor/limits.hpp#L66) | block | Motor control pack (if present) | Holds the inverter's torque current within DC-bus current/power limits |
| [`DcBusLimits`](inc/damp/motor/limits.hpp#L31) | block | Motor control pack (if present) | DC-bus current and voltage limits for an inverter |
| [`DcBusState`](inc/damp/motor/limits.hpp#L42) | block | Motor control pack (if present) | DC-bus state and the torque-current derate it implies |
| [`dcgain`](inc/damp/analysis/norms.hpp#L40) | function | Frequency-domain analysis (host) | Compute DC gain of a continuous-time system |
| [`DCM`](inc/damp/math/geometry.hpp#L63) | block | Scalar math & complex | Direction cosine matrix — 3×3 rotation (SO(3) wrapper over Mat3) |
| [`dead_time_current_sign`](inc/damp/motor/dead_time.hpp#L111) | function | Motor control pack (if present) | Soft sign used by complementary blanking: ∈ {−1, 0, +1} |
| [`dead_time_duty_compensation`](inc/damp/motor/dead_time.hpp#L182) | function | Motor control pack (if present) | Duty correction δd to add to a commanded half-bridge duty |
| [`dead_time_phase_voltage_compensation`](inc/damp/motor/dead_time.hpp#L206) | function | Motor control pack (if present) | Phase-voltage correction [V] dual of dead_time_duty_compensation |
| [`dead_time_pole_voltage_error`](inc/damp/motor/dead_time.hpp#L145) | function | Motor control pack (if present) | Average pole-voltage error Δv = v_avg − d V_dc |
| [`dead_time_ratio`](inc/damp/motor/dead_time.hpp#L97) | function | Motor control pack (if present) | Dead-time ratio r = T_d / T_s (fraction of the PWM period) |
| [`deadband`](inc/damp/toolbox/conditioning.hpp#L47) | function | Embedded helpers (controls-adjacent utilities) | Dead zone over `[lower, upper]`, matching Simulink®'s Dead Zone block (+1 more overload) |
| [`DeadTimeBlanking`](inc/damp/motor/dead_time.hpp#L78) | enum | Motor control pack (if present) | Blanking geometry assumed by the average error model |
| [`DeadTimeCompensator`](inc/damp/motor/dead_time.hpp#L279) | block | Motor control pack (if present) | Runtime dead-time compensator (holds blanking parameters) |
| [`Debounce`](inc/damp/toolbox/logic.hpp#L213) | block | Embedded helpers (controls-adjacent utilities) | Debounce: the output adopts in only after in differs from the current output continuously for stable_time. Rejects contact bounce and brief glitches. (Not an IEC block — the one everyone hand-rolls.) |
| [`default_tol`](inc/damp/matrix/matrix_traits.hpp#L86) | function | Linear algebra | Type-appropriate default tolerance for floating-point comparisons. float  ~7 decimal digits  → 1e-6 double ~15 decimal digits → 1e-12 |
| [`deg2rad`](inc/damp/math/math.hpp#L492) | function | Scalar math & complex | Degrees to radians, deg·π/180 |
| [`Delay`](inc/damp/filters/delay.hpp#L27) | block | Filters & signal conditioning | Discrete-time delay buffer |
| [`DeltaForward`](inc/damp/kinematics/motion_maps.hpp#L129) | block | Kinematics / pose | Result of a delta forward solve: the end-effector pose (orientation fixed to identity — deltas are 3-DOF translational) + validity |
| [`DeltaInverse`](inc/damp/kinematics/motion_maps.hpp#L121) | block | Kinematics / pose | Result of a delta inverse solve: the three actuator values + reachability |
| [`derate_window`](inc/damp/motor/thermal.hpp#L40) | function | Motor control pack (if present) | A two-breakpoint derating curve: 1 below derate_start, 0 at cutoff |
| [`det`](inc/damp/matrix/functions.hpp#L182) | function | Linear algebra | Matrix determinant det(A) |
| [`DFF`](inc/damp/toolbox/iec61131.hpp#L461) | block | Embedded helpers (controls-adjacent utilities) | D Flip-Flop (edge-triggered data latch) |
| [`dh_joint_axis_in_parent`](inc/damp/kinematics/serial_arm.hpp#L505) | function | Kinematics / pose | Standard-DH joint axis in the parent frame: $`z_{j-1} = (0,0,1)`$ |
| [`DhChain`](inc/damp/kinematics/serial_arm.hpp#L108) | block | Kinematics / pose | An N-joint DH chain (the arm geometry) |
| [`DhJoint`](inc/damp/kinematics/serial_arm.hpp#L90) | block | Kinematics / pose | One joint's standard (distal) DH parameters and motion limits |
| [`diag`](inc/damp/matlab.hpp#L372) | function | MATLAB®-style aliases (host) | Returns a square diagonal matrix from the given array (+1 more overload) |
| [`Diagonal`](inc/damp/matrix/views.hpp#L49) | block | Linear algebra | Diagonal view of a square matrix |
| [`differential_drive`](inc/damp/toolbox/io.hpp#L340) | function | Embedded helpers (controls-adjacent utilities) | Differential (tank) drive mixer: throttle + turn → left/right |
| [`DirectQuadrature`](inc/damp/transforms.hpp#L84) | block | Motor control pack (if present) | Direct-quadrature (rotor-frame) component pair |
| [`DirectQuadratureZero`](inc/damp/transforms.hpp#L238) | block | Motor control pack (if present) | Direct-quadrature-zero (rotor-frame) component triple |
| [`Discrete`](inc/damp/simulation/integrator.hpp#L83) | block | Simulation / SIL harness (host) | Discrete-time integrator (no integration, just one step) |
| [`discrete_lqe`](inc/damp/estimation/kalman.hpp#L175) | function | Observers & estimators | Discrete linear-quadratic estimator (LQE) — dual-of-LQR spelling of kalman (+1 more overload) |
| [`discrete_lqg`](inc/damp/controllers/lqg.hpp#L120) | function | Design-time synthesis (not PWM-rate) | Discrete Linear-Quadratic-Gaussian regulator design |
| [`discrete_lqgi`](inc/damp/controllers/lqgi.hpp#L124) | function | Design-time synthesis (not PWM-rate) | Discrete LQG with integral action (LQI + Kalman) for output tracking |
| [`discrete_lqi`](inc/damp/controllers/lqi.hpp#L93) | function | Design-time synthesis (not PWM-rate) | Discrete Linear-Quadratic-Integral (LQI) design for output tracking |
| [`discrete_lqr`](inc/damp/controllers/lqr.hpp#L115) | function | Design-time synthesis (not PWM-rate) | Discrete-time Linear-Quadratic Regulator design |
| [`discrete_lqr_from_continuous`](inc/damp/controllers/lqr.hpp#L253) | function | Design-time synthesis (not PWM-rate) | Design discrete LQR from continuous-time system via discretization (+1 more overload) |
| [`DiscretePIDResult`](inc/damp/controllers/pid.hpp#L43) | block | Design-time synthesis (not PWM-rate) | Fixed-rate discrete PID coefficients (canonical deploy form) |
| [`DiscretizationMethod`](inc/damp/systems/discretization.hpp#L26) | enum | LTI systems (SS / TF / ZPK / discretize) | Discretization methods for continuous-time state-space systems |
| [`discretize`](inc/damp/systems/discretization.hpp#L218) | function | LTI systems (SS / TF / ZPK / discretize) | Discretize a continuous-time state-space system |
| [`discretize_lqr_cost`](inc/damp/controllers/lqr.hpp#L181) | function | Design-time synthesis (not PWM-rate) | Discretize a continuous LQR cost integral over one sample (Van Loan) |
| [`DLATCH`](inc/damp/toolbox/iec61131.hpp#L492) | block | Embedded helpers (controls-adjacent utilities) | D Latch (level-sensitive / transparent latch) |
| [`dlqr`](inc/damp/controllers/lqr.hpp#L298) | function | Design-time synthesis (not PWM-rate) | Discrete-time LQR design (MATLAB®-style short name) |
| [`dlqr`](inc/damp/matlab.hpp#L618) | function | MATLAB®-style aliases (host) | Discrete-time Linear-Quadratic Regulator design |
| [`dlyap`](inc/damp/design/lyapunov.hpp#L130) | function | Design-time synthesis (not PWM-rate) | Solve the discrete-time Lyapunov (Stein) equation A X Aᵀ − X + Q = 0 |
| [`dlyap`](inc/damp/matlab.hpp#L1062) | function | MATLAB®-style aliases (host) | MATLAB® alias for the discrete Lyapunov solve AXAᵀ−X+Q=0 |
| [`dob`](inc/damp/estimation/dob.hpp#L128) | function | Observers & estimators | Validate and package DOB configuration into a runtime-ready design result |
| [`DOB`](inc/damp/estimation/dob.hpp#L153) | block | Observers & estimators | Lightweight SISO disturbance observer runtime |
| [`DOBConfig`](inc/damp/estimation/dob.hpp#L31) | block | Observers & estimators | Configuration for a first-order disturbance observer |
| [`DP45`](inc/damp/simulation/integrator.hpp#L757) | block | Simulation / SIL harness (host) | Dormand-Prince 5(4) adaptive integrator — the ode45 pair |
| [`dpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L268) | function | Motor control pack (if present) | DPWM zero-sequence for a chosen discontinuous scheme |
| [`DqCommand`](inc/damp/motor/foc.hpp#L331) | block | Motor control pack (if present) | Result of FOController::current_controller(): the dq voltage command plus its saturation signals |
| [`DriftMonitorResult`](inc/damp/analysis/identification.hpp#L215) | block | Frequency-domain analysis (host) | Result of model-drift monitoring used for adaptation triggers |
| [`Drive`](inc/damp/motor/drive.hpp#L88) | block | Motor control pack (if present) | The motor-agnostic drive core: encoder tracker + a Motor stage + torque limiting |
| [`drive_nameplate_limits`](inc/damp/motor/drive_design.hpp#L66) | function | Motor control pack (if present) | Build voltage/current/base-speed limits from bus, current ceiling, and machine |
| [`DriveFeel`](inc/damp/trajectory/selector_governor.hpp#L69) | block | Trajectory value types | Stick / command feel knobs (deadband, expo, PIO pole frequency) |
| [`DriveNameplateLimits`](inc/damp/motor/drive_design.hpp#L41) | block | Motor control pack (if present) | Nameplate-derived drive limits at one current operating point for base speed |
| [`DrivePair`](inc/damp/toolbox/io.hpp#L315) | block | Embedded helpers (controls-adjacent utilities) | Left/right actuator pair (differential/arcade drive output) |
| [`DriveResult`](inc/damp/motor/drive.hpp#L46) | block | Motor control pack (if present) | Output of one electrical step: the switch duties plus realized torque and flags |
| [`DriveSelectorGovernor`](inc/damp/trajectory/selector_governor.hpp#L87) | block | Trajectory value types | Two-axis drive governor: longitudinal v and yaw rate ω |
| [`DsogiPll`](inc/damp/filters/pll.hpp#L282) | block | Filters & signal conditioning | Dual-SOGI three-phase positive-sequence PLL (DSOGI-PLL) |
| [`duties_from_phase_voltages`](inc/damp/motor/modulation.hpp#L360) | function | Motor control pack (if present) | Map phase voltages + zero-sequence to clamped half-bridge duties |
| [`EemfObserver`](inc/damp/motor/sensorless_eemf.hpp#L93) | block | Motor control pack (if present) | Extended-EMF sensorless observer (salient mid-speed) |
| [`EemfResult`](inc/damp/motor/sensorless_eemf.hpp#L32) | block | Motor control pack (if present) | EEMF design payload |
| [`eig`](inc/damp/matlab.hpp#L473) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the eigenvalues of a square matrix |
| [`EigenResult`](inc/damp/matrix/eigen.hpp#L34) | block | Linear algebra | Eigenvalue computation result |
| [`EKFMeasFn`](inc/damp/estimation/ekf.hpp#L82) | concept | Observers & estimators | Concept for EKF measurement functions |
| [`EKFStateFn`](inc/damp/estimation/ekf.hpp#L54) | concept | Observers & estimators | Concept for EKF state functions |
| [`electromagnetic_torque`](inc/damp/motor/foc.hpp#L248) | function | Motor control pack (if present) | Electromagnetic torque produced by a dq current (salient PMSM) |
| [`ElectronicGear`](inc/damp/trajectory/cam.hpp#L155) | block | Trajectory value types | Electronic gearing law q = r·θ + φ (PLCopen MC_GearIn) |
| [`EncoderTracker`](inc/damp/motor/encoder_tracker.hpp#L61) | block | Motor control pack (if present) | Kinematic angle/speed tracking observer (PLL) over a wrapped absolute position |
| [`EncoderTrackerConfig`](inc/damp/motor/encoder_tracker.hpp#L23) | block | Motor control pack (if present) | Configuration for EncoderTracker |
| [`enu_from_ned`](inc/damp/estimation/ins_mechanization.hpp#L116) | function | Observers & estimators | Map a NED vector into ENU (axis permute) |
| [`ErrorStateJacobian`](inc/damp/estimation/eskf.hpp#L247) | block | Observers & estimators | ESKF prediction Jacobians F and G (nominal state updated outside) |
| [`ErrorStateKalmanFilter`](inc/damp/estimation/eskf.hpp#L294) | block | Observers & estimators | Runtime error-state KF: tracks δx and P, not the full nominal state |
| [`ErrorTolerance`](inc/damp/simulation/solver.hpp#L259) | block | Simulation / SIL harness (host) | Scalar-or-vector local-error tolerance for adaptive ODE control |
| [`esc`](inc/damp/controllers/esc.hpp#L153) | function | Design-time synthesis (not PWM-rate) | Synthesize an extremum-seeking controller |
| [`esc_lpf_alpha`](inc/damp/controllers/esc.hpp#L131) | function | Design-time synthesis (not PWM-rate) | First-order discrete low-pass coefficient for corner wc [rad/s] at Ts |
| [`esc_mppt`](inc/damp/controllers/esc.hpp#L200) | function | Design-time synthesis (not PWM-rate) | MPPT-flavored ESC: maximize a power measurement by perturbing the operating point (e.g. converter duty or reference voltage) |
| [`ESCConfig`](inc/damp/controllers/esc.hpp#L59) | block | Runtime controllers | Extremum-seeking controller configuration (discrete realization) |
| [`ESCResult`](inc/damp/controllers/esc.hpp#L103) | block | Design-time synthesis (not PWM-rate) | Design result for the extremum-seeking controller |
| [`eskf_imu`](inc/damp/estimation/eskf.hpp#L179) | function | Observers & estimators | Design Q, R, P₀ for a 6-state IMU attitude ESKF (gyro + accel) |
| [`eskf_marg`](inc/damp/estimation/eskf.hpp#L218) | function | Observers & estimators | Design Q, R, P₀ for a 6-state MARG attitude ESKF (gyro + accel + mag) |
| [`ESKFMeasFn`](inc/damp/estimation/eskf.hpp#L270) | concept | Observers & estimators | Callable that returns a measurement linearization for one update |
| [`ESKFOrientationFilter`](inc/damp/estimation/sensor_fusion.hpp#L399) | block | Observers & estimators | Turnkey 6-state attitude ESKF: owns q, b_g, and the error filter |
| [`ESKFPredictFn`](inc/damp/estimation/eskf.hpp#L259) | concept | Observers & estimators | Callable that returns ErrorStateJacobian for one predict step |
| [`ESKFResult`](inc/damp/estimation/eskf.hpp#L50) | block | Observers & estimators | Design payload for an ESKF: Q, R, P₀, and success |
| [`estim`](inc/damp/matlab.hpp#L490) | function | MATLAB®-style aliases (host) | Form state estimator from system and estimator gain |
| [`Euler`](inc/damp/math/geometry.hpp#L205) | block | Scalar math & complex | Intrinsic Euler angles for sequence Order (default ZYX aerospace YPR) |
| [`EulerOrder`](inc/damp/math/geometry.hpp#L39) | enum | Scalar math & complex | Intrinsic Euler angle sequence (axis order of successive principal rotations) |
| [`eval_frf`](inc/damp/systems/state_space.hpp#L166) | function | LTI systems (SS / TF / ZPK / discretize) | Evaluate frequency response of a state-space system |
| [`EventSimConfig`](inc/damp/simulation/hybrid.hpp#L67) | block | Simulation / SIL harness (host) | Configuration for event-driven hybrid / multi-rate style runs |
| [`Exact`](inc/damp/simulation/integrator.hpp#L117) | block | Simulation / SIL harness (host) | Exact integrator for LTI systems |
| [`exact_lti_step`](inc/damp/simulation/hybrid.hpp#L315) | function | Simulation / SIL harness (host) | Exact open-loop step of a single continuous LTI system |
| [`exp`](inc/damp/math/math.hpp#L232) | function | Scalar math & complex | Exponential function |
| [`ExperimentSafety`](inc/damp/estimation/experiment_safety.hpp#L62) | block | Observers & estimators | Experiment lifecycle + command protection for on-target commissioning |
| [`ExperimentSafetyConfig`](inc/damp/estimation/experiment_safety.hpp#L39) | block | Observers & estimators | Configuration for ExperimentSafety |
| [`expm`](inc/damp/matrix/functions.hpp#L323) | function | Linear algebra | Matrix exponential via scaling and squaring with Padé approximant of degree 13 |
| [`expo`](inc/damp/toolbox/conditioning.hpp#L252) | function | Embedded helpers (controls-adjacent utilities) | Exponential response curve `y = (1−k)·x + k·x³` (RC "expo") |
| [`ExtendedKalmanFilter`](inc/damp/estimation/ekf.hpp#L112) | block | Observers & estimators | Extended Kalman Filter for nonlinear discrete-time systems |
| [`extract_channel`](inc/damp/simulation/plot_plotly.hpp#L223) | function | Simulation / SIL harness (host) | Extract the i-th element from each vector entry into a std::vector<double> |
| [`extract_frf_features`](inc/damp/estimation/frequency_response.hpp#L740) | function | Observers & estimators | Extract peaks and valleys in one pass over the FRF table |
| [`extract_modes`](inc/damp/estimation/frequency_response.hpp#L715) | function | Observers & estimators | Extract resonant peaks (local \|G\| maxima) from an FRF table |
| [`extract_valleys`](inc/damp/estimation/frequency_response.hpp#L729) | function | Observers & estimators | Extract anti-resonance valleys (local \|G\| minima) from an FRF table |
| [`Extrapolation`](inc/damp/toolbox/lookup.hpp#L36) | enum | Embedded helpers (controls-adjacent utilities) | Out-of-range behaviour for a Lut1D / Lut2D / Lut3D query beyond its breakpoints |
| [`ExtremumSeekingController`](inc/damp/controllers/esc.hpp#L238) | block | Runtime controllers | Extremum-seeking controller runtime (model-free online optimizer) |
| [`ExtremumType`](inc/damp/controllers/esc.hpp#L47) | enum | Runtime controllers | Whether ESC climbs to a maximum or descends to a minimum of the objective |
| [`eye`](inc/damp/matlab.hpp#L413) | function | MATLAB®-style aliases (host) | Create an identity matrix of size n x n |
| [`F_TRIG`](inc/damp/toolbox/iec61131.hpp#L153) | block | Embedded helpers (controls-adjacent utilities) | F_TRIG (Falling Edge Trigger) |
| [`factorial`](inc/damp/trajectory/trajectory_types.hpp#L206) | function | Trajectory value types | k! (exact for modest k in float/double; used up through nonic BVP / jets) |
| [`falling_factorial`](inc/damp/trajectory/trajectory_types.hpp#L217) | function | Trajectory value types | Falling factorial i·(i−1)···(i−k+1) = i! / (i−k)! — the k-th derivative coefficient of tⁱ. Zero when i < k |
| [`FallingEdge`](inc/damp/toolbox/logic.hpp#L44) | block | Embedded helpers (controls-adjacent utilities) | Falling-edge detector: true on the tick x goes true → false |
| [`feedback`](inc/damp/systems/state_space.hpp#L311) | function | LTI systems (SS / TF / ZPK / discretize) | Negative feedback: y = sys1(u − sys2(y)) (+1 more overload) |
| [`FetLossModel`](inc/damp/motor/thermal.hpp#L161) | block | Motor control pack (if present) | First-order inverter FET loss model (conduction + switching) |
| [`field_for_speed`](inc/damp/motor/sepex.hpp#L134) | function | Motor control pack (if present) | Maximum field current that keeps e_a ≤ V_{a,max} at speed |
| [`field_weakening_id`](inc/damp/motor/field_weakening.hpp#L58) | function | Motor control pack (if present) | Feedforward field-weakening d-axis current from the voltage ellipse |
| [`FieldExcitation`](inc/damp/motor/sepex.hpp#L240) | enum | Motor control pack (if present) | How the field winding is powered |
| [`FieldWeakening`](inc/damp/motor/field_weakening.hpp#L141) | block | Motor control pack (if present) | Field-weakening current-reference regulator (voltage-feedback or feedforward) |
| [`FieldWeakeningConfig`](inc/damp/motor/field_weakening.hpp#L89) | block | Motor control pack (if present) | Configuration for FieldWeakening |
| [`FieldWeakeningPolicy`](inc/damp/motor/field_weakening.hpp#L205) | concept | Motor control pack (if present) | Concept for a pluggable field-weakening / current-reference policy |
| [`finish_extremum`](inc/damp/estimation/frequency_response.hpp#L561) | function | Observers & estimators | Half-power (peaks) or double-power (valleys) bandwidth → ζ |
| [`finite_non_negative`](inc/damp/math/math.hpp#L448) | function | Scalar math & complex | True if x is finite and non-negative (x ≥ 0) |
| [`finite_positive`](inc/damp/math/math.hpp#L438) | function | Scalar math & complex | True if x is finite and strictly greater than zero |
| [`fir1`](inc/damp/filters/fir.hpp#L218) | function | Filters & signal conditioning | Window-method FIR design (normalized frequency, single cutoff) (+1 more overload) |
| [`fir1_hz`](inc/damp/filters/fir.hpp#L372) | function | Filters & signal conditioning | Window-method FIR design from frequencies in Hz (+1 more overload) |
| [`fir_window`](inc/damp/filters/fir.hpp#L408) | function | Filters & signal conditioning | Alias for fir1 (descriptive name) (+1 more overload) |
| [`FirDesignResult`](inc/damp/filters/fir.hpp#L67) | block | Filters & signal conditioning | Window-method FIR design result |
| [`fires_at`](inc/damp/simulation/multirate.hpp#L73) | function | Simulation / SIL harness (host) | True when base tick tick is a firing instant for divisor |
| [`FirFilter`](inc/damp/filters/fir.hpp#L450) | block | Filters & signal conditioning | Direct-form FIR filter runtime (tapped delay line + dot product) |
| [`FirstOrderCoeffs`](inc/damp/filters/iir_design.hpp#L30) | block | Filters & signal conditioning | DSP coefficients for first-order IIR filter |
| [`FirstOrderPlantEstimator`](inc/damp/estimation/parameter_estimation.hpp#L126) | block | Observers & estimators | Online estimator for a first-order plant's gain and time constant |
| [`FirstOrderPlantEstimatorConfig`](inc/damp/estimation/parameter_estimation.hpp#L85) | block | Observers & estimators | Configuration for the first-order grey-box estimator |
| [`FirType`](inc/damp/filters/fir.hpp#L38) | enum | Filters & signal conditioning | Ideal frequency-selective response for fir1 / fir_window |
| [`FirWindow`](inc/damp/filters/fir.hpp#L50) | enum | Filters & signal conditioning | Window applied to the ideal truncated impulse response |
| [`fit_second_order_step`](inc/damp/estimation/successive_compensator.hpp#L85) | function | Observers & estimators | Fit one underdamped second-order mode from a step / ring-down capture |
| [`FitMetrics`](inc/damp/analysis/identification.hpp#L138) | block | Frequency-domain analysis (host) | Scalar fit metrics used for model comparison and gating |
| [`five_bar_symmetric`](inc/damp/kinematics/scara.hpp#L241) | function | Kinematics / pose | Build a symmetric five-bar parallel SCARA |
| [`FiveBar`](inc/damp/kinematics/scara.hpp#L92) | block | Kinematics / pose | Planar five-bar parallel manipulator (parallel SCARA) |
| [`FiveBarForward`](inc/damp/kinematics/scara.hpp#L76) | block | Kinematics / pose | Forward-kinematics result: the end-effector point (x, y) + validity |
| [`FiveBarGeometry`](inc/damp/kinematics/scara.hpp#L57) | block | Kinematics / pose | Symmetric five-bar geometry (two base motors, equal proximal/distal links) |
| [`FiveBarInverse`](inc/damp/kinematics/scara.hpp#L69) | block | Kinematics / pose | Inverse-kinematics result: the two motor angles [rad] + reachability |
| [`fixed_solve`](inc/damp/simulation/solver.hpp#L243) | function | Simulation / SIL harness (host) | Fixed-step solve in one call |
| [`FixedStepSolver`](inc/damp/simulation/solver.hpp#L118) | block | Simulation / SIL harness (host) | Fixed-step ODE solver (+1 more overload) |
| [`floor`](inc/damp/math/math.hpp#L317) | function | Scalar math & complex | Floor — largest integer ≤ x |
| [`flux_for_speed`](inc/damp/motor/acim.hpp#L163) | function | Motor control pack (if present) | Rotor-flux command that keeps referred flux inside the voltage circle |
| [`flux_from_id`](inc/damp/motor/acim.hpp#L99) | function | Motor control pack (if present) | Steady-state rotor flux from magnetizing current: λ_r = Lₘ i_d [Wb] |
| [`flux_from_Kv`](inc/damp/motor/foc.hpp#L184) | function | Motor control pack (if present) | PM flux linkage from the datasheet velocity constant Kᵥ |
| [`flux_from_torque_constant`](inc/damp/motor/foc.hpp#L127) | function | Motor control pack (if present) | PM flux linkage from a motor's torque constant (amplitude-invariant) |
| [`FluxSensorlessAdapter`](inc/damp/motor/sensorless.hpp#L104) | block | Motor control pack (if present) | Adapter: SensorlessEstimator → SensorlessObserver |
| [`fmod`](inc/damp/math/math.hpp#L361) | function | Scalar math & complex | Floating-point remainder, x − y·trunc(x/y) (sign of x), matching std::fmod's truncated-quotient convention |
| [`foc_autotune`](inc/damp/motor/foc_autotune.hpp#L225) | function | Motor control pack (if present) | Short alias for the nameplate FOC cascade designer |
| [`foc_cascade_tune`](inc/damp/motor/foc_autotune.hpp#L132) | function | Motor control pack (if present) | Nameplate / post-calibration FOC cascade tune (no experiment) |
| [`foc_cascade_tune_velocity_amigo`](inc/damp/motor/foc_autotune.hpp#L194) | function | Motor control pack (if present) | Velocity-loop PI from a biased-relay ultimate point (current loop already closed) |
| [`foc_velocity_amigo`](inc/damp/motor/foc_cascade_experiment.hpp#L236) | function | Motor control pack (if present) | Build a validated velocity-AMIGO experiment payload |
| [`foc_velocity_frf`](inc/damp/motor/foc_cascade_experiment.hpp#L558) | function | Motor control pack (if present) | Build a validated velocity-FRF experiment payload |
| [`FocCascadeExperimentOutput`](inc/damp/motor/foc_cascade_experiment.hpp#L70) | block | Motor control pack (if present) | Per-tick output of a FOC cascade velocity experiment |
| [`FocCascadeExperimentStatus`](inc/damp/motor/foc_cascade_experiment.hpp#L56) | enum | Motor control pack (if present) | Lifecycle of a FOC cascade velocity experiment |
| [`FocCascadeTuneResult`](inc/damp/motor/foc_autotune.hpp#L78) | block | Motor control pack (if present) | FOC cascade design result (current + velocity + position) |
| [`FocElectricalParams`](inc/damp/motor/foc_cascade_experiment.hpp#L89) | block | Motor control pack (if present) | Electrical nameplate inputs shared by velocity experiment drivers |
| [`FOController`](inc/damp/motor/foc.hpp#L346) | block | Motor control pack (if present) | Field-oriented current controller (dq PI + decoupling + SVPWM duties) |
| [`FocResult`](inc/damp/motor/foc.hpp#L312) | block | Motor control pack (if present) | Result of one FOController::step(), carrying the actuator command plus the saturation/measurement signals an outer (velocity/position) loop needs to propagate anti-windup back up a cascade |
| [`FocVelocityAmigoConfig`](inc/damp/motor/foc_cascade_experiment.hpp#L179) | block | Motor control pack (if present) | Configuration for FocVelocityAmigoExperiment |
| [`FocVelocityAmigoExperiment`](inc/damp/motor/foc_cascade_experiment.hpp#L275) | block | Motor control pack (if present) | Runtime step + relay experiment → FOC cascade velocity AMIGO gains |
| [`FocVelocityAmigoResult`](inc/damp/motor/foc_cascade_experiment.hpp#L221) | block | Motor control pack (if present) | Design payload for the velocity-AMIGO experiment (validated config) |
| [`FocVelocityFrfConfig`](inc/damp/motor/foc_cascade_experiment.hpp#L506) | block | Motor control pack (if present) | Configuration for FocVelocityFrfExperiment |
| [`FocVelocityFrfExperiment`](inc/damp/motor/foc_cascade_experiment.hpp#L578) | block | Motor control pack (if present) | Runtime stepped-sine FRF experiment → FOC cascade velocity PI gains |
| [`FocVelocityFrfResult`](inc/damp/motor/foc_cascade_experiment.hpp#L541) | block | Motor control pack (if present) | Validated design payload for the velocity FRF experiment |
| [`FOPDTModel`](inc/damp/analysis/identification.hpp#L76) | block | Frequency-domain analysis (host) | First-order plus dead-time candidate model |
| [`forward_substitute`](inc/damp/matrix/solve.hpp#L37) | function | Linear algebra | Forward substitution for Lx = b |
| [`ForwardEuler`](inc/damp/simulation/integrator.hpp#L169) | block | Simulation / SIL harness (host) | Forward Euler integrator |
| [`foster_thermal_ss`](inc/damp/motor/thermal.hpp#L72) | function | Motor control pack (if present) | Continuous state-space model of a Foster RC thermal network |
| [`FourSwitchBuckBoostMap`](inc/damp/power/four_switch_buck_boost.hpp#L33) | block | Power electronics pack (if present) | Four-switch buck–boost region map — voltage → Q1/Q3 duties (layer 3) |
| [`frame_orientations`](inc/damp/kinematics/serial_arm.hpp#L546) | function | Kinematics / pose | Orientations of DH frames {0,…,N} at configuration q |
| [`FrequencyPoint`](inc/damp/analysis/frequency.hpp#L38) | block | Frequency-domain analysis (host) | Single-point frequency response result |
| [`FrequencyResponseEstimator`](inc/damp/estimation/frequency_response.hpp#L106) | block | Observers & estimators | Lock-in FRF estimator over a fixed table of NFreq bins |
| [`FrequencyResponsePoint`](inc/damp/analysis/identification.hpp#L52) | block | Frequency-domain analysis (host) | One frequency-response point used by sweep/chirp-based identification |
| [`FrequencyResponseSummary`](inc/damp/analysis/identification.hpp#L64) | block | Frequency-domain analysis (host) | Summary statistics for a frequency-response data set |
| [`frf_pid_autotune`](inc/damp/estimation/commissioning_frontends.hpp#L152) | function | Observers & estimators | PID gains from an open-loop FRF table + margin / bandwidth targets |
| [`frf_s_or_z`](inc/damp/analysis/frequency.hpp#L974) | function | Frequency-domain analysis (host) | FRF evaluation point: jω (continuous) or e^{jωTs} (discrete) |
| [`FRFEstimateResult`](inc/damp/analysis/identification.hpp#L116) | block | Frequency-domain analysis (host) | Frequency-response estimate quality summary |
| [`FrfFeatureKind`](inc/damp/estimation/frequency_response.hpp#L373) | enum | Observers & estimators | Peak (resonance) or valley (anti-resonance) on an FRF magnitude curve |
| [`FrfFeatureResult`](inc/damp/estimation/frequency_response.hpp#L453) | block | Observers & estimators | Combined peak + valley reduction of an FRF table |
| [`FrfMargins`](inc/damp/estimation/frequency_response.hpp#L485) | block | Observers & estimators | Open-loop margins read from a discrete FRF table (on-target) |
| [`FrfPidAutotuneConfig`](inc/damp/estimation/commissioning_frontends.hpp#L109) | block | Observers & estimators | Configuration for FRF-based PID gain selection |
| [`FrfPidAutotuneResult`](inc/damp/estimation/commissioning_frontends.hpp#L125) | block | Observers & estimators | FRF PID autotune hand-off (gains + success) |
| [`FrfPoint`](inc/damp/estimation/frequency_response.hpp#L60) | block | Observers & estimators | One measured frequency-response point (on-target FRF table entry) |
| [`frobenius_norm`](inc/damp/matrix/functions.hpp#L96) | function | Linear algebra | Frobenius norm ‖A‖F = √(Σᵢⱼ \|aᵢⱼ\|²) |
| [`full_qr`](inc/damp/matrix/decomposition.hpp#L291) | function | Linear algebra | Full QR factorization via Householder reflections (real or complex T) |
| [`FullBridge1phMap`](inc/damp/power/full_bridge.hpp#L56) | block | Power electronics pack (if present) | 1φ full-bridge map — voltage command → two high-side duties (layer 3) |
| [`FullBridgePwm`](inc/damp/power/full_bridge.hpp#L33) | enum | Power electronics pack (if present) | 1φ full-bridge carrier style |
| [`FullQR`](inc/damp/matrix/decomposition.hpp#L273) | block | Linear algebra | Result of a full (complete) QR factorization |
| [`FwMethod`](inc/damp/motor/field_weakening.hpp#L79) | enum | Motor control pack (if present) | Field-weakening law selection |
| [`gain_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L202) | function | Frequency-domain analysis (host) | Find gain margin using unwrapped phase trajectory |
| [`GearFollower`](inc/damp/trajectory/cam.hpp#L540) | block | Trajectory value types | Runtime electronic gear q = k_q(r θ + φ) + o_q with chain-rule rates |
| [`geomspace`](inc/damp/analysis/linspace.hpp#L88) | function | Frequency-domain analysis (host) | Geometric sequence from start to end (endpoint values, not exponents) |
| [`Goertzel`](inc/damp/filters/spectral.hpp#L46) | block | Filters & signal conditioning | Generalized Goertzel single-bin DFT — amplitude/phase at one frequency |
| [`gram`](inc/damp/matlab.hpp#L1109) | function | MATLAB®-style aliases (host) | MATLAB® alias for the controllability/observability Gramian of a system |
| [`gravity_nav`](inc/damp/estimation/ins_mechanization.hpp#L93) | function | Observers & estimators | Gravity vector in the local-level frame [m/s²] |
| [`GreyBoxCandidateModels`](inc/damp/analysis/identification.hpp#L172) | block | Frequency-domain analysis (host) | Holds all candidate models produced during an identification pass |
| [`GreyBoxIdentificationResult`](inc/damp/analysis/identification.hpp#L188) | block | Frequency-domain analysis (host) | Selected model used by downstream model-based controller design |
| [`h_pattern`](inc/damp/toolbox/io.hpp#L377) | function | Embedded helpers (controls-adjacent utilities) | H-pattern (dual-lever / tank) drive: one lever per track, no mixing |
| [`HalfBridge1phMap`](inc/damp/power/half_bridge.hpp#L40) | block | Power electronics pack (if present) | Half-bridge pole map — voltage command → one high-side duty (layer 3) |
| [`HallDecoder`](inc/damp/motor/commutation.hpp#L69) | block | Motor control pack (if present) | 120° Hall-effect decoder: sector, coarse angle, direction, and speed |
| [`HallState`](inc/damp/motor/commutation.hpp#L40) | block | Motor control pack (if present) | One three-phase Hall reading decoded to an electrical sector and its kinematics |
| [`hankel_singular_values`](inc/damp/design/model_reduction.hpp#L596) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for hankelsv |
| [`hankelsv`](inc/damp/design/model_reduction.hpp#L439) | function | Design-time synthesis (not PWM-rate) | Hankel singular values of a stable state-space system |
| [`HardSwitchingConverter`](inc/damp/power/concepts.hpp#L94) | concept | Power electronics pack (if present) | Voltage-shaped command → duties (topology / plant map plug) |
| [`harmonic_suppressor`](inc/damp/controllers/harmonic_suppression.hpp#L80) | function | Design-time synthesis (not PWM-rate) | Synthesize a multi-resonant harmonic suppressor |
| [`HarmonicAnalyzer`](inc/damp/filters/spectral.hpp#L173) | block | Filters & signal conditioning | Harmonic analyzer — a Goertzel bank over a fundamental and K−1 harmonics |
| [`HarmonicSuppressor`](inc/damp/controllers/harmonic_suppression.hpp#L121) | block | Runtime controllers | Multi-resonant harmonic suppressor — a parallel bank of PR resonators |
| [`HarmonicSuppressorResult`](inc/damp/controllers/harmonic_suppression.hpp#L45) | block | Design-time synthesis (not PWM-rate) | Design result for a multi-resonant harmonic suppressor |
| [`HasSwitchDuties`](inc/damp/power/concepts.hpp#L76) | concept | Power electronics pack (if present) | Duck-typed switch command: duty vector + clip flag |
| [`heading_from_baseline_nav`](inc/damp/estimation/ins_eskf.hpp#L536) | function | Observers & estimators | Heading [rad] from a nav-frame baseline vector (dual-antenna difference) |
| [`Heun`](inc/damp/simulation/integrator.hpp#L840) | block | Simulation / SIL harness (host) | Heun's method (Improved Euler, RK2) integrator |
| [`hfi`](inc/damp/motor/sensorless_hfi.hpp#L79) | function | Motor control pack (if present) | Design pulsating HFI from nameplate / PWM headroom |
| [`HfInjector`](inc/damp/motor/sensorless_hfi.hpp#L142) | block | Motor control pack (if present) | High-frequency voltage injector (estimated dq) |
| [`HfiObserver`](inc/damp/motor/sensorless_hfi.hpp#L193) | block | Motor control pack (if present) | Pulsating-HFI saliency tracker (pre-polarity: locks θ or θ+π) |
| [`HfiPolarityIdent`](inc/damp/motor/sensorless_hfi.hpp#L318) | block | Motor control pack (if present) | Polarity ID via short +d voltage pulse after saliency lock |
| [`HfiResult`](inc/damp/motor/sensorless_hfi.hpp#L39) | block | Motor control pack (if present) | Pulsating-HFI design payload |
| [`HighPass`](inc/damp/filters/highpass.hpp#L40) | block | Filters & signal conditioning | First-order high-pass (washout) filter runtime |
| [`highpass_1st`](inc/damp/filters/iir_design.hpp#L169) | function | Filters & signal conditioning | First-order high-pass filter design (Tustin / bilinear) |
| [`highpass_2nd`](inc/damp/filters/iir_design.hpp#L684) | function | Filters & signal conditioning | Second-order high-pass filter (RBJ) |
| [`highshelf`](inc/damp/filters/iir_design.hpp#L765) | function | Filters & signal conditioning | High-shelf EQ filter: boost or cut everything above fc |
| [`hinfnorm`](inc/damp/matlab.hpp#L1131) | function | MATLAB®-style aliases (host) | MATLAB® alias for the H∞ system norm norm(sys,Inf) / hinfnorm(sys) |
| [`holonomic_drive`](inc/damp/toolbox/io.hpp#L495) | function | Embedded helpers (controls-adjacent utilities) | Four-wheel holonomic drive mixer (mecanum / omni) |
| [`HolonomicOutput`](inc/damp/toolbox/io.hpp#L463) | block | Embedded helpers (controls-adjacent utilities) | Four-wheel actuator set (mecanum/holonomic drive output), X-config |
| [`HolonomicSelectorGovernor`](inc/damp/trajectory/selector_governor.hpp#L189) | block | Trajectory value types | Three-axis holonomic governor: forward v_y, strafe v_x, yaw ω |
| [`HybridEventAction`](inc/damp/simulation/hybrid.hpp#L102) | block | Simulation / SIL harness (host) | Action applied when an event fires (after Exact advance to the event time) |
| [`HybridSensorlessFrontEnd`](inc/damp/motor/sensorless_hfi.hpp#L399) | block | Motor control pack (if present) | Hybrid HFI + mid-speed front-end |
| [`HybridSimulationResult`](inc/damp/simulation/hybrid.hpp#L80) | block | Simulation / SIL harness (host) | Result of a piecewise-LTI hybrid simulation |
| [`hypot`](inc/damp/math/math.hpp#L76) | function | Scalar math & complex | Euclidean distance hypot(x, y) = √(x² + y²), without overflow |
| [`Hysteresis`](inc/damp/toolbox/conditioning.hpp#L332) | block | Embedded helpers (controls-adjacent utilities) | Hysteresis comparator (Schmitt trigger): bool output with separate on/off thresholds to reject chatter |
| [`ia_from_torque`](inc/damp/motor/sepex.hpp#L101) | function | Motor control pack (if present) | Armature current for a torque command at a field current |
| [`IdentificationHandoff`](inc/damp/analysis/identification.hpp#L202) | block | Frequency-domain analysis (host) | Aggregated handoff payload from identification to model-based control design |
| [`IdentificationModelKind`](inc/damp/analysis/identification.hpp#L161) | enum | Frequency-domain analysis (host) | Enumerates which reduced model family is selected for downstream design |
| [`IdentifiedModelLike`](inc/damp/analysis/identification.hpp#L237) | concept | Frequency-domain analysis (host) | Concept for identified models consumed by downstream design modules |
| [`IfCommand`](inc/damp/motor/scalar_control.hpp#L156) | block | Motor control pack (if present) | Result of one IfController step: an open-loop angle and a dq current reference |
| [`IfController`](inc/damp/motor/scalar_control.hpp#L182) | block | Motor control pack (if present) | Forced-current (I-F) open-loop startup controller |
| [`impedance`](inc/damp/analysis/frequency.hpp#L999) | function | Frequency-domain analysis (host) | Compute impedance frequency response from a SISO admittance system |
| [`impedance_direct`](inc/damp/analysis/frequency.hpp#L1039) | function | Frequency-domain analysis (host) | Compute impedance frequency response from a SISO impedance transfer function |
| [`ImpedanceResult`](inc/damp/analysis/frequency.hpp#L872) | block | Frequency-domain analysis (host) | Result of impedance frequency response evaluation |
| [`impulse`](inc/damp/analysis/time_response.hpp#L185) | function | Frequency-domain analysis (host) | Impulse response of a (MIMO) state-space system (+1 more overload) |
| [`impulse`](inc/damp/estimation/excitation/impulse.hpp#L70) | function | Observers & estimators | Build an impulse design payload |
| [`Impulse`](inc/damp/estimation/excitation/impulse.hpp#L86) | block | Observers & estimators | Rectangular impulse (finite-width Dirac stand-in) |
| [`ImpulseConfig`](inc/damp/estimation/excitation/impulse.hpp#L32) | block | Observers & estimators | Configuration for a one-shot rectangular impulse (Dirac stand-in) |
| [`impulseplot`](inc/damp/simulation/plot_plotly.hpp#L733) | function | Simulation / SIL harness (host) | Plot an impulse response, one trace per input/output pair |
| [`ImpulseResult`](inc/damp/estimation/excitation/impulse.hpp#L52) | block | Observers & estimators | Impulse design payload |
| [`ImuMount`](inc/damp/estimation/serial_arm_pose.hpp#L88) | block | Observers & estimators | Fixed extrinsic of an IMU relative to its DH / link frame |
| [`ImuSample`](inc/damp/estimation/ins_mechanization.hpp#L125) | block | Observers & estimators | IMU sample in the body frame |
| [`in_jerk_terminal`](inc/damp/trajectory/online_otg.hpp#L153) | function | Trajectory value types | True when residual a is closing on rem and inside the jmax stop budget |
| [`InductorCurrentPIResult`](inc/damp/design/pid_design.hpp#L835) | block | Design-time synthesis (not PWM-rate) | Fixed-rate inductor current-loop PI (topology-agnostic) |
| [`infinity_norm`](inc/damp/matrix/functions.hpp#L35) | function | Linear algebra | Infinity norm ‖A‖∞: maximum absolute row sum |
| [`initial`](inc/damp/analysis/time_response.hpp#L219) | function | Frequency-domain analysis (host) | Initial-condition (free) response of a (MIMO) state-space system |
| [`input_shaper`](inc/damp/trajectory/input_shaper.hpp#L113) | function | Trajectory value types | Synthesize an input shaper for a second-order mode |
| [`InputOutputSample`](inc/damp/analysis/identification.hpp#L24) | block | Frequency-domain analysis (host) | One time-aligned sample used for closed-loop or open-loop identification |
| [`InputShaper`](inc/damp/trajectory/input_shaper.hpp#L195) | block | Trajectory value types | Input-shaper runtime — convolves a command stream with the shaper impulses |
| [`InputShaperBank`](inc/damp/trajectory/input_shaper.hpp#L256) | block | Trajectory value types | Multi-axis input-shaper bank — one shaper per axis, shared buffer length |
| [`InputShaperResult`](inc/damp/trajectory/input_shaper.hpp#L70) | block | Trajectory value types | Input-shaper design result: impulse amplitudes and sample delays |
| [`ins_aid_position`](inc/damp/estimation/ins_eskf.hpp#L503) | function | Observers & estimators | Predicted antenna / marker position with body lever-arm |
| [`ins_apply_correction`](inc/damp/estimation/ins_eskf.hpp#L489) | function | Observers & estimators | Inject filter correction into x and reset the error state |
| [`ins_error_jacobian`](inc/damp/estimation/ins_eskf.hpp#L307) | function | Observers & estimators | Discrete error-state transition and noise input for one IMU step |
| [`ins_eskf_design`](inc/damp/estimation/ins_eskf.hpp#L219) | function | Observers & estimators | Build Q / R / P0 for a 15-state INS ESKF from IMU noise densities |
| [`ins_heading`](inc/damp/estimation/ins_eskf.hpp#L521) | function | Observers & estimators | Horizontal heading [rad] of a body axis expressed in the nav frame |
| [`ins_inject`](inc/damp/estimation/ins_eskf.hpp#L433) | function | Observers & estimators | Inject δx into the nominal INS state (right-multiplicative q) |
| [`ins_predict`](inc/damp/estimation/ins_eskf.hpp#L466) | function | Observers & estimators | Mechanize nominal state and propagate the error covariance |
| [`ins_propagate_covariance`](inc/damp/estimation/ins_eskf.hpp#L403) | function | Observers & estimators | P ← F P Fᵀ + Q using sparse F (exact for ins_error_jacobian structure) |
| [`ins_update_heading`](inc/damp/estimation/ins_eskf.hpp#L638) | function | Observers & estimators | Sparse dual-antenna / yaw heading update (scalar) |
| [`ins_update_orientation`](inc/damp/estimation/ins_eskf.hpp#L676) | function | Observers & estimators | Sparse orientation update from a measured body→nav quaternion |
| [`ins_update_pose`](inc/damp/estimation/ins_eskf.hpp#L701) | function | Observers & estimators | Generic pose alias: position then orientation (two sparse updates) |
| [`ins_update_pose_heading`](inc/damp/estimation/ins_eskf.hpp#L720) | function | Observers & estimators | Pose alias: position + dual-antenna heading (two sparse updates) |
| [`ins_update_position`](inc/damp/estimation/ins_eskf.hpp#L551) | function | Observers & estimators | Sparse position update: y = p + Rℓ (+1 more overload) |
| [`ins_update_velocity`](inc/damp/estimation/ins_eskf.hpp#L597) | function | Observers & estimators | Sparse velocity update: y = v at the IMU origin (no lever-arm) |
| [`ins_update_zupt`](inc/damp/estimation/ins_eskf.hpp#L616) | function | Observers & estimators | Zero-velocity update (ZUPT): velocity measurement of zero |
| [`InsEskfResult`](inc/damp/estimation/ins_eskf.hpp#L180) | block | Observers & estimators | Design result for the 15-state INS error-state filter |
| [`InsNavigator`](inc/damp/estimation/ins_eskf.hpp#L754) | block | Observers & estimators | Runtime strapdown navigator: nominal InsState + 15-state ESKF |
| [`InsState`](inc/damp/estimation/ins_mechanization.hpp#L139) | block | Observers & estimators | Nominal strapdown navigation state (external to the ESKF error state) |
| [`instantaneous_power`](inc/damp/transforms.hpp#L605) | function | Motor control pack (if present) | Instantaneous active and reactive power from dq quantities (+1 more overload) |
| [`InstantaneousPower`](inc/damp/transforms.hpp#L566) | block | Motor control pack (if present) | Instantaneous active and reactive power |
| [`integrate_jerk_axis`](inc/damp/trajectory/online_otg.hpp#L372) | function | Trajectory value types | Integrate jerk → accel → value (semi-implicit), land/snap on value_target (+1 more overload) |
| [`IntegrationResult`](inc/damp/simulation/integrator.hpp#L31) | block | Simulation / SIL harness (host) | Result of an integration step |
| [`InteriorPointSolver`](inc/damp/design/qp.hpp#L1163) | block | Design-time synthesis (not PWM-rate) | Interior-point QP solver policy for damp::MPC |
| [`Interpolation`](inc/damp/toolbox/lookup.hpp#L42) | enum | Embedded helpers (controls-adjacent utilities) | In-range interpolant for Lut1D (2-D stays bilinear; use SplineSurface for C¹ grids) |
| [`inverse`](inc/damp/matrix/solve.hpp#L248) | function | Linear algebra | Matrix inverse A⁻¹ via mat::solve(A, I) |
| [`inverse_clarke_transform`](inc/damp/transforms.hpp#L360) | function | Motor control pack (if present) | Inverse Clarke transform (αβ → abc) |
| [`inverse_clarke_zero_transform`](inc/damp/transforms.hpp#L303) | function | Motor control pack (if present) | Inverse zero-retaining Clarke transform (αβ0 → abc) |
| [`inverse_deadband`](inc/damp/toolbox/conditioning.hpp#L79) | function | Embedded helpers (controls-adjacent utilities) | Inverse dead zone: add an offset to overcome a physical dead zone (valve overlap, static friction, motor stiction), with independent negative/positive offsets (+1 more overload) |
| [`inverse_deadband_sine`](inc/damp/toolbox/conditioning.hpp#L152) | function | Embedded helpers (controls-adjacent utilities) | Sine-shaped inverse dead zone (full boost once \|x\|≥w) |
| [`inverse_deadband_soft`](inc/damp/toolbox/conditioning.hpp#L128) | function | Embedded helpers (controls-adjacent utilities) | Smooth inverse dead zone via soft-sign (tanh-like Coulomb boost) |
| [`inverse_lerp`](inc/damp/toolbox/scaling.hpp#L48) | function | Embedded helpers (controls-adjacent utilities) | Inverse of lerp: the fraction t such that `lerp(a, b, t) == x` |
| [`inverse_park_clarke_transform`](inc/damp/transforms.hpp#L473) | function | Motor control pack (if present) | Fused inverse Park-Clarke transform (dq → abc) |
| [`inverse_park_clarke_zero_transform`](inc/damp/transforms.hpp#L540) | function | Motor control pack (if present) | Fused inverse Park-Clarke transform with zero (dq0 → abc) |
| [`inverse_park_transform`](inc/damp/transforms.hpp#L417) | function | Motor control pack (if present) | Inverse Park transform (dq → αβ) |
| [`inverse_park_zero_transform`](inc/damp/transforms.hpp#L508) | function | Motor control pack (if present) | Inverse Park transform with zero passthrough (dq0 → αβ0) |
| [`inverse_symmetrical_components`](inc/damp/transforms.hpp#L684) | function | Motor control pack (if present) | Inverse symmetrical-component transform (012 → abc) |
| [`iq_from_torque`](inc/damp/motor/acim.hpp#L131) | function | Motor control pack (if present) | q-axis current for a torque command at a given rotor flux (+1 more overload) |
| [`is_closed_loop_stable_discrete`](inc/damp/design/stability.hpp#L245) | function | Design-time synthesis (not PWM-rate) | Check closed-loop stability for discrete system with state feedback |
| [`is_controllable`](inc/damp/design/stability.hpp#L176) | function | Design-time synthesis (not PWM-rate) | Check if a system is controllable |
| [`is_fault`](inc/damp/toolbox/conditioning.hpp#L382) | function | Embedded helpers (controls-adjacent utilities) | True for a wire fault (FaultLow/FaultHigh) — i.e. not a real reading at all |
| [`is_matrix_element`](inc/damp/matrix/matrix_traits.hpp#L41) | block | Linear algebra | True if T may be a Matrix element type |
| [`is_observable`](inc/damp/design/stability.hpp#L194) | function | Design-time synthesis (not PWM-rate) | Check if a system is observable |
| [`is_salient`](inc/damp/motor/ipm.hpp#L62) | function | Motor control pack (if present) | True when \|Lq − Ld\| is meaningful saliency (same test as MTPA) |
| [`is_spherical_wrist`](inc/damp/kinematics/serial_arm.hpp#L564) | function | Kinematics / pose | Spherical-wrist (Pieper) criterion for a 6R chain: axes 4-5-6 intersect, i.e. `a₄ = a₅ = a₆ = 0` and `d₅ = 0` (within eps). `a₄` (joints[3].a) is required so the wrist centre sits at the joint-4/5/6 concurrency point |
| [`is_stabilizable`](inc/damp/design/riccati.hpp#L44) | function | Design-time synthesis (not PWM-rate) | Check if (A, B) is a stabilizable pair |
| [`is_stable_continuous`](inc/damp/analysis/poles.hpp#L58) | function | Frequency-domain analysis (host) | Check continuous-time stability |
| [`is_stable_discrete`](inc/damp/design/stability.hpp#L215) | function | Design-time synthesis (not PWM-rate) | Check if a discrete-time system matrix A is stable |
| [`is_usable`](inc/damp/toolbox/conditioning.hpp#L539) | function | Embedded helpers (controls-adjacent utilities) | True when a channel reading may still be used (in-span or saturated) |
| [`is_valid`](inc/damp/toolbox/conditioning.hpp#L377) | function | Embedded helpers (controls-adjacent utilities) | True only for the in-span status |
| [`isfinite`](inc/damp/math/math.hpp#L425) | function | Scalar math & complex | Finiteness test — false for NaN and ±∞ |
| [`iso_c_drive`](inc/damp/toolbox/io.hpp#L412) | function | Embedded helpers (controls-adjacent utilities) | ISO-C (single-stick, coordinated / curvature) drive |
| [`iso_s_drive`](inc/damp/toolbox/io.hpp#L393) | function | Embedded helpers (controls-adjacent utilities) | ISO-S (single-stick, speed-summed) drive — arcade travel stick |
| [`isstable`](inc/damp/matlab.hpp#L956) | function | MATLAB®-style aliases (host) | Continuous-time stability predicate on a state matrix (+2 more overloads) |
| [`jerk_stop_budget`](inc/damp/trajectory/online_otg.hpp#L125) | function | Trajectory value types | Velocity change while ramping residual accel to zero at jmax |
| [`jerk_toward_accel`](inc/damp/trajectory/online_otg.hpp#L196) | function | Trajectory value types | Jerk that slews a → a_des in one step, sat to ±jmax, no reverse past a_des |
| [`JerkLimitedAxis`](inc/damp/trajectory/online_otg.hpp#L486) | block | Trajectory value types | One online OTG axis with optional override jerks (MaxOv stack array) |
| [`JerkRequest`](inc/damp/trajectory/online_otg.hpp#L89) | block | Trajectory value types | One requester's proposed jerk for this tick (no private integrator) |
| [`Jet`](inc/damp/trajectory/trajectory_types.hpp#L95) | block | Trajectory value types | Truncated jet of a scalar motion sample: d[k] = s⁽ᵏ⁾ |
| [`joint_angles_from_frame_orientations`](inc/damp/kinematics/serial_arm.hpp#L525) | function | Kinematics / pose | Recover revolute joint angles from world orientations of DH frames {0,…,N} |
| [`JointLimits`](inc/damp/trajectory/cartesian_move.hpp#L88) | block | Trajectory value types | Per-joint velocity and acceleration limits for a task-space move |
| [`JointType`](inc/damp/kinematics/serial_arm.hpp#L74) | enum | Kinematics / pose | Joint actuation type for a DH joint |
| [`JordanBlock`](inc/damp/design/pole_placement.hpp#L831) | block | Design-time synthesis (not PWM-rate) | One Jordan mini-block of a desired closed-loop spectrum |
| [`JordanObjective`](inc/damp/design/pole_placement.hpp#L1165) | enum | Design-time synthesis (not PWM-rate) | Robustness objective for place_jordan_optimal (the paper's two methods) |
| [`JunctionEstimator`](inc/damp/motor/thermal.hpp#L220) | block | Motor control pack (if present) | FET junction-temperature estimator: case temperature plus a thermal model |
| [`kalman`](inc/damp/estimation/kalman.hpp#L103) | function | Observers & estimators | Steady-state Kalman filter design |
| [`KalmanFilter`](inc/damp/estimation/kalman.hpp#L249) | block | Observers & estimators | Runtime Kalman filter for embedded systems |
| [`KalmanResult`](inc/damp/estimation/kalman.hpp#L45) | block | Observers & estimators | Steady-state Kalman filter design result |
| [`lag`](inc/damp/controllers/lead_lag.hpp#L181) | function | Design-time synthesis (not PWM-rate) | Design a lag compensator from desired low-frequency gain boost |
| [`lambda_tuning`](inc/damp/design/pid_design.hpp#L398) | function | Design-time synthesis (not PWM-rate) | Lambda tuning for FOPDT model |
| [`Latch`](inc/damp/toolbox/logic.hpp#L62) | block | Embedded helpers (controls-adjacent utilities) | Set/reset latch. @tparam SetDominant which input wins when both are asserted (default: set-dominant, e.g. a trip overriding a clear) |
| [`lead`](inc/damp/controllers/lead_lag.hpp#L132) | function | Design-time synthesis (not PWM-rate) | Design a lead compensator from desired phase boost at a target frequency |
| [`lead_lag`](inc/damp/controllers/lead_lag.hpp#L234) | function | Design-time synthesis (not PWM-rate) | Design a lead-lag compensator (cascade of lead + lag sections) |
| [`lead_lag_direct`](inc/damp/controllers/lead_lag.hpp#L264) | function | Design-time synthesis (not PWM-rate) | Direct lead-lag specification from zero/pole locations |
| [`LeadLagController`](inc/damp/controllers/lead_lag.hpp#L281) | block | Runtime controllers | Discrete lead-lag compensator |
| [`LeadLagResult`](inc/damp/controllers/lead_lag.hpp#L59) | block | Design-time synthesis (not PWM-rate) | Lead-lag compensator design result |
| [`LeadLagSeriesResult`](inc/damp/controllers/lead_lag.hpp#L204) | block | Design-time synthesis (not PWM-rate) | Cascaded lead+lag design result (2nd-order StateSpace + success) |
| [`leakage_factor`](inc/damp/motor/acim.hpp#L72) | function | Motor control pack (if present) | Leakage factor σ = 1 − Lₘ²/(L_s L_r) |
| [`lerp`](inc/damp/toolbox/scaling.hpp#L35) | function | Embedded helpers (controls-adjacent utilities) | Linear interpolation between a and b by fraction t |
| [`LeverrierResult`](inc/damp/systems/zpk.hpp#L500) | block | LTI systems (SS / TF / ZPK / discretize) | Faddeev–LeVerrier characteristic polynomial and adjoint coefficient matrices |
| [`LIMIT`](inc/damp/toolbox/iec61131.hpp#L598) | function | Embedded helpers (controls-adjacent utilities) | LIMIT (IEC 61131-3 selection function): clamp in to [mn, mx] |
| [`linear_modulation_voltage`](inc/damp/motor/modulation.hpp#L130) | function | Motor control pack (if present) | Peak phase voltage at the linear SVPWM hexagon limit |
| [`linear_screw`](inc/damp/toolbox/actuator.hpp#L166) | function | Embedded helpers (controls-adjacent utilities) | Build a ServoAxis for a linear axis driven by a leadscrew/belt |
| [`LinearDelta`](inc/damp/kinematics/motion_maps.hpp#L285) | block | Kinematics / pose | Linear delta robot — per-carriage closed-form inverse, sphere- trilateration forward. Towers at 90°, 210°, 330° |
| [`LinearDeltaGeometry`](inc/damp/kinematics/motion_maps.hpp#L269) | block | Kinematics / pose | Linear delta geometry (three vertical carriages, fixed-length rods) |
| [`LinearizationResult`](inc/damp/design/linearization.hpp#L41) | block | Design-time synthesis (not PWM-rate) | Result of nonlinear operating-point linearization |
| [`linearize`](inc/damp/design/linearization.hpp#L135) | function | Design-time synthesis (not PWM-rate) | Linearize nonlinear dynamics and output maps about an operating point (+1 more overload) |
| [`LinearPath`](inc/damp/trajectory/cartesian_move.hpp#L62) | block | Trajectory value types | A straight-line path `p(s) = start + s·dir`, `s ∈ [0, length]` |
| [`LinkAttitudeEskf`](inc/damp/estimation/serial_arm_pose.hpp#L143) | block | Observers & estimators | One-frame attitude for undercarriage / house (accel-only ESKF) |
| [`linmod`](inc/damp/matlab.hpp#L329) | function | MATLAB®-style aliases (host) | MATLAB®-style nonlinear linearization about an operating point |
| [`locate_zero_crossing_exact`](inc/damp/simulation/hybrid.hpp#L335) | function | Simulation / SIL harness (host) | Locate a state zero-crossing of g(x) on an Exact LTI segment |
| [`log`](inc/damp/matrix/functions.hpp#L511) | function | Linear algebra | Principal matrix logarithm via inverse scaling and squaring |
| [`log`](inc/damp/math/math.hpp#L249) | function | Scalar math & complex | Natural logarithm |
| [`log10`](inc/damp/math/math.hpp#L380) | function | Scalar math & complex | Base-10 logarithm, log10(x) = ln(x) / ln(10) |
| [`logm`](inc/damp/matrix/functions.hpp#L550) | function | Linear algebra | @brief MATLAB®-style alias for log |
| [`LogPolicy`](inc/damp/simulation/hybrid.hpp#L58) | enum | Simulation / SIL harness (host) | How often to record samples along a hybrid trajectory |
| [`loop_metrics`](inc/damp/analysis/frequency.hpp#L480) | function | Frequency-domain analysis (host) | One-call loop analysis: compute L/S/T response and return compact metrics (+1 more overload) |
| [`loop_response`](inc/damp/analysis/frequency.hpp#L558) | function | Frequency-domain analysis (host) | Compute open-loop L, sensitivity S, complementary sensitivity T, and Nyquist data (+1 more overload) |
| [`LoopResponseResult`](inc/damp/analysis/frequency.hpp#L415) | block | Frequency-domain analysis (host) | Open-loop and closed-loop frequency response package |
| [`LoopSummary`](inc/damp/analysis/frequency.hpp#L448) | block | Frequency-domain analysis (host) | Compact loop summary metrics for quick stability/robustness checks |
| [`LowerTriangle`](inc/damp/matrix/views.hpp#L141) | block | Linear algebra | Lower triangular view of a square matrix |
| [`LowPass`](inc/damp/filters/lowpass.hpp#L49) | block | Filters & signal conditioning | Nth-order low-pass filter |
| [`lowpass_1st`](inc/damp/filters/iir_design.hpp#L136) | function | Filters & signal conditioning | First-order low-pass filter design (+1 more overload) |
| [`lowpass_2nd`](inc/damp/filters/iir_design.hpp#L221) | function | Filters & signal conditioning | Second-order low-pass filter design (+1 more overload) |
| [`lowpass_2nd_continuous`](inc/damp/filters/iir_design.hpp#L272) | function | Filters & signal conditioning | Second-order low-pass filter design (continuous-time) |
| [`lowshelf`](inc/damp/filters/iir_design.hpp#L733) | function | Filters & signal conditioning | Low-shelf EQ filter: boost or cut everything below fc |
| [`lqg`](inc/damp/matlab.hpp#L737) | function | MATLAB®-style aliases (host) | Linear-Quadratic-Gaussian regulator design |
| [`LQG`](inc/damp/controllers/lqg.hpp#L176) | block | Runtime controllers | Linear-Quadratic-Gaussian (LQG) controller |
| [`lqg_bundle`](inc/damp/design/synthesis.hpp#L322) | function | Design-time synthesis (not PWM-rate) | Synthesize the full LQG artifact bundle in one call |
| [`lqg_from_parts`](inc/damp/controllers/lqg.hpp#L149) | function | Design-time synthesis (not PWM-rate) | Assemble an LQG design from separately computed Kalman and LQR results |
| [`lqg_pr_bundle`](inc/damp/design/synthesis.hpp#L349) | function | Design-time synthesis (not PWM-rate) | Synthesize a SISO LQG + PR design with internal-model compensation |
| [`LQGAnalysisModels`](inc/damp/design/synthesis.hpp#L39) | block | Design-time synthesis (not PWM-rate) | Analysis-oriented models produced from an LQG design |
| [`LQGArtifacts`](inc/damp/design/synthesis.hpp#L89) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle: design + analysis models + runtime bundle |
| [`LQGI`](inc/damp/controllers/lqgi.hpp#L158) | block | Runtime controllers | Linear-Quadratic-Gaussian-Integral (LQGI) controller |
| [`lqgi_bundle`](inc/damp/design/synthesis.hpp#L396) | function | Design-time synthesis (not PWM-rate) | Synthesize the full LQGI servo artifact bundle in one call |
| [`LQGIAnalysisModels`](inc/damp/design/synthesis.hpp#L152) | block | Design-time synthesis (not PWM-rate) | Analysis-oriented models produced from an LQGI design |
| [`LQGIArtifacts`](inc/damp/design/synthesis.hpp#L180) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle: LQGI servo design + analysis + ready-to-run controller |
| [`LQGIResult`](inc/damp/controllers/lqgi.hpp#L34) | block | Design-time synthesis (not PWM-rate) | LQGI design result |
| [`LQGPRArtifacts`](inc/damp/design/synthesis.hpp#L133) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle for SISO LQG + PR |
| [`LQGPRRuntimeBundle`](inc/damp/design/synthesis.hpp#L100) | block | Design-time synthesis (not PWM-rate) | Runtime bundle for SISO LQG + PR internal model compensation |
| [`lqgreg`](inc/damp/matlab.hpp#L753) | function | MATLAB®-style aliases (host) | Combine separate Kalman filter and LQR designs into an LQG controller |
| [`LQGResult`](inc/damp/controllers/lqg.hpp#L46) | block | Design-time synthesis (not PWM-rate) | LQG design result |
| [`LQGRuntimeBundle`](inc/damp/design/synthesis.hpp#L50) | block | Design-time synthesis (not PWM-rate) | Runtime bundle for LQG control |
| [`lqgtrack`](inc/damp/matlab.hpp#L765) | function | MATLAB®-style aliases (host) | Linear-Quadratic-Gaussian design with integral action for tracking |
| [`lqi`](inc/damp/matlab.hpp#L724) | function | MATLAB®-style aliases (host) | Linear-Quadratic Integral design for tracking |
| [`LQI`](inc/damp/controllers/lqi.hpp#L155) | block | Runtime controllers | Linear-Quadratic-Integral (LQI) controller |
| [`lqi_bundle`](inc/damp/design/synthesis.hpp#L378) | function | Design-time synthesis (not PWM-rate) | Synthesize the full LQI servo artifact bundle in one call |
| [`LQIAnalysisModels`](inc/damp/design/synthesis.hpp#L143) | block | Design-time synthesis (not PWM-rate) | Analysis-oriented models produced from an LQI design |
| [`LQIArtifacts`](inc/damp/design/synthesis.hpp#L165) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle: LQI servo design + analysis + ready-to-run controller |
| [`LQIResult`](inc/damp/controllers/lqi.hpp#L48) | block | Design-time synthesis (not PWM-rate) | LQI design result |
| [`lqr`](inc/damp/matlab.hpp#L603) | function | MATLAB®-style aliases (host) | Continuous-time LQR design (MATLAB®'s lqr) |
| [`lqr_gain`](inc/damp/design/riccati.hpp#L714) | function | Design-time synthesis (not PWM-rate) | Optimal LQR state-feedback gain from a Riccati solution |
| [`LQRCost`](inc/damp/controllers/lqr.hpp#L142) | block | Runtime controllers | Discretized LQR cost weights (Q, R, N) for a sampled-data problem |
| [`lqrd`](inc/damp/controllers/lqr.hpp#L317) | function | Design-time synthesis (not PWM-rate) | Sampled-data LQR from continuous plant (MATLAB®-style short name) (+1 more overload) |
| [`lqrd`](inc/damp/matlab.hpp#L633) | function | MATLAB®-style aliases (host) | Design discrete LQR from continuous-time system via discretization (+1 more overload) |
| [`LQRResult`](inc/damp/controllers/lqr.hpp#L56) | block | Design-time synthesis (not PWM-rate) | Linear-Quadratic Regulator design result |
| [`lqry`](inc/damp/matlab.hpp#L680) | function | MATLAB®-style aliases (host) | Output-weighted continuous LQR (state cost Q = Cᵀ Q_y C) (+1 more overload) |
| [`lsim`](inc/damp/analysis/time_response.hpp#L257) | function | Frequency-domain analysis (host) | Forced time response of a (MIMO) state-space system to an input signal (+1 more overload) |
| [`LsimInfo`](inc/damp/analysis/time_response.hpp#L484) | block | Frequency-domain analysis (host) | Transient characteristics of an arbitrary response signal |
| [`lsiminfo`](inc/damp/analysis/time_response.hpp#L503) | function | Frequency-domain analysis (host) | Compute transient characteristics from an output/time signal |
| [`lsimplot`](inc/damp/simulation/plot_plotly.hpp#L749) | function | Simulation / SIL harness (host) | Plot a forced (lsim) simulation, one trace per output |
| [`LsimResult`](inc/damp/analysis/time_response.hpp#L66) | block | Frequency-domain analysis (host) | Result of a single-trajectory simulation: time, output, and state history |
| [`lu_decomposition`](inc/damp/matrix/decomposition.hpp#L131) | function | Linear algebra | LU decomposition with partial pivoting |
| [`lu_solve`](inc/damp/matrix/solve.hpp#L187) | function | Linear algebra | Solve AX = B via LU with partial pivoting |
| [`luenberger`](inc/damp/estimation/luenberger.hpp#L138) | function | Observers & estimators | Design a Luenberger observer by robust pole placement (matrix form) (+3 more overloads) |
| [`Luenberger`](inc/damp/estimation/luenberger.hpp#L400) | block | Observers & estimators | Luenberger state observer (runtime) |
| [`LuenbergerResult`](inc/damp/estimation/luenberger.hpp#L92) | block | Observers & estimators | Luenberger observer design result |
| [`Lut1D`](inc/damp/toolbox/lookup.hpp#L108) | block | Embedded helpers (controls-adjacent utilities) | 1-D interpolating lookup table over monotonic breakpoints |
| [`Lut2D`](inc/damp/toolbox/lookup.hpp#L200) | block | Embedded helpers (controls-adjacent utilities) | 2-D bilinear interpolating lookup table over a regular grid |
| [`Lut3D`](inc/damp/toolbox/lookup.hpp#L274) | block | Embedded helpers (controls-adjacent utilities) | 3-D trilinear interpolating lookup table over a regular grid |
| [`lut_segment`](inc/damp/toolbox/lookup.hpp#L63) | function | Embedded helpers (controls-adjacent utilities) | Index of the interpolation segment containing x |
| [`lyap`](inc/damp/design/lyapunov.hpp#L105) | function | Design-time synthesis (not PWM-rate) | Solve the continuous-time Lyapunov equation A X + X Aᵀ + Q = 0 |
| [`lyap`](inc/damp/matlab.hpp#L1052) | function | MATLAB®-style aliases (host) | MATLAB® alias for the continuous Lyapunov solve AX+XAᵀ+Q=0 |
| [`MadgwickFilter`](inc/damp/estimation/sensor_fusion.hpp#L105) | block | Observers & estimators | Madgwick gradient-descent AHRS filter |
| [`mag2db`](inc/damp/math/math.hpp#L459) | function | Scalar math & complex | Magnitude to decibels, 20·log10(mag) |
| [`MahonyFilter`](inc/damp/estimation/sensor_fusion.hpp#L163) | block | Observers & estimators | Mahony nonlinear complementary filter with PI correction |
| [`make1DOF`](inc/damp/matlab.hpp#L220) | function | MATLAB®-style aliases (host) | Force 1-DOF setpoint weights on a PID design result (b=c=1) |
| [`make2DOF`](inc/damp/matlab.hpp#L236) | function | MATLAB®-style aliases (host) | Apply 2-DOF setpoint weights on a PID design result |
| [`make_acim_servo`](inc/damp/motor/acim.hpp#L491) | function | Motor control pack (if present) | Build the common `Servo<Cascade\<T\>, Acim\<T\>>` from one flat config |
| [`make_brushed_servo`](inc/damp/motor/brushed_dc.hpp#L164) | function | Motor control pack (if present) | Build the common `Servo<Cascade\<T\>, BrushedDc\<T\>>` from one flat config |
| [`make_cartesian_move`](inc/damp/trajectory/cartesian_move.hpp#L266) | function | Trajectory value types | Deduction-friendly factory (deduces NJoints from joint_limits, callables from args) |
| [`make_ipm_servo`](inc/damp/motor/servo.hpp#L214) | function | Motor control pack (if present) | Build `Servo<Cascade, Pmsm<T, FieldWeakening>>` for the IPM daily path |
| [`make_pmsm_servo`](inc/damp/motor/servo.hpp#L190) | function | Motor control pack (if present) | Build the common `Servo<Cascade\<T\>, Pmsm<T, Policy>>` from one flat config |
| [`make_stepper_servo`](inc/damp/motor/stepper.hpp#L170) | function | Motor control pack (if present) | Build the common `Servo<Cascade\<T\>, Stepper\<T\>>` from one flat config |
| [`make_topp_move`](inc/damp/trajectory/topp.hpp#L361) | function | Trajectory value types | Deduction-friendly factory (deduces NJoints from joint_limits; pick NGrid) |
| [`margin`](inc/damp/matlab.hpp#L889) | function | MATLAB®-style aliases (host) | Gain and phase margins of a SISO loop over a frequency grid |
| [`MarginResult`](inc/damp/matlab.hpp#L868) | block | MATLAB®-style aliases (host) | Gain/phase margins and their crossover frequencies |
| [`margins_from_frf`](inc/damp/estimation/frequency_response.hpp#L771) | function | Observers & estimators | Gain / phase margins from a fixed on-target FRF table |
| [`Matrix`](inc/damp/matrix/core.hpp#L74) | block | Linear algebra | Fixed-size, stack-allocated matrix for linear algebra operations |
| [`MatrixLike`](inc/damp/matrix/matrix_traits.hpp#L70) | concept | Linear algebra | Concept for any type that provides 2D matrix-like element access |
| [`MatrixLikeOf`](inc/damp/matrix/matrix_traits.hpp#L80) | concept | Linear algebra | Concept for a MatrixLike type with specific dimensions |
| [`max`](inc/damp/matrix/core.hpp#L823) | function | Linear algebra | Addition of two MatrixLike types (with broadcasting support) (+1 more overload) |
| [`max_accel_for_velocity_error`](inc/damp/trajectory/online_otg.hpp#L139) | function | Trajectory value types | Max \|a\| still unloadable under jmax for remaining velocity error |
| [`max_iq`](inc/damp/motor/spm.hpp#L61) | function | Motor control pack (if present) | Max positive i_q for an SPM (i_d = 0) at electrical speed |
| [`max_torque_at_speed`](inc/damp/motor/ipm.hpp#L124) | function | Motor control pack (if present) | Max \|Te\| at a mechanical speed under V and I limits (MTPA + open-loop FW) (+1 more overload) |
| [`MeasJacobian`](inc/damp/estimation/ekf.hpp#L69) | block | Observers & estimators | Measurement prediction result from the user's observation function |
| [`mecanum_drive`](inc/damp/toolbox/io.hpp#L514) | function | Embedded helpers (controls-adjacent utilities) | Mecanum (H-layout, 45° rollers) mixer — alias for holonomic_drive |
| [`mech_turns_from_elec`](inc/damp/motor/sensorless.hpp#L78) | function | Motor control pack (if present) | Mechanical single-turn angle [turns] from electrical angle [rad] |
| [`mech_turns_per_s_from_elec`](inc/damp/motor/sensorless.hpp#L88) | function | Motor control pack (if present) | Mechanical speed [turns/s] from electrical velocity [elec rad/s] |
| [`MechanicalEstimator`](inc/damp/motor/mechanical_estimator.hpp#L128) | block | Motor control pack (if present) | Cheap-predict mechanical estimator for position, speed, and load torque |
| [`MechanicalEstimatorConfig`](inc/damp/motor/mechanical_estimator.hpp#L92) | block | Motor control pack (if present) | Configuration for MechanicalEstimator |
| [`mechanize_step_from_corrected`](inc/damp/estimation/ins_mechanization.hpp#L173) | function | Observers & estimators | One dead-reckoning step from an IMU sample (strapdown mechanization) |
| [`MedianFilter`](inc/damp/filters/median.hpp#L34) | block | Filters & signal conditioning | Sliding-window median filter — nonlinear spike/outlier rejection |
| [`mhe`](inc/damp/estimation/mhe.hpp#L166) | function | Observers & estimators | Synthesize a constrained moving-horizon estimator |
| [`MHE`](inc/damp/estimation/mhe.hpp#L283) | block | Observers & estimators | Runtime moving-horizon estimator (fixed per-tick iteration budget) (+1 more overload) |
| [`MHEArtifacts`](inc/damp/estimation/mhe.hpp#L110) | block | Observers & estimators | Window QP data produced by mhe(), consumed by damp::MHE |
| [`MHEConstraints`](inc/damp/estimation/mhe.hpp#L87) | block | Observers & estimators | Box constraints on the estimated states (applied to every window state) |
| [`middlebrook`](inc/damp/analysis/frequency.hpp#L1083) | function | Frequency-domain analysis (host) | Middlebrook stability analysis for cascaded source-load systems (+1 more overload) |
| [`MiddlebrookResult`](inc/damp/analysis/frequency.hpp#L896) | block | Frequency-domain analysis (host) | Result of Middlebrook minor loop gain analysis |
| [`min_accel`](inc/damp/trajectory/polynomial.hpp#L213) | function | Trajectory value types | Minimum-acceleration (cubic) rest-to-rest move p0 → pT over duration T |
| [`min_jerk`](inc/damp/trajectory/polynomial.hpp#L207) | function | Trajectory value types | Minimum-jerk (quintic) rest-to-rest move p0 → pT over duration T (Flash–Hogan) |
| [`min_snap`](inc/damp/trajectory/polynomial.hpp#L219) | function | Trajectory value types | Minimum-snap (septic) rest-to-rest move p0 → pT over duration T (Mellinger–Kumar) |
| [`minimal_realization`](inc/damp/design/minreal.hpp#L428) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for minreal |
| [`MinimalRealizationResult`](inc/damp/design/minreal.hpp#L84) | block | Design-time synthesis (not PWM-rate) | Minimal realization result |
| [`minmax`](inc/damp/backend.hpp#L164) | function | Core, configuration & backend vocabulary | Ordered {min, max} pair returned by value (+1 more overload) |
| [`minreal`](inc/damp/design/minreal.hpp#L284) | function | Design-time synthesis (not PWM-rate) | Minimal realization — cancel uncontrollable and unobservable modes |
| [`minreal`](inc/damp/matlab.hpp#L250) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::minreal |
| [`minreal_zpk`](inc/damp/systems/zpk.hpp#L867) | function | LTI systems (SS / TF / ZPK / discretize) | Cancel matching pole-zero pairs on a ZPK model |
| [`ModeExtractorConfig`](inc/damp/estimation/frequency_response.hpp#L411) | block | Observers & estimators | Configuration for FRF peak / valley extraction |
| [`ModeExtractorResult`](inc/damp/estimation/frequency_response.hpp#L428) | block | Observers & estimators | Result of reducing an FRF table to at most MaxModes peaks (resonances) |
| [`ModelReductionMethod`](inc/damp/design/model_reduction.hpp#L58) | enum | Design-time synthesis (not PWM-rate) | Method for eliminating states in modred |
| [`modified_sine_rise`](inc/damp/trajectory/cam.hpp#L327) | function | Trajectory value types | Modified-sine rise (Norton): low peak accel, common industrial law |
| [`modified_trapezoid_rise`](inc/damp/trajectory/cam.hpp#L368) | function | Trajectory value types | Modified-trapezoid rise (constant-accel flanks + harmonic transitions) |
| [`modred`](inc/damp/design/model_reduction.hpp#L518) | function | Design-time synthesis (not PWM-rate) | Model reduction by truncation or DC-matched residualization |
| [`modulation_duty_cycles`](inc/damp/motor/modulation.hpp#L399) | function | Motor control pack (if present) | Carrier-based VSI duty cycles from an αβ voltage command |
| [`modulation_index`](inc/damp/motor/modulation.hpp#L168) | function | Motor control pack (if present) | Modulation index relative to the linear SVPWM circle |
| [`modulation_zero_sequence`](inc/damp/motor/modulation.hpp#L326) | function | Motor control pack (if present) | Zero-sequence for any PwmScheme |
| [`Modulator`](inc/damp/motor/modulation.hpp#L496) | block | Motor control pack (if present) | Thin scheme-holding wrapper for FOC / deploy paths |
| [`Motor`](inc/damp/motor/drive.hpp#L64) | concept | Motor control pack (if present) | The pluggable electrical stage of a drive for one machine type |
| [`motor_constant`](inc/damp/motor/foc.hpp#L206) | function | Motor control pack (if present) | Motor constant Kₘ (torque per √copper-loss) — a figure of merit |
| [`MovingAverage`](inc/damp/filters/moving_average.hpp#L62) | block | Filters & signal conditioning | Moving-average (boxcar) filter — also a DC-preserving harmonic-notch comb |
| [`MPC`](inc/damp/controllers/mpc.hpp#L907) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its artifacts: MPC controller{art}; (+1 more overload) |
| [`mpc`](inc/damp/controllers/offset_free_mpc.hpp#L132) | function | Design-time synthesis (not PWM-rate) | Synthesize an offset-free constrained MPC (controller + estimator) |
| [`MPC`](inc/damp/controllers/mpc.hpp#L770) | block | Runtime controllers | Runtime constrained MPC controller (fixed per-tick iteration budget) |
| [`MPCAnalysisModels`](inc/damp/controllers/mpc.hpp#L666) | block | Runtime controllers | LTI models of the unconstrained MPC loop, for margin/robustness analysis |
| [`MPCArtifacts`](inc/damp/controllers/mpc.hpp#L235) | block | Runtime controllers | Condensed-QP data produced by mpc(), consumed by damp::MPC |
| [`MPCConstraints`](inc/damp/controllers/mpc.hpp#L197) | block | Runtime controllers | Box constraints for mpc() |
| [`MPCHorizonSuggestion`](inc/damp/controllers/mpc.hpp#L561) | block | Runtime controllers | Advisory NP/NC horizon values from suggest_mpc_horizon() |
| [`MPCWeights`](inc/damp/controllers/mpc.hpp#L162) | block | Runtime controllers | Cost weights for mpc() |
| [`mstogi`](inc/damp/filters/sogi.hpp#L117) | function | Filters & signal conditioning | Mixed Second/Third-Order Generalized Integrator (MSTOGI) (+1 more overload) |
| [`MSTOGI`](inc/damp/filters/sogi.hpp#L242) | block | Filters & signal conditioning | Runtime MSTOGI with exact resonator and forward-Euler washout |
| [`mtpa_id_from_iq`](inc/damp/motor/mtpa.hpp#L49) | function | Motor control pack (if present) | MTPA d-axis current on the trajectory for a given q-axis current |
| [`mtpa_reference`](inc/damp/motor/mtpa.hpp#L85) | function | Motor control pack (if present) | MTPA dq current reference for a commanded torque |
| [`MtpaReference`](inc/damp/motor/mtpa.hpp#L127) | block | Motor control pack (if present) | Maximum-torque-per-ampere current-reference generator (PMSM / IPMSM / SynRM) |
| [`multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L133) | function | Observers & estimators | Build a multi-sine design payload from a configuration |
| [`MultiPRController`](inc/damp/controllers/pr.hpp#L369) | block | Runtime controllers | Multi-harmonic PR Controller |
| [`MultiRateConfig`](inc/damp/simulation/multirate.hpp#L81) | block | Simulation / SIL harness (host) | Configuration for multi-rate closed-loop runs |
| [`MultiSine`](inc/damp/estimation/excitation/multi_sine.hpp#L200) | block | Observers & estimators | Sum-of-tones multi-sine runtime generator |
| [`MultiSineConfig`](inc/damp/estimation/excitation/multi_sine.hpp#L62) | block | Observers & estimators | Configuration for fixed-component multi-sine excitation |
| [`MultiSineResult`](inc/damp/estimation/excitation/multi_sine.hpp#L104) | block | Observers & estimators | Multi-sine design payload |
| [`MUX`](inc/damp/toolbox/iec61131.hpp#L617) | function | Embedded helpers (controls-adjacent utilities) | MUX (IEC 61131-3 multiplexer): select input k of N (0-based) |
| [`nameplate_design`](inc/damp/motor/drive_design.hpp#L232) | function | Motor control pack (if present) | Nameplate → MTPA limits + foc_cascade_tune + servo + field-weakening configs |
| [`nameplate_limits`](inc/damp/motor/ipm.hpp#L301) | function | Motor control pack (if present) | Voltage/current/base-speed/MTPA peak limits from bus and machine nameplate |
| [`NameplateDesignResult`](inc/damp/motor/acim.hpp#L527) | block | Motor control pack (if present) | Nameplate → cascade tune + filled motor::AcimServoConfig (Design Is Deploy entry) (+2 more overloads) |
| [`NameplateLimits`](inc/damp/motor/acim.hpp#L192) | block | Motor control pack (if present) | ACIM nameplate limits (base speed, peak torque, rated power) (+2 more overloads) |
| [`NavFrame`](inc/damp/estimation/ins_mechanization.hpp#L77) | enum | Observers & estimators | Local-level navigation frame for strapdown mechanization |
| [`nearbyint`](inc/damp/math/math.hpp#L345) | function | Scalar math & complex | Round to nearest integer; ties to even (IEEE default / `FE_TONEAREST`) |
| [`ned_from_enu`](inc/damp/estimation/ins_mechanization.hpp#L106) | function | Observers & estimators | Map an ENU vector into NED (axis permute) |
| [`negative_sequence_ab`](inc/damp/filters/pll.hpp#L244) | function | Filters & signal conditioning | Instantaneous negative-sequence αβ from a quadrature signal pair |
| [`nest_a_des`](inc/damp/trajectory/online_otg.hpp#L249) | function | Trajectory value types | Nest a_des under hard a limits and the jmax unload budget for rem |
| [`nichols`](inc/damp/analysis/frequency.hpp#L683) | function | Frequency-domain analysis (host) | Build Nichols points from existing Bode data (+2 more overloads) |
| [`nicholsplot`](inc/damp/simulation/plot_plotly.hpp#L936) | function | Simulation / SIL harness (host) | Plot a Nichols chart (open-loop phase vs magnitude) with M-circle grid |
| [`NicholsPoint`](inc/damp/analysis/frequency.hpp#L654) | block | Frequency-domain analysis (host) | Single-point Nichols chart sample (open-loop phase vs magnitude) |
| [`NicholsResult`](inc/damp/analysis/frequency.hpp#L666) | block | Frequency-domain analysis (host) | Nichols chart locus across a frequency sweep |
| [`NlmsFilter`](inc/damp/filters/fir.hpp#L546) | block | Filters & signal conditioning | Normalized LMS adaptive FIR filter |
| [`NoFieldWeakening`](inc/damp/motor/field_weakening.hpp#L220) | block | Motor control pack (if present) | Null field-weakening policy — passes the base reference through unchanged |
| [`norm`](inc/damp/matlab.hpp#L1122) | function | MATLAB®-style aliases (host) | MATLAB® alias for the H2 system norm norm(sys,2) |
| [`norm_h2`](inc/damp/analysis/norms.hpp#L76) | function | Frequency-domain analysis (host) | H2 norm of a state-space system |
| [`norm_hinf`](inc/damp/analysis/norms.hpp#L134) | function | Frequency-domain analysis (host) | H∞ norm of a state-space system: sup_ω σ̄(G(jω)) |
| [`notch`](inc/damp/filters/iir_design.hpp#L639) | function | Filters & signal conditioning | Second-order band-reject (notch) filter |
| [`NpcMap`](inc/damp/power/npc.hpp#L38) | block | Power electronics pack (if present) | Single-phase diode NPC leg map — voltage command → 4 duties (layer 3) |
| [`null`](inc/damp/matlab.hpp#L459) | function | MATLAB®-style aliases (host) | MATLAB® short alias for an orthonormal null-space basis |
| [`null_space`](inc/damp/matrix/svd.hpp#L350) | function | Linear algebra | Orthonormal basis for the null space {x : A·x = 0} via SVD |
| [`NullSpace`](inc/damp/matrix/svd.hpp#L331) | block | Linear algebra | Orthonormal basis for the null space (kernel) of a matrix |
| [`nyquist`](inc/damp/analysis/frequency.hpp#L497) | function | Frequency-domain analysis (host) | Compute Nyquist data for a SISO state-space system (+1 more overload) |
| [`nyquistplot`](inc/damp/simulation/plot_plotly.hpp#L821) | function | Simulation / SIL harness (host) | Plot a Nyquist locus with the -1 critical point marked |
| [`NyquistPoint`](inc/damp/analysis/frequency.hpp#L369) | block | Frequency-domain analysis (host) | Single-point Nyquist response data |
| [`NyquistResult`](inc/damp/analysis/frequency.hpp#L381) | block | Frequency-domain analysis (host) | Nyquist response data across a frequency sweep |
| [`observability_gramian`](inc/damp/design/stability.hpp#L146) | function | Design-time synthesis (not PWM-rate) | Continuous/discrete observability Gramian W_o |
| [`observability_matrix`](inc/damp/design/stability.hpp#L82) | function | Design-time synthesis (not PWM-rate) | Compute the observability matrix [C; CA; CA²; ...; CA^(N-1)] |
| [`obsv`](inc/damp/matlab.hpp#L274) | function | MATLAB®-style aliases (host) | MATLAB® short alias for observability_matrix (+1 more overload) |
| [`OffDelayTimer`](inc/damp/toolbox/logic.hpp#L133) | block | Embedded helpers (controls-adjacent utilities) | Off-delay timer: output goes true immediately when in is true and stays true until in has been false continuously for delay |
| [`OffsetFreeMPC`](inc/damp/controllers/offset_free_mpc.hpp#L258) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its artifacts: OffsetFreeMPC controller{art}; |
| [`OffsetFreeMPC`](inc/damp/controllers/offset_free_mpc.hpp#L196) | block | Runtime controllers | Runtime offset-free MPC: constrained MPC + disturbance-augmented Kalman filter |
| [`OffsetFreeMPCArtifacts`](inc/damp/controllers/offset_free_mpc.hpp#L84) | block | Runtime controllers | Combined MPC + disturbance-augmented Kalman design, consumed by damp::OffsetFreeMPC |
| [`omniwheel_drive`](inc/damp/toolbox/io.hpp#L520) | function | Embedded helpers (controls-adjacent utilities) | Omni-wheel (X-layout) mixer — alias for holonomic_drive |
| [`OnDelayTimer`](inc/damp/toolbox/logic.hpp#L96) | block | Embedded helpers (controls-adjacent utilities) | On-delay timer: output goes true once in has been held true continuously for delay; drops immediately when in goes false |
| [`one_norm`](inc/damp/matrix/functions.hpp#L64) | function | Linear algebra | One-norm ‖A‖₁: maximum absolute column sum |
| [`operator_a_des`](inc/damp/trajectory/online_otg.hpp#L300) | function | Trajectory value types | PIO a_des (legacy name used by drive governor docs) |
| [`OptimalJordanPlacement`](inc/damp/design/pole_placement.hpp#L1172) | block | Design-time synthesis (not PWM-rate) | Result of optimized arbitrary pole placement (place_jordan_optimal) |
| [`otg_jerk`](inc/damp/trajectory/online_otg.hpp#L315) | function | Trajectory value types | OTG jerk from any a_des generator, with terminal unload toward value_target |
| [`OutputFeedbackController`](inc/damp/concepts.hpp#L88) | concept | Core, configuration & backend vocabulary | Vector output-feedback controller: u = control(r, y), self-contained tick |
| [`pade`](inc/damp/matlab.hpp#L998) | function | MATLAB®-style aliases (host) | First-order Padé approximation of pure delay e^{−sT} (+1 more overload) |
| [`pade2`](inc/damp/matlab.hpp#L1020) | function | MATLAB®-style aliases (host) | Second-order Padé approximation of pure delay e^{−sT} (+1 more overload) |
| [`pade_delay_1st`](inc/damp/filters/iir_design.hpp#L378) | function | Filters & signal conditioning | First-order Pade approximation of time delay (+1 more overload) |
| [`pade_delay_2nd`](inc/damp/filters/iir_design.hpp#L422) | function | Filters & signal conditioning | Second-order Pade approximation of time delay (+1 more overload) |
| [`panel_color`](inc/damp/simulation/plot_plotly.hpp#L314) | function | Simulation / SIL harness (host) | Plotly default colorway entry for index i within a single subplot |
| [`panel_legend`](inc/damp/simulation/plot_plotly.hpp#L297) | function | Simulation / SIL harness (host) | Per-panel legend box just right of a stacked subplot band |
| [`panel_legend_id`](inc/damp/simulation/plot_plotly.hpp#L283) | function | Simulation / SIL harness (host) | Legend id for subplot row (1-based): "legend", "legend2", … |
| [`panel_x_domain`](inc/damp/simulation/plot_plotly.hpp#L273) | function | Simulation / SIL harness (host) | Horizontal domain for stacked multi-panel x-axes |
| [`parallel`](inc/damp/systems/state_space.hpp#L253) | function | LTI systems (SS / TF / ZPK / discretize) | Parallel connection (shared input, summed outputs) |
| [`ParameterDriftMonitor`](inc/damp/estimation/parameter_estimation.hpp#L252) | block | Observers & estimators | Residual-based drift detector: has the plant moved away from the model? |
| [`ParameterEstimator`](inc/damp/concepts.hpp#L142) | concept | Core, configuration & backend vocabulary | Online grey-box parameter estimator: physical parameters + gating |
| [`park_transform`](inc/damp/transforms.hpp#L390) | function | Motor control pack (if present) | Park transform (αβ → dq) |
| [`park_zero_transform`](inc/damp/transforms.hpp#L494) | function | Motor control pack (if present) | Park transform with zero passthrough (αβ0 → dq0) |
| [`peaking`](inc/damp/filters/iir_design.hpp#L708) | function | Filters & signal conditioning | Peaking (bell) EQ filter: boost or cut a band around f0 |
| [`Periodic`](inc/damp/toolbox/timing.hpp#L125) | block | Embedded helpers (controls-adjacent utilities) | Periodic trigger — fires once per elapsed period |
| [`periodic_cubic_spline`](inc/damp/trajectory/spline.hpp#L363) | function | Trajectory value types | Periodic cubic (C² across the cycle wrap) |
| [`periodic_quintic_spline`](inc/damp/trajectory/spline.hpp#L371) | function | Trajectory value types | Periodic quintic (C⁴ across the cycle wrap) |
| [`periodic_spline`](inc/damp/trajectory/spline.hpp#L279) | function | Trajectory value types | Periodic multi-waypoint spline (cyclic wrap of derivatives 1…Order−1) |
| [`phase_margin_from_damping_ratio`](inc/damp/design/pid_design.hpp#L540) | function | Design-time synthesis (not PWM-rate) | Approximate phase margin from damping ratio |
| [`phase_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L166) | function | Frequency-domain analysis (host) | Find phase margin using unwrapped phase trajectory |
| [`PhaseCalibrationCommand`](inc/damp/motor/calibration.hpp#L54) | block | Motor control pack (if present) | One step's output from PhaseParameterCalibrator |
| [`PhaseCalibrationConfig`](inc/damp/motor/calibration.hpp#L34) | block | Motor control pack (if present) | Configuration for online phase resistance/inductance commissioning |
| [`PhaseDrive`](inc/damp/motor/commutation.hpp#L143) | enum | Motor control pack (if present) | Per-phase drive state for one half-bridge under six-step commutation |
| [`PhaseParameterCalibrator`](inc/damp/motor/calibration.hpp#L111) | block | Motor control pack (if present) | Online phase R/L identification by recursive least squares (PRBS injected) |
| [`pi_pole_placement_first_order`](inc/damp/design/pid_design.hpp#L769) | function | Design-time synthesis (not PWM-rate) | PI gains that place the closed-loop poles of a first-order plant (+1 more overload) |
| [`pid`](inc/damp/controllers/pid.hpp#L222) | function | Design-time synthesis (not PWM-rate) | 2-DOF continuous PID controller design |
| [`pid`](inc/damp/matlab.hpp#L161) | function | MATLAB®-style aliases (host) | MATLAB®-style parallel-form continuous PID constructor |
| [`pid_from_bandwidth`](inc/damp/design/pid_design.hpp#L432) | function | Design-time synthesis (not PWM-rate) | Design PID from desired bandwidth and phase margin |
| [`pid_from_performance_spec`](inc/damp/design/pid_design.hpp#L603) | function | Design-time synthesis (not PWM-rate) | Design PID directly from settling-time and overshoot targets |
| [`pid_pole_placement`](inc/damp/design/pid_design.hpp#L640) | function | Design-time synthesis (not PWM-rate) | Direct PID pole placement for a first-order-plus-dead-time model (+1 more overload) |
| [`PIDController`](inc/damp/controllers/pid.hpp#L312) | block | Runtime controllers | Fixed-rate discrete 2-DOF PID (canonical runtime) (+2 more overloads) |
| [`PIDMode`](inc/damp/controllers/pid.hpp#L247) | enum | Runtime controllers | Compile-time selection of the PID control-law structure |
| [`PIDPerformanceSpec`](inc/damp/design/pid_design.hpp#L585) | block | Design-time synthesis (not PWM-rate) | Time-domain performance targets for quick PID synthesis |
| [`PIDResult`](inc/damp/controllers/pid.hpp#L97) | block | Design-time synthesis (not PWM-rate) | 2-DOF continuous-time PID controller design result |
| [`PIDRuntimeMode`](inc/damp/controllers/pid.hpp#L271) | enum | Runtime controllers | Runtime operating mode for PIDController / ContinuousPID |
| [`pidstd`](inc/damp/matlab.hpp#L189) | function | MATLAB®-style aliases (host) | Standard-form continuous PID constructor (1-DOF) |
| [`pidstd2`](inc/damp/matlab.hpp#L207) | function | MATLAB®-style aliases (host) | Standard-form continuous 2-DOF PID constructor |
| [`pidtune`](inc/damp/matlab.hpp#L787) | function | MATLAB®-style aliases (host) | PID controller tuning using frequency domain method |
| [`PIDType`](inc/damp/design/pid_design.hpp#L47) | enum | Design-time synthesis (not PWM-rate) | PID controller type selection for tuning methods |
| [`PiecewiseLTI`](inc/damp/simulation/hybrid.hpp#L94) | block | Simulation / SIL harness (host) | Plant: NModes continuous LTI pieces (same state/input/output dimensions) |
| [`pinv`](inc/damp/matlab.hpp#L446) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the Moore–Penrose pseudoinverse |
| [`pio_a_des`](inc/damp/trajectory/online_otg.hpp#L269) | function | Trajectory value types | PIO a_des generator: a_des = α rem + nest + rest min |
| [`place`](inc/damp/design/pole_placement.hpp#L88) | function | Design-time synthesis (not PWM-rate) | Robust multi-input pole placement (Kautsky–Nichols–Van Dooren, real poles) (+5 more overloads) |
| [`place`](inc/damp/matlab.hpp#L582) | function | MATLAB®-style aliases (host) | Robust multi-input pole placement (MATLAB®'s place) |
| [`place_discrete`](inc/damp/design/pole_placement.hpp#L581) | function | Design-time synthesis (not PWM-rate) | Discrete plant + z-plane poles → discrete \(K\) (+3 more overloads) |
| [`place_jordan`](inc/damp/design/pole_placement.hpp#L1138) | function | Design-time synthesis (not PWM-rate) | Exact pole placement with an arbitrary Jordan structure (Schmid–Ntogramatzidis–Nguyen–Pandey / Klein–Moore parametric form) |
| [`place_jordan_optimal`](inc/damp/design/pole_placement.hpp#L1237) | function | Design-time synthesis (not PWM-rate) | Robust / minimum-gain arbitrary pole placement (Schmid et al., Methods 1–2) |
| [`place_observer`](inc/damp/design/pole_placement.hpp#L659) | function | Design-time synthesis (not PWM-rate) | Continuous observer gain: place eigenvalues of \(A - LC\) (s-plane poles) (+5 more overloads) |
| [`place_observer_discrete`](inc/damp/design/pole_placement.hpp#L769) | function | Design-time synthesis (not PWM-rate) | Discrete plant + z-plane poles → discrete observer gain \(L\) (+3 more overloads) |
| [`plot_bode`](inc/damp/simulation/plot_plotly.hpp#L422) | function | Simulation / SIL harness (host) | Plot Bode magnitude and phase as subplots |
| [`plot_line`](inc/damp/simulation/plot_plotly.hpp#L491) | function | Simulation / SIL harness (host) | Simple line plot of time vs value |
| [`plot_simulation`](inc/damp/simulation/plot_plotly.hpp#L343) | function | Simulation / SIL harness (host) | Plot simulation results with subplots for states, outputs, and inputs |
| [`plot_step`](inc/damp/simulation/plot_plotly.hpp#L621) | function | Simulation / SIL harness (host) | Plot step response data |
| [`plot_xy`](inc/damp/simulation/plot_plotly.hpp#L556) | function | Simulation / SIL harness (host) | Multi-series planar scatter (true 2-D path / locus) (+1 more overload) |
| [`plots_root_for`](inc/damp/simulation/plot_plotly.hpp#L118) | function | Simulation / SIL harness (host) | Walk up from dir until a directory named "plots" is found |
| [`Pmsm`](inc/damp/motor/pmsm.hpp#L71) | block | Motor control pack (if present) | Field-oriented PMSM/BLDC electrical stage |
| [`pmsm_config`](inc/damp/motor/servo.hpp#L169) | function | Motor control pack (if present) | Split a PmsmServoConfig into a PmsmConfig for the motor stage |
| [`pmsm_nameplate_design`](inc/damp/motor/drive_design.hpp#L134) | function | Motor control pack (if present) | Nameplate → limits + foc_cascade_tune + filled motor::PmsmServoConfig |
| [`PmsmConfig`](inc/damp/motor/pmsm.hpp#L45) | block | Motor control pack (if present) | Configuration for a Pmsm electrical stage |
| [`PmsmNameplateDesignResult`](inc/damp/motor/drive_design.hpp#L99) | block | Motor control pack (if present) | Full nameplate design result: limits + cascade tune + deployable config |
| [`PmsmServoConfig`](inc/damp/motor/servo.hpp#L105) | block | Motor control pack (if present) | Flat configuration for a PMSM cascade servo (make_pmsm_servo) |
| [`PolarMap`](inc/damp/kinematics/motion_maps.hpp#L95) | block | Kinematics / pose | Polar / R-θ mapping (radius + angle ↔ Cartesian X/Y) |
| [`pole`](inc/damp/matlab.hpp#L853) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the open-loop poles of a system |
| [`PoleInfo`](inc/damp/analysis/poles.hpp#L75) | block | Frequency-domain analysis (host) | Natural frequency and damping ratio for each pole |
| [`poles`](inc/damp/analysis/poles.hpp#L41) | function | Frequency-domain analysis (host) | Compute open-loop poles (eigenvalues of A matrix) |
| [`PoleZeroMap`](inc/damp/analysis/poles.hpp#L115) | block | Frequency-domain analysis (host) | Poles and zeros of a system, for pole-zero plotting |
| [`poly_horner`](inc/damp/toolbox/scaling.hpp#L134) | function | Embedded helpers (controls-adjacent utilities) | Evaluate a polynomial at x by Horner's method |
| [`poly_rise`](inc/damp/trajectory/polynomial.hpp#L242) | function | Trajectory value types | Rest–rest polynomial rise of height h over unit parameter u ∈ [0,1] |
| [`poly_roots`](inc/damp/analysis/poles.hpp#L128) | function | Frequency-domain analysis (host) | Roots of a polynomial given in ascending powers (MATLAB® `roots`, reversed order) |
| [`poly_roots`](inc/damp/systems/zpk.hpp#L118) | function | LTI systems (SS / TF / ZPK / discretize) | Compute roots of an ascending-power real polynomial |
| [`poly_trajectory`](inc/damp/trajectory/polynomial.hpp#L149) | function | Trajectory value types | Synthesize a fixed-duration polynomial matching endpoint derivative jets (+1 more overload) |
| [`PolynomialTrajectory`](inc/damp/trajectory/polynomial.hpp#L265) | block | Trajectory value types | Runtime evaluator for a precomputed polynomial trajectory |
| [`PolyRootsResult`](inc/damp/systems/zpk.hpp#L89) | block | LTI systems (SS / TF / ZPK / discretize) | Roots of a real polynomial given in ascending powers |
| [`PolyTrajectory`](inc/damp/trajectory/polynomial.hpp#L73) | block | Trajectory value types | A synthesized polynomial trajectory: the coefficients of p(t) = Σ cᵢ·tⁱ over t ∈ [0, T], plus the duration |
| [`Pose`](inc/damp/kinematics/pose.hpp#L60) | block | Kinematics / pose | Rigid-body pose: a translation and an orientation (unit quaternion) |
| [`PositionAxes`](inc/damp/estimation/serial_arm_pose.hpp#L57) | block | Observers & estimators | Which world Cartesian components of a frame origin to fuse |
| [`positive_sequence_ab`](inc/damp/filters/pll.hpp#L223) | function | Filters & signal conditioning | Instantaneous positive-sequence αβ from a quadrature signal pair |
| [`pow`](inc/damp/matrix/functions.hpp#L570) | function | Linear algebra | Integer matrix power via binary exponentiation (+1 more overload) |
| [`pow`](inc/damp/math/math.hpp#L268) | function | Scalar math & complex | Power function, base^exponent (+1 more overload) |
| [`pr`](inc/damp/controllers/pr.hpp#L153) | function | Design-time synthesis (not PWM-rate) | Design a Proportional-Resonant controller |
| [`pr_harmonics`](inc/damp/controllers/pr.hpp#L185) | function | Design-time synthesis (not PWM-rate) | Design multiple-harmonic PR controller gains |
| [`prbs`](inc/damp/estimation/excitation/prbs.hpp#L158) | function | Observers & estimators | Build a PRBS design payload from a configuration |
| [`PRBS`](inc/damp/estimation/excitation/prbs.hpp#L181) | block | Observers & estimators | Maximal-length PRBS runtime generator |
| [`PRBSConfig`](inc/damp/estimation/excitation/prbs.hpp#L81) | block | Observers & estimators | Configuration for maximal-length pseudo-random binary excitation |
| [`PRBSResult`](inc/damp/estimation/excitation/prbs.hpp#L125) | block | Observers & estimators | PRBS design payload |
| [`PRController`](inc/damp/controllers/pr.hpp#L222) | block | Runtime controllers | Discrete Proportional-Resonant Controller |
| [`predicted_abs_value`](inc/damp/trajectory/online_otg.hpp#L440) | function | Trajectory value types | Predict \|value\| after one step of jerk j (selector ranking) |
| [`project_affine`](inc/damp/controllers/action_governor.hpp#L231) | function | Design-time synthesis (not PWM-rate) | Project u_des onto the polyhedron A u ≤ b (Euclidean) |
| [`project_box`](inc/damp/controllers/action_governor.hpp#L78) | function | Design-time synthesis (not PWM-rate) | Euclidean projection of u onto the axis-aligned box [umin, umax] (+1 more overload) |
| [`ProportionalCurrentMap`](inc/damp/toolbox/io.hpp#L139) | block | Embedded helpers (controls-adjacent utilities) | Affine map from a normalized command to a proportional current (mA) |
| [`PRResult`](inc/damp/controllers/pr.hpp#L72) | block | Design-time synthesis (not PWM-rate) | Proportional-Resonant controller design result |
| [`pseudo_inverse`](inc/damp/matrix/svd.hpp#L303) | function | Linear algebra | Moore–Penrose pseudoinverse A⁺ via SVD |
| [`PulseTimer`](inc/damp/toolbox/logic.hpp#L170) | block | Embedded helpers (controls-adjacent utilities) | Pulse timer (non-retriggerable): a rising edge of in emits a fixed |
| [`pwm_center_aligned_ccr`](inc/damp/motor/inverter.hpp#L94) | function | Motor control pack (if present) | Single CCR for up-down centre-aligned PWM (ARR = peak count) |
| [`pwm_center_aligned_ccrs`](inc/damp/motor/inverter.hpp#L133) | function | Motor control pack (if present) | Three single-CCR levels for up-down centre-aligned PWM |
| [`pwm_center_aligned_compare`](inc/damp/motor/inverter.hpp#L69) | function | Motor control pack (if present) | Dual-edge centre-aligned compare levels (full carrier period in ticks) |
| [`pwm_center_aligned_compare_deadtime`](inc/damp/motor/inverter.hpp#L151) | function | Motor control pack (if present) | Dual-edge compares with delayed turn-on (software dead-time on rise) |
| [`pwm_center_aligned_compares`](inc/damp/motor/inverter.hpp#L113) | function | Motor control pack (if present) | Dual-edge compares for three legs (full-period tick count) |
| [`pwm_center_aligned_edges`](inc/damp/motor/inverter.hpp#L42) | function | Motor control pack (if present) | Centre-aligned ideal high window for duty d in period Ts |
| [`PwmScheme`](inc/damp/motor/modulation.hpp#L82) | enum | Motor control pack (if present) | Carrier-based three-phase VSI modulation scheme |
| [`pzmap`](inc/damp/analysis/poles.hpp#L154) | function | Frequency-domain analysis (host) | Pole-zero map of a SISO transfer function (MATLAB® `pzmap(tf)`) (+2 more overloads) |
| [`pzplot`](inc/damp/simulation/plot_plotly.hpp#L848) | function | Simulation / SIL harness (host) | Plot a pole-zero map on the complex plane (poles as ×, zeros as ○) |
| [`QPResult`](inc/damp/design/qp.hpp#L82) | block | Design-time synthesis (not PWM-rate) | Result of a dense QP solve |
| [`QPStatus`](inc/damp/design/qp.hpp#L51) | enum | Design-time synthesis (not PWM-rate) | Termination status of a QP solve |
| [`qr_decompose`](inc/damp/matrix/decomposition.hpp#L223) | function | Linear algebra | Thin QR via modified Gram–Schmidt |
| [`QRDecomposition`](inc/damp/matrix/decomposition.hpp#L191) | block | Linear algebra | Thin QR factorization A = QR (modified Gram–Schmidt) |
| [`QuadMode`](inc/damp/toolbox/encoder.hpp#L47) | enum | Embedded helpers (controls-adjacent utilities) | Quadrature decode resolution (edges counted per A/B cycle) |
| [`quadprog`](inc/damp/matlab.hpp#L1073) | function | MATLAB®-style aliases (host) | MATLAB® alias for the dense inequality-constrained QP solve (+1 more overload) |
| [`quadratic_form`](inc/damp/matrix/core.hpp#L924) | function | Linear algebra | Symmetric congruence (quadratic) form  S = M X Mᵀ |
| [`QuadratureDecoder`](inc/damp/toolbox/encoder.hpp#L69) | block | Embedded helpers (controls-adjacent utilities) | Software A/B quadrature decoder with optional index |
| [`Quaternion`](inc/damp/math/geometry.hpp#L381) | block | Scalar math & complex | Unit quaternion rotation (w, x, y, z) (Hamilton product) |
| [`quintic_rise`](inc/damp/trajectory/cam.hpp#L266) | function | Trajectory value types | Quintic rest–rest rise (degree 5; C² ends) |
| [`quintic_spline`](inc/damp/trajectory/spline.hpp#L251) | function | Trajectory value types | Quintic (C⁴ — jerk- and snap-continuous) spline; clamped end velocity + accel |
| [`R_TRIG`](inc/damp/toolbox/iec61131.hpp#L118) | block | Embedded helpers (controls-adjacent utilities) | R_TRIG (Rising Edge Trigger) |
| [`rad2deg`](inc/damp/math/math.hpp#L481) | function | Scalar math & complex | Radians to degrees, rad·180/π |
| [`ramp`](inc/damp/estimation/excitation/ramp.hpp#L94) | function | Observers & estimators | Build a ramp design payload from a configuration |
| [`Ramp`](inc/damp/estimation/excitation/ramp.hpp#L110) | block | Observers & estimators | Rate-limited ramp runtime generator |
| [`RampConfig`](inc/damp/estimation/excitation/ramp.hpp#L32) | block | Observers & estimators | Configuration for a slew-rate-limited ramp excitation |
| [`RampResult`](inc/damp/estimation/excitation/ramp.hpp#L64) | block | Observers & estimators | Ramp design payload |
| [`RangeMonitor`](inc/damp/toolbox/conditioning.hpp#L430) | block | Embedded helpers (controls-adjacent utilities) | Analog-input range/fault monitor (NAMUR NE43 pattern) |
| [`rank`](inc/damp/design/stability.hpp#L163) | function | Design-time synthesis (not PWM-rate) | Compute rank of a matrix via Gaussian elimination with partial pivoting |
| [`rank`](inc/damp/matrix/functions.hpp#L259) | function | Linear algebra | Matrix rank via Gaussian elimination with partial pivoting |
| [`rank_from_svd`](inc/damp/matrix/svd.hpp#L271) | function | Linear algebra | Numerical rank from a precomputed SVD result |
| [`reduced_luenberger`](inc/damp/estimation/luenberger.hpp#L271) | function | Observers & estimators | Design a reduced-order (Gopinath) observer by pole placement (matrix form) (+3 more overloads) |
| [`ReducedLuenberger`](inc/damp/estimation/luenberger.hpp#L507) | block | Observers & estimators | Reduced-order (Gopinath) state observer (runtime) |
| [`ReducedLuenbergerResult`](inc/damp/estimation/luenberger.hpp#L217) | block | Observers & estimators | Reduced-order (Gopinath) observer design result |
| [`RedundantAnalogPair`](inc/damp/toolbox/io.hpp#L90) | block | Embedded helpers (controls-adjacent utilities) | Dual-channel analog input with same-slope or opposite-slope cross-check |
| [`ReferredRL`](inc/damp/motor/sepex.hpp#L164) | block | Motor control pack (if present) | Field plant as seen from the rotary-transformer primary (referred) |
| [`reg`](inc/damp/matlab.hpp#L519) | function | MATLAB®-style aliases (host) | Form dynamic regulator from system, state-feedback gain, and estimator gain |
| [`relay_autotune`](inc/damp/estimation/relay_autotune.hpp#L213) | function | Observers & estimators | Build a validated relay-autotune design payload |
| [`RelayAutotuneConfig`](inc/damp/estimation/relay_autotune.hpp#L120) | block | Observers & estimators | Configuration for the relay-feedback autotuning experiment |
| [`RelayAutotuneOutput`](inc/damp/estimation/relay_autotune.hpp#L238) | block | Observers & estimators | Per-tick output of RelayAutotuner::step |
| [`RelayAutotuner`](inc/damp/estimation/relay_autotune.hpp#L257) | block | Observers & estimators | Runtime relay-feedback autotuner |
| [`RelayAutotuneResult`](inc/damp/estimation/relay_autotune.hpp#L180) | block | Observers & estimators | Relay-autotuner design payload |
| [`RelayAutotuneStatus`](inc/damp/estimation/relay_autotune.hpp#L222) | enum | Observers & estimators | Lifecycle state of a RelayAutotuner |
| [`repetitive`](inc/damp/controllers/repetitive.hpp#L171) | function | Design-time synthesis (not PWM-rate) | Synthesize a repetitive controller with a scalar robustness filter Q |
| [`repetitive_binomial`](inc/damp/controllers/repetitive.hpp#L223) | function | Design-time synthesis (not PWM-rate) | Synthesize a repetitive controller with a binomial zero-phase FIR Q |
| [`RepetitiveConfig`](inc/damp/controllers/repetitive.hpp#L70) | block | Runtime controllers | Repetitive-controller tuning + period (with optional zero-phase FIR Q) |
| [`RepetitiveController`](inc/damp/controllers/repetitive.hpp#L273) | block | Runtime controllers | Plug-in repetitive controller runtime (fixed-size internal model) |
| [`RepetitiveResult`](inc/damp/controllers/repetitive.hpp#L134) | block | Design-time synthesis (not PWM-rate) | Design result for the repetitive controller |
| [`rescale`](inc/damp/toolbox/scaling.hpp#L69) | function | Embedded helpers (controls-adjacent utilities) | Affine map of x from the input range to the output range |
| [`ResistiveLossModel`](inc/damp/motor/thermal.hpp#L196) | block | Motor control pack (if present) | Minimal conduction-only loss model for a weak datasheet |
| [`resonance_autotune_from_frf`](inc/damp/estimation/commissioning_frontends.hpp#L76) | function | Observers & estimators | Map FRF peaks and valleys to a fixed bank of biquad compensators |
| [`ResonanceAutotuneResult`](inc/damp/estimation/commissioning_frontends.hpp#L44) | block | Observers & estimators | Result of resonance / anti-resonance compensator design from an FRF |
| [`ResonantMode`](inc/damp/estimation/frequency_response.hpp#L384) | block | Observers & estimators | One extracted FRF extremum (peak or valley) |
| [`revolute_angle_about_axis`](inc/damp/kinematics/serial_arm.hpp#L491) | function | Kinematics / pose | Signed revolute angle about a known axis from two body→world orientations |
| [`RisingEdge`](inc/damp/toolbox/logic.hpp#L30) | block | Embedded helpers (controls-adjacent utilities) | Rising-edge detector: true on the tick x goes false → true |
| [`RK23`](inc/damp/simulation/integrator.hpp#L902) | block | Simulation / SIL harness (host) | Bogacki-Shampine 2(3) adaptive integrator |
| [`RK3`](inc/damp/simulation/integrator.hpp#L969) | block | Simulation / SIL harness (host) | Classical 3rd-order Runge-Kutta (RK3) integrator |
| [`RK4`](inc/damp/simulation/integrator.hpp#L683) | block | Simulation / SIL harness (host) | Classical 4th-order Runge-Kutta (RK4) integrator |
| [`rlocus`](inc/damp/analysis/poles.hpp#L262) | function | Frequency-domain analysis (host) | Root locus of a SISO state-space plant over an explicit gain grid (+1 more overload) |
| [`rlocusplot`](inc/damp/simulation/plot_plotly.hpp#L1033) | function | Simulation / SIL harness (host) | Plot a root locus (closed-loop poles vs gain) on the complex plane |
| [`rls`](inc/damp/estimation/rls.hpp#L154) | function | Observers & estimators | Build scalar RLS design payload |
| [`Rls`](inc/damp/estimation/rls.hpp#L193) | block | Observers & estimators | Scalar runtime RLS estimator |
| [`rls_vector`](inc/damp/estimation/rls.hpp#L170) | function | Observers & estimators | Build vector RLS design payload |
| [`RlsConfig`](inc/damp/estimation/rls.hpp#L40) | block | Observers & estimators | Common RLS configuration |
| [`RlsResult`](inc/damp/estimation/rls.hpp#L82) | block | Observers & estimators | Scalar RLS design payload |
| [`RlsState`](inc/damp/estimation/rls.hpp#L68) | block | Observers & estimators | Scalar RLS runtime state |
| [`RlsVector`](inc/damp/estimation/rls.hpp#L287) | block | Observers & estimators | Vector runtime RLS estimator (NP parameters) |
| [`RlsVectorResult`](inc/damp/estimation/rls.hpp#L125) | block | Observers & estimators | Vector RLS design payload for N parameters |
| [`RlsVectorState`](inc/damp/estimation/rls.hpp#L111) | block | Observers & estimators | Vector RLS runtime state for N parameters |
| [`RobustExactDifferentiator`](inc/damp/filters/differentiator.hpp#L63) | block | Filters & signal conditioning | First-order robust exact differentiator (super-twisting differentiator) |
| [`RootLocusResult`](inc/damp/analysis/poles.hpp#L203) | block | Frequency-domain analysis (host) | Root-locus data: closed-loop poles along a gain grid |
| [`rotary_gearbox`](inc/damp/toolbox/actuator.hpp#L149) | function | Embedded helpers (controls-adjacent utilities) | Build a ServoAxis for a rotary joint behind a gearbox |
| [`rotary_transformer_primary_voltage`](inc/damp/motor/sepex.hpp#L182) | function | Motor control pack (if present) | Primary voltage that produces secondary field voltage vf_sec under ideal RT |
| [`RotaryDelta`](inc/damp/kinematics/motion_maps.hpp#L161) | block | Kinematics / pose | Rotary delta robot — closed-form inverse, quadratic-intersection forward |
| [`RotaryDeltaGeometry`](inc/damp/kinematics/motion_maps.hpp#L140) | block | Kinematics / pose | Rotary delta geometry (three base servos, parallelogram arms) |
| [`RotaryTransformerField`](inc/damp/motor/sepex.hpp#L250) | block | Motor control pack (if present) | Rotary-transformer field path parameters (primary-referred control) |
| [`rotational_load_ss`](inc/damp/motor/mechanical_estimator.hpp#L54) | function | Motor control pack (if present) | Continuous state-space model of a 1-DOF rotational drivetrain with an augmented load-torque state |
| [`rotor_time_constant`](inc/damp/motor/acim.hpp#L90) | function | Motor control pack (if present) | Rotor time constant T_r = L_r / R_r [s] |
| [`rows`](inc/damp/estimation/ins_eskf.hpp#L348) | function | Observers & estimators | Sparse y = F x for the 15-state INS first-order structure (G = I path) |
| [`RowVec`](inc/damp/matrix/rowvec.hpp#L29) | block | Linear algebra | Row vector specialization of Matrix<1, N, T> |
| [`RowView`](inc/damp/matrix/views.hpp#L165) | block | Linear algebra | Non-owning row view of a matrix |
| [`RS`](inc/damp/toolbox/iec61131.hpp#L91) | block | Embedded helpers (controls-adjacent utilities) | RS Latch (Reset-Set Latch) |
| [`sample_rise_law`](inc/damp/trajectory/cam.hpp#L445) | function | Trajectory value types | Sample a unit-parameter rise law onto NPts knots on [0, rise_length] |
| [`scaled_deadband`](inc/damp/toolbox/conditioning.hpp#L233) | function | Embedded helpers (controls-adjacent utilities) | Center dead zone that rescales the surviving range back to full span |
| [`scara_arm`](inc/damp/kinematics/scara.hpp#L262) | function | Kinematics / pose | Build a series SCARA (RRPR) as a 4-joint DH chain |
| [`schroeder_multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L180) | function | Observers & estimators | Build a multi-sine design payload with Schroeder phases applied |
| [`scurve`](inc/damp/trajectory/scurve.hpp#L196) | function | Trajectory value types | Synthesize a minimum-time jerk-limited (7-segment double-S) profile from (Xi, Vi) to (Xf, Vf) under asymmetric kinematic limits (+1 more overload) |
| [`ScurveProfile`](inc/damp/trajectory/scurve.hpp#L72) | block | Trajectory value types | A synthesized jerk-limited (double-S) profile: a sequence of constant-jerk segments, evaluated exactly (cubic in t within a segment) |
| [`ScurveSegment`](inc/damp/trajectory/scurve.hpp#L56) | block | Trajectory value types | One constant-jerk segment, valid for t ∈ [t0, t0 + duration), with the position/velocity/acceleration cached at the segment start |
| [`ScurveTrajectory`](inc/damp/trajectory/scurve.hpp#L307) | block | Trajectory value types | Runtime evaluator for a precomputed jerk-limited (double-S) profile |
| [`SecondOrderCoeffs`](inc/damp/filters/iir_design.hpp#L65) | block | Filters & signal conditioning | DSP coefficients for second-order IIR filter |
| [`SEL`](inc/damp/toolbox/iec61131.hpp#L606) | function | Embedded helpers (controls-adjacent utilities) | SEL (IEC 61131-3 binary selection): g ? in1 : in0 |
| [`select_jerk`](inc/damp/trajectory/online_otg.hpp#L452) | function | Trajectory value types | Pick the jerk that most reduces \|value\| if any override beats the operator |
| [`select_nearest`](inc/damp/kinematics/serial_arm.hpp#L455) | function | Kinematics / pose | Pick the solution branch nearest a reference configuration |
| [`SensorlessEstimator`](inc/damp/motor/sensorless_flux.hpp#L82) | block | Motor control pack (if present) | Sensorless rotor flux/position estimator for a PMSM, with optional sensor fusion |
| [`SensorlessObserver`](inc/damp/motor/sensorless.hpp#L62) | concept | Motor control pack (if present) | Sensorless observer concept (electrical frame) |
| [`Sepex`](inc/damp/motor/sepex.hpp#L290) | block | Motor control pack (if present) | Separately-excited DC electrical stage (armature + field bridges) |
| [`SepexConfig`](inc/damp/motor/sepex.hpp#L262) | block | Motor control pack (if present) | Configuration for a Sepex electrical stage |
| [`SepexServoConfig`](inc/damp/motor/sepex.hpp#L468) | block | Motor control pack (if present) | Flat configuration for a SEPEX cascade servo (make_sepex_servo) |
| [`SepicMap`](inc/damp/power/sepic.hpp#L29) | block | Power electronics pack (if present) | SEPIC CCM plant map — HardSwitchingConverter (layer 3) |
| [`septic_rise`](inc/damp/trajectory/cam.hpp#L299) | function | Trajectory value types | Septic rest–rest rise (degree 7; zero end jerk; C³ into dwells) |
| [`SequenceComponents`](inc/damp/transforms.hpp#L623) | block | Motor control pack (if present) | Symmetrical (sequence) components of a three-phase phasor set |
| [`serial_arm`](inc/damp/kinematics/serial_arm.hpp#L578) | function | Kinematics / pose | Validate a serial-arm DH chain and flag a spherical wrist |
| [`serial_arm_pose_design`](inc/damp/estimation/serial_arm_pose.hpp#L272) | function | Observers & estimators | Diagonal process / prior for joint-space pose filter |
| [`SerialArm`](inc/damp/kinematics/serial_arm.hpp#L209) | block | Kinematics / pose | Serial N-DOF manipulator runtime (joint-space N free; task space SE(3)) |
| [`SerialArmConfig`](inc/damp/kinematics/serial_arm.hpp#L176) | block | Kinematics / pose | Validated serial-arm configuration (the design payload) |
| [`SerialArmPoseFilter`](inc/damp/estimation/serial_arm_pose.hpp#L353) | block | Observers & estimators | Joint-centric arm pose: optional base attitude + peel gyros → q + FK |
| [`SerialArmPoseResult`](inc/damp/estimation/serial_arm_pose.hpp#L225) | block | Observers & estimators | Design payload for the joint-space arm filter (+ optional attitude) |
| [`series`](inc/damp/systems/state_space.hpp#L196) | function | LTI systems (SS / TF / ZPK / discretize) | Series connection: sys2 follows sys1 (u → sys1 → sys2 → y) (+1 more overload) |
| [`Servo`](inc/damp/motor/servo.hpp#L52) | block | Motor control pack (if present) | A drive plus an outer torque law — convenience for the cascade daily path |
| [`ServoAxis`](inc/damp/toolbox/actuator.hpp#L107) | block | Embedded helpers (controls-adjacent utilities) | One servoactuator transmission: SI joint unit ⟷ drive (motor) units |
| [`ServoBank`](inc/damp/toolbox/actuator.hpp#L229) | block | Embedded helpers (controls-adjacent utilities) | A bank of ServoAxis transmissions: maps a synchronized multi-axis TrajectoryState array straight to drive commands in one call |
| [`ServoCommand`](inc/damp/toolbox/actuator.hpp#L84) | block | Embedded helpers (controls-adjacent utilities) | A drive-native servoactuator setpoint: position, velocity, torque |
| [`sgn`](inc/damp/math/math.hpp#L397) | function | Scalar math & complex | Sign function — −1 if val < 0, 1 if val > 0, 0 if val == 0 |
| [`shape_stick`](inc/damp/trajectory/selector_governor.hpp#L56) | function | Trajectory value types | Joystick / RC stick feel: scaled center dead zone then RC expo |
| [`ShaperType`](inc/damp/trajectory/input_shaper.hpp#L56) | enum | Trajectory value types | Input-shaper family |
| [`sigma`](inc/damp/analysis/frequency.hpp#L796) | function | Frequency-domain analysis (host) | Singular-value frequency response of a (possibly MIMO) state-space system (+1 more overload) |
| [`sigmaplot`](inc/damp/simulation/plot_plotly.hpp#L990) | function | Simulation / SIL harness (host) | Plot singular-value frequency response (log frequency, dB) |
| [`SigmaPoint`](inc/damp/analysis/frequency.hpp#L739) | block | Frequency-domain analysis (host) | Singular values of G(jω) at one frequency |
| [`SigmaResult`](inc/damp/analysis/frequency.hpp#L753) | block | Frequency-domain analysis (host) | Singular-value frequency response over a grid |
| [`SignalSource`](inc/damp/concepts.hpp#L127) | concept | Core, configuration & backend vocabulary | Self-clocked signal source: u = step(), finished when done() |
| [`SignalStatus`](inc/damp/toolbox/conditioning.hpp#L368) | enum | Embedded helpers (controls-adjacent utilities) | Classification of an analog input against its valid/fault bands |
| [`simc`](inc/damp/design/pid_design.hpp#L338) | function | Design-time synthesis (not PWM-rate) | SIMC (Skogestad Internal Model Control) tuning for FOPDT models |
| [`simulate`](inc/damp/simulation/simulate.hpp#L132) | function | Simulation / SIL harness (host) | Simulate a nonlinear plant with a controller in closed loop |
| [`simulate_cached_zoh`](inc/damp/simulation/cached_zoh.hpp#L148) | function | Simulation / SIL harness (host) | Simulate a continuous LTI plant with fixed ZOH step and held input (+1 more overload) |
| [`simulate_cached_zoh_multirate`](inc/damp/simulation/cached_zoh.hpp#L227) | function | Simulation / SIL harness (host) | Multi-rate ZOH: held control samples, dense plant samples for plotting |
| [`simulate_discrete`](inc/damp/simulation/simulate.hpp#L478) | function | Simulation / SIL harness (host) | Simulate a discrete-time system with a controller |
| [`simulate_discrete_nonlinear`](inc/damp/simulation/simulate.hpp#L429) | function | Simulation / SIL harness (host) | Simulate a discrete-time nonlinear plant with a controller |
| [`simulate_lti`](inc/damp/simulation/simulate.hpp#L319) | function | Simulation / SIL harness (host) | Simulate a continuous LTI system with a controller |
| [`simulate_lti_siso`](inc/damp/simulation/simulate.hpp#L377) | function | Simulation / SIL harness (host) | SIL: continuous SISO LTI plant under a sampled SISO controller with constant reference |
| [`simulate_multirate`](inc/damp/simulation/multirate.hpp#L128) | function | Simulation / SIL harness (host) | Multi-rate closed-loop simulation with MultiRateConfig (+1 more overload) |
| [`simulate_piecewise_lti`](inc/damp/simulation/hybrid.hpp#L182) | function | Simulation / SIL harness (host) | Simulate a piecewise-LTI plant with Exact integration between events |
| [`simulate_sampled`](inc/damp/simulation/simulate.hpp#L205) | function | Simulation / SIL harness (host) | Simulate a continuous plant under a discrete (sampled) controller — multi-rate |
| [`simulate_state_feedback`](inc/damp/simulation/simulate.hpp#L264) | function | Simulation / SIL harness (host) | Simulate a nonlinear plant with state-feedback controller |
| [`simulate_two_rate`](inc/damp/simulation/multirate.hpp#L244) | function | Simulation / SIL harness (host) | Two-rate cascade with MultiRateConfig (+1 more overload) |
| [`simulate_two_rate_lti`](inc/damp/simulation/multirate.hpp#L388) | function | Simulation / SIL harness (host) | Two-rate cascade on a continuous LTI plant (A, B, C) |
| [`SimulationResult`](inc/damp/simulation/simulate.hpp#L103) | block | Simulation / SIL harness (host) | Result of a closed-loop simulation |
| [`sin`](inc/damp/matrix/functions.hpp#L723) | function | Linear algebra | Matrix sine via scaling and double-angle reconstruction |
| [`sin`](inc/damp/math/math.hpp#L182) | function | Scalar math & complex | Sine |
| [`sincos`](inc/damp/matrix/functions.hpp#L653) | function | Linear algebra | Compute sin(A) and cos(A) together via scaling and double-angle reconstruction |
| [`sincos`](inc/damp/math/math.hpp#L200) | function | Scalar math & complex | Combined sine and cosine, {sin(x), cos(x)} |
| [`SinCosCalibration`](inc/damp/motor/resolver.hpp#L48) | block | Motor control pack (if present) | Per-channel gain/offset trims for a SinCosDecoder (hardware calibration) |
| [`SinCosDecoder`](inc/damp/motor/resolver.hpp#L78) | block | Motor control pack (if present) | Resolver / sin-cos-encoder decoder: (sin, cos) → angle, speed, multi-turn position |
| [`SinCosState`](inc/damp/motor/resolver.hpp#L57) | block | Motor control pack (if present) | Decoded state from one SinCosDecoder step |
| [`SinglePhasePLL`](inc/damp/filters/pll.hpp#L43) | block | Filters & signal conditioning | Single-Phase PLL |
| [`sinh`](inc/damp/matrix/functions.hpp#L793) | function | Linear algebra | Matrix hyperbolic sine sinh(A) = (exp(A) − exp(−A))/2 |
| [`siso_ref`](inc/damp/simulation/simulate.hpp#L90) | function | Simulation / SIL harness (host) | Build a SisoReferenceAdapter for a SISO `control(r,y)` controller |
| [`SisoController`](inc/damp/concepts.hpp#L74) | concept | Core, configuration & backend vocabulary | Scalar output-feedback controller: u = control(r, y), fixed rate |
| [`SisoReferenceAdapter`](inc/damp/simulation/simulate.hpp#L71) | block | Simulation / SIL harness (host) | Adapt a SISO `control(r, y)` controller for callables that want `u = f(y)` or `u = f(t, y)` (the older simulate contracts) |
| [`six_step_duty_cycles`](inc/damp/motor/modulation.hpp#L465) | function | Motor control pack (if present) | Classical six-step (full-wave) duty pattern from the αβ angle |
| [`six_step_fundamental_voltage`](inc/damp/motor/modulation.hpp#L149) | function | Motor control pack (if present) | Six-step (square-wave) fundamental peak phase voltage |
| [`SixStepCommutator`](inc/damp/motor/commutation.hpp#L171) | block | Motor control pack (if present) | Six-step (trapezoidal) commutator: electrical sector → phase energization |
| [`SixStepDrive`](inc/damp/motor/commutation.hpp#L150) | block | Motor control pack (if present) | One commutation step: the drive state of each of the three phases {a, b, c} |
| [`SlewLimiter`](inc/damp/toolbox/conditioning.hpp#L268) | block | Embedded helpers (controls-adjacent utilities) | Slew-rate limiter: bound how fast the output may follow the target |
| [`sliding_jerk_toward_value`](inc/damp/trajectory/online_otg.hpp#L347) | function | Trajectory value types | PIO operator path: pio_a_des then otg_jerk |
| [`slip_frequency`](inc/damp/motor/acim.hpp#L147) | function | Motor control pack (if present) | Slip frequency [elec rad/s] for IFOC |
| [`smc`](inc/damp/controllers/smc.hpp#L59) | function | Design-time synthesis (not PWM-rate) | Bundle hand-picked SMC parameters into an SMCResult |
| [`SMCController`](inc/damp/controllers/smc.hpp#L136) | block | Runtime controllers | First-order sliding-mode controller (SMC) for a SISO plant |
| [`SMCResult`](inc/damp/controllers/smc.hpp#L32) | block | Design-time synthesis (not PWM-rate) | Tuning parameters for a first-order sliding-mode controller |
| [`smith_predictor_from_fopdt`](inc/damp/controllers/smith_predictor.hpp#L141) | function | Design-time synthesis (not PWM-rate) | Design a Smith predictor from FOPDT parameters |
| [`smith_predictor_from_pid`](inc/damp/controllers/smith_predictor.hpp#L207) | function | Design-time synthesis (not PWM-rate) | Pack a Smith predictor from an existing continuous PID and FOPDT plant |
| [`SmithPredictor`](inc/damp/controllers/smith_predictor.hpp#L267) | block | Runtime controllers | Smith predictor runtime: SISO primary + FO model + pure delay |
| [`SmithPredictorResult`](inc/damp/controllers/smith_predictor.hpp#L85) | block | Design-time synthesis (not PWM-rate) | Smith predictor design result (discrete PID + FO model + delay samples) |
| [`smo`](inc/damp/motor/sensorless_smo.hpp#L67) | function | Motor control pack (if present) | Design SMO gains from SPM electrical parameters |
| [`SmoObserver`](inc/damp/motor/sensorless_smo.hpp#L97) | block | Motor control pack (if present) | αβ sliding-mode back-EMF sensorless observer |
| [`SmoResult`](inc/damp/motor/sensorless_smo.hpp#L33) | block | Motor control pack (if present) | SMO design payload |
| [`soft_sign`](inc/damp/toolbox/conditioning.hpp#L103) | function | Embedded helpers (controls-adjacent utilities) | Smooth unit direction x / √(x²+ε²) ∈ (−1,1) |
| [`sogi`](inc/damp/filters/sogi.hpp#L39) | function | Filters & signal conditioning | Second-Order Generalized Integrator (SOGI) design (+1 more overload) |
| [`SOGI`](inc/damp/filters/sogi.hpp#L177) | block | Filters & signal conditioning | Runtime SOGI wrapper around design::sogi(w0, alpha, Ts) |
| [`SogiFll`](inc/damp/filters/sogi.hpp#L330) | block | Filters & signal conditioning | SOGI with a Frequency-Locked Loop — self-tuning single-tone tracker |
| [`solve`](inc/damp/matrix/solve.hpp#L82) | function | Linear algebra | Solve lower-triangular system LX = B via forward substitution (+2 more overloads) |
| [`solve_miqp`](inc/damp/design/qp.hpp#L1215) | function | Design-time synthesis (not PWM-rate) | Mixed-integer QP by branch and bound over the active-set relaxation |
| [`solve_qp`](inc/damp/design/qp.hpp#L525) | function | Design-time synthesis (not PWM-rate) | Solve a strictly convex inequality-constrained QP (dual active-set) (+1 more overload) |
| [`solve_qp_admm`](inc/damp/design/qp.hpp#L916) | function | Design-time synthesis (not PWM-rate) | One-shot ADMM QP solve (cold start) |
| [`solve_qp_interior_point`](inc/damp/design/qp.hpp#L960) | function | Design-time synthesis (not PWM-rate) | Interior-point QP solver (primal-dual path following) |
| [`SolveResult`](inc/damp/simulation/solver.hpp#L59) | block | Simulation / SIL harness (host) | Result of an ODE solve operation |
| [`SOPDTModel`](inc/damp/analysis/identification.hpp#L89) | block | Frequency-domain analysis (host) | Second-order plus dead-time candidate model |
| [`specific_force_at_rest`](inc/damp/estimation/ins_mechanization.hpp#L212) | function | Observers & estimators | Specific force [m/s²] a stationary IMU measures for the given orientation |
| [`spline`](inc/damp/trajectory/spline.hpp#L155) | function | Trajectory value types | Synthesize a multi-waypoint spline through points at times |
| [`SplineProfile`](inc/damp/trajectory/spline.hpp#L69) | block | Trajectory value types | A synthesized multi-waypoint spline: per-segment polynomial coefficients (ascending power, in segment-local time) plus the knot times |
| [`SplineSurface`](inc/damp/toolbox/lookup.hpp#L393) | block | Embedded helpers (controls-adjacent utilities) | 2-D interpolating surface with smooth (Catmull-Rom spline) blending |
| [`SplineTrajectory`](inc/damp/trajectory/spline.hpp#L391) | block | Trajectory value types | Runtime player for a multi-waypoint spline (design::SplineProfile) |
| [`sqrt`](inc/damp/matrix/functions.hpp#L446) | function | Linear algebra | Matrix square root via Denman–Beavers iteration |
| [`sqrt`](inc/damp/math/complex.hpp#L331) | function | Scalar math & complex | Compute complex square root (constexpr) (+1 more overload) |
| [`sqrtm`](inc/damp/matrix/functions.hpp#L483) | function | Linear algebra | @brief MATLAB®-style alias for sqrt |
| [`SR`](inc/damp/toolbox/iec61131.hpp#L63) | block | Embedded helpers (controls-adjacent utilities) | SR Latch (Set-dominant Set-Reset Latch) |
| [`ss`](inc/damp/matlab.hpp#L97) | function | MATLAB®-style aliases (host) | MATLAB®-style state-space model constructor |
| [`ss2tf`](inc/damp/systems/zpk.hpp#L843) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for SISO to_transfer_function |
| [`ss2zpk`](inc/damp/systems/zpk.hpp#L813) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for SISO to_zpk (state-space) |
| [`stability_margin_continuous`](inc/damp/design/stability.hpp#L268) | function | Design-time synthesis (not PWM-rate) | Compute stability margin for continuous system |
| [`stability_margin_discrete`](inc/damp/design/stability.hpp#L297) | function | Design-time synthesis (not PWM-rate) | Compute stability margin for discrete system |
| [`state_mpc`](inc/damp/controllers/mpc.hpp#L502) | function | Design-time synthesis (not PWM-rate) | Synthesize a constrained linear MPC (condensed dense QP, Δu form) (+1 more overload) |
| [`StateEstimator`](inc/damp/concepts.hpp#L114) | concept | Core, configuration & backend vocabulary | State estimator: x̂ = estimate(y, u) — the fused per-tick form |
| [`StateFeedback`](inc/damp/controllers/lqr.hpp#L431) | block | Runtime controllers | Runtime full-state feedback law u = −Kx |
| [`StateFeedbackController`](inc/damp/concepts.hpp#L101) | concept | Core, configuration & backend vocabulary | State-feedback law: u = control(r, x), reference-first |
| [`StateJacobian`](inc/damp/estimation/ekf.hpp#L41) | block | Observers & estimators | State prediction result from the user's dynamics function |
| [`StateSpace`](inc/damp/systems/state_space.hpp#L125) | block | LTI systems (SS / TF / ZPK / discretize) | State-space representation for linear time-invariant systems (discrete or continuous) |
| [`StateSpaceZPKResult`](inc/damp/systems/zpk.hpp#L570) | block | LTI systems (SS / TF / ZPK / discretize) | SISO state-space → ZPK conversion result with runtime zero count |
| [`steady_state_vdq`](inc/damp/motor/pmsm_equations.hpp#L47) | function | Motor control pack (if present) | Steady-state dq voltage at an operating point (rotor-synchronous frame) |
| [`steady_state_voltage_magnitude`](inc/damp/motor/pmsm_equations.hpp#L64) | function | Motor control pack (if present) | \|Vdq\| magnitude at a steady-state operating point |
| [`SteadyStateKalmanFilter`](inc/damp/estimation/kalman.hpp#L426) | block | Observers & estimators | Steady-state (fixed-gain) Kalman estimator for LQG-class designs |
| [`steer_map`](inc/damp/toolbox/io.hpp#L452) | function | Embedded helpers (controls-adjacent utilities) | Curve-driven steering map — arbitrary `(throttle, turn)` geometry |
| [`steinhart_hart`](inc/damp/toolbox/thermistor.hpp#L118) | function | Embedded helpers (controls-adjacent utilities) | Fit the Steinhart-Hart coefficients from three calibration points |
| [`step`](inc/damp/analysis/time_response.hpp#L153) | function | Frequency-domain analysis (host) | Step response of a (MIMO) state-space system (+1 more overload) |
| [`step`](inc/damp/estimation/excitation/step.hpp#L91) | function | Observers & estimators | Build a single-step design payload |
| [`Step`](inc/damp/estimation/excitation/step.hpp#L108) | block | Observers & estimators | Single step with optional \|dy/dt\| settle detection |
| [`step_train`](inc/damp/estimation/excitation/step_train.hpp#L94) | function | Observers & estimators | Build a step-train design payload from a configuration |
| [`StepConfig`](inc/damp/estimation/excitation/step.hpp#L35) | block | Observers & estimators | Configuration for a single step with optional settle detection |
| [`StepInfo`](inc/damp/analysis/time_response.hpp#L354) | block | Frequency-domain analysis (host) | Step-response characteristics of a single output signal |
| [`stepinfo`](inc/damp/analysis/time_response.hpp#L374) | function | Frequency-domain analysis (host) | Compute step-response characteristics from an output/time signal (+1 more overload) |
| [`stepped_sine`](inc/damp/estimation/excitation/stepped_sine.hpp#L104) | function | Observers & estimators | Build a stepped-sine design payload from a configuration |
| [`SteppedSine`](inc/damp/estimation/excitation/stepped_sine.hpp#L125) | block | Observers & estimators | Stepped-sine excitation — one pure tone at a time across a frequency table |
| [`SteppedSineConfig`](inc/damp/estimation/excitation/stepped_sine.hpp#L33) | block | Observers & estimators | Configuration for single-tone stepped-sine excitation |
| [`SteppedSineResult`](inc/damp/estimation/excitation/stepped_sine.hpp#L73) | block | Observers & estimators | Stepped-sine design payload |
| [`Stepper`](inc/damp/motor/stepper.hpp#L71) | block | Motor control pack (if present) | Two-phase field-oriented stepper electrical stage |
| [`StepperConfig`](inc/damp/motor/stepper.hpp#L49) | block | Motor control pack (if present) | Configuration for a Stepper electrical stage |
| [`StepperServoConfig`](inc/damp/motor/stepper.hpp#L150) | block | Motor control pack (if present) | Flat configuration for a stepper cascade servo (make_stepper_servo) |
| [`stepplot`](inc/damp/simulation/plot_plotly.hpp#L717) | function | Simulation / SIL harness (host) | Plot a step response, one trace per input/output pair |
| [`StepResponseSummary`](inc/damp/analysis/identification.hpp#L36) | block | Frequency-domain analysis (host) | Compact statistics extracted from a step experiment |
| [`StepResult`](inc/damp/estimation/excitation/step.hpp#L67) | block | Observers & estimators | Single-step design payload |
| [`StepTrain`](inc/damp/estimation/excitation/step_train.hpp#L110) | block | Observers & estimators | Alternating +/- step train runtime generator |
| [`StepTrainConfig`](inc/damp/estimation/excitation/step_train.hpp#L32) | block | Observers & estimators | Configuration for alternating +/- step excitation |
| [`StepTrainResult`](inc/damp/estimation/excitation/step_train.hpp#L64) | block | Observers & estimators | Step-train design payload |
| [`stewart`](inc/damp/kinematics/stewart.hpp#L317) | function | Kinematics / pose | Validate a hand-entered Stewart geometry and confirm the home pose is reachable |
| [`stewart_symmetric`](inc/damp/kinematics/stewart.hpp#L352) | function | Kinematics / pose | Tier-2 builder for the common symmetric hexagonal layout |
| [`StewartConfig`](inc/damp/kinematics/stewart.hpp#L113) | block | Kinematics / pose | Validated Stewart configuration (the design payload) |
| [`StewartForward`](inc/damp/kinematics/stewart.hpp#L102) | block | Kinematics / pose | Result of a forward (Newton–Raphson) solve |
| [`StewartGeometry`](inc/damp/kinematics/stewart.hpp#L65) | block | Kinematics / pose | Rig geometry: the six fixed base anchors `bᵢ`, the six moving-platform anchors `pᵢ`, the actuator stroke limits, and the nominal home height |
| [`StewartInverse`](inc/damp/kinematics/stewart.hpp#L95) | block | Kinematics / pose | Result of an inverse solve: the six leg lengths + a stroke-window flag |
| [`StewartPlatform`](inc/damp/kinematics/stewart.hpp#L165) | block | Kinematics / pose | Gough–Stewart platform runtime — closed-form inverse, Newton forward |
| [`Stopwatch`](inc/damp/toolbox/timing.hpp#L43) | block | Embedded helpers (controls-adjacent utilities) | Free-running elapsed-time accumulator |
| [`stsmc`](inc/damp/controllers/stsmc.hpp#L87) | function | Design-time synthesis (not PWM-rate) | Synthesize super-twisting gains from a disturbance-derivative bound |
| [`stsmc_gains`](inc/damp/controllers/stsmc.hpp#L123) | function | Design-time synthesis (not PWM-rate) | Super-twisting controller from gains you specify directly |
| [`STSMCController`](inc/damp/controllers/stsmc.hpp#L179) | block | Runtime controllers | Super-twisting controller (second-order sliding mode) |
| [`STSMCResult`](inc/damp/controllers/stsmc.hpp#L30) | block | Design-time synthesis (not PWM-rate) | Super-twisting (second-order sliding-mode) controller design result |
| [`subtract`](inc/damp/systems/state_space.hpp#L394) | function | LTI systems (SS / TF / ZPK / discretize) | Differencing connection: outputs y = y₁ − y₂ |
| [`SuccessiveCompensatorCommissioner`](inc/damp/estimation/successive_compensator.hpp#L208) | block | Observers & estimators | Successive biquad compensator bank (notch and/or band-pass) |
| [`SuccessiveCompensatorConfig`](inc/damp/estimation/successive_compensator.hpp#L50) | block | Observers & estimators | Configuration for successive compensator commissioning |
| [`suggest_cascade_bandwidths`](inc/damp/motor/foc_autotune.hpp#L60) | function | Motor control pack (if present) | Nested cascade bandwidths from a current-loop target (rad/s) |
| [`suggest_mpc_horizon`](inc/damp/controllers/mpc.hpp#L585) | function | Design-time synthesis (not PWM-rate) | Suggest MPC horizons from an explicit settling-time target (+1 more overload) |
| [`summarize_loop_response`](inc/damp/analysis/frequency.hpp#L460) | function | Frequency-domain analysis (host) | Summarize loop_response() results into one compact metrics struct |
| [`svd`](inc/damp/matrix/svd.hpp#L246) | function | Linear algebra | Full singular value decomposition A = U·Σ·Vᴴ (one-sided Jacobi) |
| [`svd`](inc/damp/matlab.hpp#L437) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the singular value decomposition |
| [`SVDResult`](inc/damp/matrix/svd.hpp#L228) | block | Linear algebra | Result of a full singular value decomposition A = U·Σ·Vᴴ |
| [`svm_duty_cycles`](inc/damp/motor/modulation.hpp#L442) | function | Motor control pack (if present) | Space-vector PWM duty cycles from an αβ voltage command |
| [`SvmDuties`](inc/damp/motor/modulation.hpp#L106) | block | Motor control pack (if present) | Result of a duty-map: half-bridge duties plus an over-modulation flag |
| [`svpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L203) | function | Motor control pack (if present) | Min-max zero-sequence injection for space-vector PWM |
| [`SweepIdentificationResult`](inc/damp/analysis/identification.hpp#L127) | block | Frequency-domain analysis (host) | End-to-end output of a sweep-based identification pass |
| [`Switch`](inc/damp/toolbox/io.hpp#L282) | block | Embedded helpers (controls-adjacent utilities) | Debounced maintained switch (toggle/selector contact) with change flag |
| [`SwitchedController`](inc/damp/controllers/composition.hpp#L117) | block | Runtime controllers | Bumpless-ish switch between normal, experiment, and backup SISO laws |
| [`SwitchMode`](inc/damp/controllers/composition.hpp#L94) | enum | Runtime controllers | Which path owns the plant command |
| [`symmetrical_components`](inc/damp/transforms.hpp#L653) | function | Motor control pack (if present) | Forward symmetrical-component (Fortescue) transform (abc → 012) |
| [`SymplecticEuler`](inc/damp/simulation/integrator.hpp#L229) | block | Simulation / SIL harness (host) | Semi-implicit (symplectic) Euler for mechanical systems |
| [`synthetic_arm_imus`](inc/damp/estimation/serial_arm_pose.hpp#L879) | function | Observers & estimators | Build IMU samples for all frames (needs three configs); optional mounts |
| [`synthetic_link_imu`](inc/damp/estimation/serial_arm_pose.hpp#L813) | function | Observers & estimators | Synthetic IMU at a sense point (gyro + specific force) for SIL |
| [`Tachometer`](inc/damp/toolbox/encoder.hpp#L144) | block | Embedded helpers (controls-adjacent utilities) | Pulse-based speed (tachometer) with frequency/period crossover |
| [`tan`](inc/damp/math/math.hpp#L215) | function | Scalar math & complex | Tangent |
| [`terminal_unload_jerk`](inc/damp/trajectory/online_otg.hpp#L167) | function | Trajectory value types | Continuous land jerk toward (v , a → 0): j = −a²/(2 rem) |
| [`tf`](inc/damp/matlab.hpp#L61) | function | MATLAB®-style aliases (host) | MATLAB®-style transfer function constructor (+1 more overload) |
| [`tf2zpk`](inc/damp/systems/zpk.hpp#L471) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for to_zpk (transfer function) |
| [`TFF`](inc/damp/toolbox/iec61131.hpp#L517) | block | Embedded helpers (controls-adjacent utilities) | T Flip-Flop (toggle on rising edge) |
| [`thermal_kalman_plant`](inc/damp/toolbox/thermal.hpp#L195) | function | Embedded helpers (controls-adjacent utilities) | Prepare a discretized thermal plant for design::kalman |
| [`ThermalKalmanObserver`](inc/damp/toolbox/thermal.hpp#L248) | block | Embedded helpers (controls-adjacent utilities) | Runtime thermal-network Kalman observer |
| [`ThermalLimiter`](inc/damp/motor/thermal.hpp#L308) | block | Motor control pack (if present) | Derates the current command from a temperature (Tj for FETs, winding for the motor) |
| [`ThermalLimits`](inc/damp/motor/thermal.hpp#L272) | block | Motor control pack (if present) | A derating curve plus a hard fault threshold |
| [`ThermalLossModel`](inc/damp/motor/thermal.hpp#L142) | concept | Motor control pack (if present) | A loss model usable by JunctionEstimator |
| [`ThermalState`](inc/damp/motor/thermal.hpp#L292) | block | Motor control pack (if present) | State from a ThermalLimiter evaluation |
| [`Thermistor`](inc/damp/toolbox/thermistor.hpp#L166) | block | Embedded helpers (controls-adjacent utilities) | NTC thermistor linearization (resistance → temperature) |
| [`ThermistorCoeffs`](inc/damp/toolbox/thermistor.hpp#L41) | block | Embedded helpers (controls-adjacent utilities) | Fitted NTC coefficients in Steinhart-Hart form |
| [`thipwm_zero_sequence`](inc/damp/motor/modulation.hpp#L229) | function | Motor control pack (if present) | Third-harmonic injection (THIPWM) zero-sequence — 1/6 of fundamental |
| [`three_level_modulator_1ph`](inc/damp/power/multilevel.hpp#L102) | function | Power electronics pack (if present) | 1φ unit modulator from pole-to-mid voltage command + NP balance |
| [`three_level_np_zero_sequence`](inc/damp/power/multilevel.hpp#L73) | function | Power electronics pack (if present) | Neutral-point zero-sequence offset on unit modulator m |
| [`three_level_pd_dwells`](inc/damp/power/multilevel.hpp#L44) | function | Power electronics pack (if present) | Phase-disposition (PD) dwell map for one pole |
| [`three_level_pole_to_mid`](inc/damp/power/multilevel.hpp#L128) | function | Power electronics pack (if present) | Pole-to-mid: balanced bus → P=+Vdc/2, O=0, N=−Vdc/2 |
| [`three_level_pole_to_neg`](inc/damp/power/multilevel.hpp#L120) | function | Power electronics pack (if present) | Average pole-to-negative voltage: N=0, O=v_mid, P=Vdc |
| [`ThreeLevelDwells`](inc/damp/power/multilevel.hpp#L34) | block | Power electronics pack (if present) | P/O/N dwell fractions for one pole (sum to 1 in the linear region) |
| [`ThreePhasePLL`](inc/damp/filters/pll.hpp#L141) | block | Filters & signal conditioning | Synchronous-reference-frame (SRF) PLL for balanced three-phase input |
| [`Timeout`](inc/damp/toolbox/timing.hpp#L78) | block | Embedded helpers (controls-adjacent utilities) | One-shot timeout |
| [`TimeResponse`](inc/damp/analysis/time_response.hpp#L51) | block | Frequency-domain analysis (host) | Multi-channel time-domain response sampled on a time grid |
| [`to_coeffs`](inc/damp/filters/iir_design.hpp#L467) | function | Filters & signal conditioning | Convert StateSpace system to first-order DSP coefficients (+3 more overloads) |
| [`to_double_vector`](inc/damp/simulation/plot_plotly.hpp#L236) | function | Simulation / SIL harness (host) | Convert std::vector\<T\> to std::vector<double> |
| [`to_std_vector`](inc/damp/simulation/plot_plotly.hpp#L211) | function | Simulation / SIL harness (host) | Convert a ColVec<N,T> to std::vector<double> for plotlypp |
| [`to_transfer_function`](inc/damp/systems/zpk.hpp#L826) | function | LTI systems (SS / TF / ZPK / discretize) | SISO state-space → transfer function (Leverrier), size NX+1 / NX+1 |
| [`to_zpk`](inc/damp/systems/zpk.hpp#L429) | function | LTI systems (SS / TF / ZPK / discretize) | Convert a SISO transfer function to zero-pole-gain form (+1 more overload) |
| [`TOF`](inc/damp/toolbox/iec61131.hpp#L226) | block | Embedded helpers (controls-adjacent utilities) | TOF Timer (Timer Off Delay) |
| [`Toggle`](inc/damp/toolbox/logic.hpp#L245) | block | Embedded helpers (controls-adjacent utilities) | Toggle (T flip-flop): output flips on each rising edge of in |
| [`TON`](inc/damp/toolbox/iec61131.hpp#L183) | block | Embedded helpers (controls-adjacent utilities) | TON Timer (Timer On Delay) |
| [`Tone`](inc/damp/estimation/excitation/multi_sine.hpp#L29) | block | Observers & estimators | One sinusoidal component in a multi-sine excitation |
| [`ToppMove`](inc/damp/trajectory/topp.hpp#L123) | block | Trajectory value types | Time-optimal task-space move (path-preserving, pointwise minimum-time) |
| [`ToppProfile`](inc/damp/trajectory/topp.hpp#L66) | block | Trajectory value types | The scalar time-optimal path-timing produced by TOPP |
| [`torque`](inc/damp/motor/acim.hpp#L120) | function | Motor control pack (if present) | Electromagnetic torque from rotor flux and q-axis current (+3 more overloads) |
| [`torque_constant`](inc/damp/motor/acim.hpp#L111) | function | Motor control pack (if present) | Torque constant at a rotor flux (amplitude-invariant Clarke/Park) (+1 more overload) |
| [`torque_constant_from_flux`](inc/damp/motor/foc.hpp#L103) | function | Motor control pack (if present) | Torque constant Kₜ of a PMSM (amplitude-invariant convention) |
| [`torque_constant_from_Kv`](inc/damp/motor/foc.hpp#L165) | function | Motor control pack (if present) | Torque constant from the datasheet velocity constant Kᵥ |
| [`torque_speed_envelope`](inc/damp/motor/ipm.hpp#L208) | function | Motor control pack (if present) | IPM torque–speed envelope table (host / Lut1D bake) (+1 more overload) |
| [`TorqueSpeedPoint`](inc/damp/motor/spm.hpp#L94) | block | Motor control pack (if present) | One sample of the SPM torque–speed envelope (MTPA i_d=0 branch) |
| [`TotemPoleMap`](inc/damp/power/totem_pole.hpp#L56) | block | Power electronics pack (if present) | Totem-pole bridgeless map — boost plant + slow-leg polarity fan-out |
| [`TP`](inc/damp/toolbox/iec61131.hpp#L271) | block | Embedded helpers (controls-adjacent utilities) | TP Timer (Timer Pulse) |
| [`TrajectoryBank`](inc/damp/trajectory/polynomial.hpp#L334) | block | Trajectory value types | Multi-axis coordination: time-scale each axis's profile to the slowest so a multi-DOF move starts and finishes synchronized ("linear" / coordinated joint moves — the feedforward reference for a manipulator) |
| [`TrajectoryBoundary`](inc/damp/trajectory/trajectory_types.hpp#L178) | block | Trajectory value types | Boundary conditions at one endpoint of a polynomial trajectory: a position and its time derivatives through jerk |
| [`TrajectoryLimits`](inc/damp/trajectory/trajectory_types.hpp#L45) | block | Trajectory value types | Asymmetric kinematic limits for a trapezoidal or S-curve motion profile |
| [`TrajectoryState`](inc/damp/trajectory/trajectory_types.hpp#L76) | block | Trajectory value types | A point on a motion profile: commanded position, velocity, acceleration |
| [`TransferFunction`](inc/damp/systems/transfer_function.hpp#L74) | block | LTI systems (SS / TF / ZPK / discretize) | SISO polynomial transfer function G(s) = num(s)/den(s) |
| [`Transform4`](inc/damp/math/geometry.hpp#L826) | block | Scalar math & complex | 4×4 homogeneous transform (SE(3) interop / DH export) |
| [`transient_inductance`](inc/damp/motor/acim.hpp#L81) | function | Motor control pack (if present) | Stator transient inductance L_σ = σ L_s [H] |
| [`Translation3`](inc/damp/kinematics/pose.hpp#L37) | block | Kinematics / pose | A 3-D translation — a thin Vec3 with domain-named conveniences |
| [`TransposeView`](inc/damp/matrix/views.hpp#L344) | block | Linear algebra | Non-owning transpose view of a matrix (zero-copy) |
| [`Trapezoidal`](inc/damp/simulation/integrator.hpp#L623) | block | Simulation / SIL harness (host) | Trapezoidal (Tustin) integrator |
| [`trapezoidal`](inc/damp/trajectory/trapezoidal.hpp#L206) | function | Trajectory value types | Synthesize the minimum-time asymmetric trapezoidal profile from (Xi, Vi) to (Xf, Vf) under the given limits (+1 more overload) |
| [`TrapezoidalProfile`](inc/damp/trajectory/trapezoidal.hpp#L56) | block | Trajectory value types | Planned trapezoidal profile: the segment durations, reached values, and boundary state needed to evaluate the trajectory at any time |
| [`TrapezoidalTrajectory`](inc/damp/trajectory/trapezoidal.hpp#L288) | block | Trajectory value types | Runtime evaluator for a precomputed trapezoidal profile |
| [`TRBDF2`](inc/damp/simulation/integrator.hpp#L497) | block | Simulation / SIL harness (host) | TR-BDF2 composite integrator — the stiff adaptive pair (ode23tb) |
| [`truncate_residual`](inc/damp/trajectory/online_otg.hpp#L231) | function | Trajectory value types | Kill residuals smaller than q (float dust → 0); leave larger values alone |
| [`TtypeMap`](inc/damp/power/ttype.hpp#L42) | block | Power electronics pack (if present) | Single-phase T-type (TNPC) leg map — voltage command → 4 duties (layer 3) |
| [`two_mass_rotational_ss`](inc/damp/motor/two_mass_estimator.hpp#L98) | function | Motor control pack (if present) | Continuous two-mass rotational drivetrain (compliant shaft) |
| [`two_mass_ss`](inc/damp/motor/two_mass_estimator.hpp#L135) | function | Motor control pack (if present) | Alias of two_mass_rotational_ss (roadmap short name) |
| [`two_mass_torsional_wn`](inc/damp/motor/two_mass_estimator.hpp#L145) | function | Motor control pack (if present) | Undamped two-mass torsional natural frequency [rad/s] |
| [`two_norm`](inc/damp/matrix/functions.hpp#L128) | function | Linear algebra | Spectral norm ‖A‖₂ = σₘₐₓ(A) |
| [`two_point_cal`](inc/damp/toolbox/scaling.hpp#L111) | function | Embedded helpers (controls-adjacent utilities) | Fit an AffineCal through two `(raw, engineering)` points |
| [`TwoMassEstimator`](inc/damp/motor/two_mass_estimator.hpp#L209) | block | Motor control pack (if present) | Dual-encoder two-mass mechanical estimator (shaft compliance) |
| [`TwoMassEstimatorConfig`](inc/damp/motor/two_mass_estimator.hpp#L166) | block | Motor control pack (if present) | Configuration for TwoMassEstimator |
| [`TwoRateSimulationResult`](inc/damp/simulation/multirate.hpp#L94) | block | Simulation / SIL harness (host) | Result of a two-rate cascade simulation |
| [`tyreus_luyben`](inc/damp/design/pid_design.hpp#L175) | function | Design-time synthesis (not PWM-rate) | Tyreus-Luyben tuning from ultimate gain and ultimate period |
| [`UKFMeasFn`](inc/damp/estimation/ukf.hpp#L65) | concept | Observers & estimators | Concept for UKF measurement functions |
| [`UKFStateFn`](inc/damp/estimation/ukf.hpp#L54) | concept | Observers & estimators | Concept for UKF state (process) functions |
| [`unbounded_bound`](inc/damp/design/qp.hpp#L67) | function | Design-time synthesis (not PWM-rate) | Sentinel bound treated as "no constraint" on that row |
| [`unload_accel_jerk`](inc/damp/trajectory/online_otg.hpp#L207) | function | Trajectory value types | Pure unload of residual a toward 0 under jmax (no reverse through 0) |
| [`UnscentedKalmanFilter`](inc/damp/estimation/ukf.hpp#L126) | block | Observers & estimators | Unscented (sigma-point) Kalman Filter for nonlinear discrete-time systems |
| [`UnscentedParams`](inc/damp/estimation/ukf.hpp#L81) | block | Observers & estimators | Tuning parameters for the scaled unscented transform |
| [`unwrap_phase_deg`](inc/damp/analysis/frequency.hpp#L125) | function | Frequency-domain analysis (host) | Unwrap phase data in degrees to avoid +/-180 discontinuities |
| [`UpperTriangle`](inc/damp/matrix/views.hpp#L114) | block | Linear algebra | Upper triangular view of a square matrix |
| [`ValidationResult`](inc/damp/analysis/identification.hpp#L150) | block | Frequency-domain analysis (host) | Validation status for a model against held-out or replayed data |
| [`VfController`](inc/damp/motor/scalar_control.hpp#L113) | block | Motor control pack (if present) | Constant Volts-per-Hertz open-loop scalar controller |
| [`voltage_circle_radius`](inc/damp/motor/foc.hpp#L268) | function | Motor control pack (if present) | Radius of the SVPWM voltage circle (max synthesizable \|V_dq\|) |
| [`WarmStartActiveSetSolver`](inc/damp/design/qp.hpp#L654) | block | Design-time synthesis (not PWM-rate) | Warm-started active-set solver policy — the damp::MPC default |
| [`with_attitude`](inc/damp/estimation/serial_arm_pose.hpp#L313) | function | Observers & estimators | Enable one-frame attitude ESKF on an existing joint design |
| [`with_schroeder_phases`](inc/damp/estimation/excitation/multi_sine.hpp#L160) | function | Observers & estimators | Assign Schroeder low-crest-factor phases to a multi-sine tone table |
| [`wrap`](inc/damp/math/math.hpp#L509) | function | Scalar math & complex | Wrap x into the half-open interval [min, max) (period max − min) |
| [`wrap_leader`](inc/damp/trajectory/cam.hpp#L72) | function | Trajectory value types | Wrap leader position into [0, L) (true modulo; supports reverse) |
| [`wrap_pi`](inc/damp/math/math.hpp#L544) | function | Scalar math & complex | Wrap an angle to [−π, π) |
| [`wrap_two_pi`](inc/damp/math/math.hpp#L562) | function | Scalar math & complex | Wrap an angle to [0, 2π) |
| [`wrapped_delta`](inc/damp/toolbox/encoder.hpp#L41) | function | Embedded helpers (controls-adjacent utilities) | Signed difference between two unsigned counter readings, wrap-safe |
| [`write_html`](inc/damp/simulation/plot_plotly.hpp#L145) | function | Simulation / SIL harness (host) | Write a plotlypp figure to HTML with damp's shell (prefer over writeHtml) |
| [`write_hybrid_trace_npy`](inc/damp/simulation/npy_export.hpp#L100) | function | Simulation / SIL harness (host) | Pack HybridSimulationResult into (N, 2+NX+NU+NY) float64: t, mode, x..., u..., y |
| [`write_npy_f64`](inc/damp/simulation/npy_export.hpp#L40) | function | Simulation / SIL harness (host) | Write a C-order float64 array as a .npy file |
| [`XyPlotOpts`](inc/damp/simulation/plot_plotly.hpp#L535) | block | Simulation / SIL harness (host) | Optional axis framing for plot_xy |
| [`XySeries`](inc/damp/simulation/plot_plotly.hpp#L517) | block | Simulation / SIL harness (host) | One series for a planar (x, y) scatter / path plot |
| [`ZetaMap`](inc/damp/power/zeta.hpp#L29) | block | Power electronics pack (if present) | Zeta CCM plant map — HardSwitchingConverter (layer 3) |
| [`ziegler_nichols`](inc/damp/design/pid_design.hpp#L75) | function | Design-time synthesis (not PWM-rate) | Ziegler-Nichols tuning from ultimate gain and ultimate period |
| [`ziegler_nichols_step`](inc/damp/design/pid_design.hpp#L124) | function | Design-time synthesis (not PWM-rate) | Ziegler-Nichols step response method (reaction curve) |
| [`ZPK`](inc/damp/systems/zpk.hpp#L177) | block | LTI systems (SS / TF / ZPK / discretize) | Zero-pole-gain (ZPK) representation of a SISO LTI system |
| [`zpk`](inc/damp/matlab.hpp#L120) | function | MATLAB®-style aliases (host) | MATLAB®-style zero-pole-gain model constructor |
| [`zpk2tf`](inc/damp/systems/zpk.hpp#L480) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for ZPK::to_transfer_function |
| [`ZPKResult`](inc/damp/systems/zpk.hpp#L399) | block | LTI systems (SS / TF / ZPK / discretize) | Result of converting a transfer function to zero-pole-gain form |
