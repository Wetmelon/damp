// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file luenberger.hpp
 * @brief Luenberger state observer: deterministic pole-placement state estimation.
 *
 * The deterministic counterpart to the Kalman filter. Given a plant
 *
 *     x[k+1] = A x[k] + B u[k]
 *     y[k]   = C x[k] + D u[k]
 *
 * an observer reconstructs the state from inputs and outputs:
 *
 *     x̂[k+1] = A x̂[k] + B u[k] + L (y[k] − C x̂[k] − D u[k])
 *
 * The estimation error e = x − x̂ evolves as e[k+1] = (A − L C) e[k], so the
 * gain L is chosen to place the eigenvalues of (A − L C) at desired locations.
 * This is the dual of state-feedback pole placement: eig(A − LC) = eig(Aᵀ − Cᵀ Lᵀ),
 * so Lᵀ is a state-feedback gain for the pair (Aᵀ, Cᵀ) and
 *
 *     L = place(Aᵀ, Cᵀ, p)ᵀ
 *
 * via the robust Kautsky–Nichols–Van Dooren routine (design::place). Because
 * place() is multi-input, the observer is multi-output (NY ≤ NX) for free; the
 * spectrum is assigned exactly and place() minimizes the eigenvector conditioning,
 * so the placed poles are less perturbation-sensitive than the single-output
 * Ackermann formula would give. Two requirements inherited from place(): complex
 * poles must be supplied as adjacent conjugate pairs, and the desired spectrum must
 * be non-defective — nudge exact repeats apart, or place a Jordan structure directly.
 *
 * A reduced-order (Gopinath) observer is also provided: it estimates only the
 * unmeasured states and reads the measured state directly from y — ideal when
 * the measurement is essentially noise-free (e.g. encoder counts).
 *
 * Example: observer for a double integrator (position measured, velocity estimated)
 * @code
 * #include "damp/estimation/luenberger.hpp"
 * #include "damp/systems/state_space.hpp"
 *
 * using namespace damp;
 *
 * constexpr StateSpace sys{
 *     .A = Matrix<2,2>{{1.0, 0.1}, {0.0, 1.0}},   // 10 Hz discrete double integrator
 *     .B = Matrix<2,1>{{0.005}, {0.1}},
 *     .C = Matrix<1,2>{{1.0, 0.0}},               // measure position only
 *     .Ts = 0.1
 * };
 *
 * // Place observer error poles inside the unit circle (faster than the plant).
 * constexpr auto result = design::luenberger(sys, ColVec<2>{0.3, 0.4});
 * static_assert(result.success);
 *
 * Luenberger observer(sys, result.as<float>());   // deploy in float
 * // In the loop: observer.step(y, u); auto xhat = observer.state();
 * @endcode
 *
 * @see "An Introduction to Observers" (Luenberger, 1971), IEEE TAC
 * @see kalman.hpp for the stochastic estimator (and MIMO outputs)
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/pole_placement.hpp" // design::place (dual placement)
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @struct LuenbergerResult
 * @brief Luenberger observer design result.
 *
 * Contains the observer gain L and the achieved observer error poles. Use
 * .as<float>() to convert for embedded deployment.
 *
 * @tparam NX Number of states
 * @tparam NY Number of outputs (NY ≤ NX)
 * @tparam T  Scalar type
 *
 * @see "An Introduction to Observers" (Luenberger, 1971)
 */
template<size_t NX, size_t NY, typename T = double>
struct LuenbergerResult {
    Matrix<NX, NY, T>            L{};            ///< Luenberger gain: x̂ ← x̂ + L (y − ŷ)
    ColVec<NX, damp::complex<T>> e{};            ///< Luenberger error poles (eigenvalues of A − LC)
    bool                         success{false}; ///< true if the system is observable and L was placed

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return LuenbergerResult<NX, NY, U>{
            L.template as<U>(),
            e.template as<damp::complex<U>>(),
            success
        };
    }

    /**
     * @brief Check whether the observer error dynamics are stable (discrete-time).
     * @return true if all error poles lie strictly inside the unit circle.
     */
    [[nodiscard]] constexpr bool is_stable() const {
        for (size_t i = 0; i < NX; ++i) {
            if (e[i].abs() >= T{1}) {
                return false;
            }
        }
        return true;
    }
};

/**
 * @brief Design a Luenberger observer by robust pole placement (matrix form).
 *
 * Places the eigenvalues of (A − L C) at @p desired_poles via the dual of
 * design::place: L = place(Aᵀ, Cᵀ, p)ᵀ. Multi-output (NY ≤ NX). Returns
 * success = false if the pair is unobservable or the spectrum cannot be placed
 * (defective target, or complex poles not in adjacent conjugate pairs).
 *
 * @note Compare with MATLAB®'s L = place(A', C', p)'.
 *
 * @param A             State transition matrix (NX × NX)
 * @param C             Output matrix (NY × NX)
 * @param desired_poles Desired error poles (complex; adjacent conjugate pairs)
 * @return LuenbergerResult with gain L and the (echoed) desired poles e
 *
 * @see design::place for the underlying robust placement and its pole requirements
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr LuenbergerResult<NX, NY, T> luenberger(
    const Matrix<NX, NX, T>&            A,
    const Matrix<NY, NX, T>&            C,
    const ColVec<NX, damp::complex<T>>& desired_poles
) {
    static_assert(NY <= NX, "An observer cannot place more poles than states (NY <= NX).");

    LuenbergerResult<NX, NY, T> result{};
    result.e = desired_poles;

    // Luenberger placement is the dual of state feedback: eig(A − LC) = eig(Aᵀ − Cᵀ Lᵀ),
    // so Lᵀ = place(Aᵀ, Cᵀ, p). place() is the robust KNV routine and is multi-input,
    // so this handles multi-output observers that the single-output Ackermann cannot.
    damp::array<damp::complex<T>, NX> poles{};
    for (size_t i = 0; i < NX; ++i) {
        poles[i] = desired_poles[i];
    }
    const auto Lt = place<NX, NY, T>(A.transpose(), C.transpose(), poles);
    if (!Lt) {
        return result; // unobservable, defective target, or bad conjugate pairing
    }
    result.L = Lt->transpose();
    result.success = true;
    return result;
}

/**
 * @brief Design a Luenberger observer from real desired poles.
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr LuenbergerResult<NX, NY, T> luenberger(
    const Matrix<NX, NX, T>& A,
    const Matrix<NY, NX, T>& C,
    const ColVec<NX, T>&     desired_poles
) {
    ColVec<NX, damp::complex<T>> poles{};
    for (size_t i = 0; i < NX; ++i) {
        poles(i, 0) = damp::complex<T>{desired_poles[i], T{0}};
    }
    return luenberger<NX, NY, T>(A, C, poles);
}

/**
 * @brief Design a Luenberger observer from a StateSpace system (complex poles).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr LuenbergerResult<NX, NY, T> luenberger(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const ColVec<NX, damp::complex<T>>&      desired_poles
) {
    return luenberger<NX, NY, T>(sys.A, sys.C, desired_poles);
}

/**
 * @brief Design a Luenberger observer from a StateSpace system (real poles).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr LuenbergerResult<NX, NY, T> luenberger(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const ColVec<NX, T>&                     desired_poles
) {
    return luenberger<NX, NY, T>(sys.A, sys.C, desired_poles);
}

/**
 * @struct ReducedLuenbergerResult
 * @brief Reduced-order (Gopinath) observer design result.
 *
 * Estimates only the NX−1 unmeasured states from a single measured output,
 * reading the measured state directly from y. Carries the runtime-ready
 * matrices; use .as<float>() for embedded deployment.
 *
 * @tparam NX Number of plant states
 * @tparam NU Number of inputs
 * @tparam T  Scalar type
 *
 * @see "On the Synthesis of Minimal-Order Observers" (Gopinath, 1971)
 */
template<size_t NX, size_t NU, typename T = double>
struct ReducedLuenbergerResult {
    static constexpr size_t NM = NX - 1; ///< Number of estimated (unmeasured) states

    Matrix<NM, NM, T>            F{};            ///< Internal dynamics (Abb − L Aab); error poles
    Matrix<NM, 1, T>             L{};            ///< Reduced gain (x̂_b = z + L y)
    Matrix<NM, 1, T>             Gy{};           ///< Internal-state update gain on y
    Matrix<NM, NU, T>            Gu{};           ///< Internal-state update gain on u
    Matrix<NX, NX, T>            Tinv{};         ///< State reconstruction: x = Tinv [y; x̂_b]
    ColVec<NM, damp::complex<T>> e{};            ///< Placed error poles
    bool                         success{false}; ///< true if observable and L was placed

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return ReducedLuenbergerResult<NX, NU, U>{
            F.template as<U>(),
            L.template as<U>(),
            Gy.template as<U>(),
            Gu.template as<U>(),
            Tinv.template as<U>(),
            e.template as<damp::complex<U>>(),
            success
        };
    }

    /**
     * @brief Check whether the reduced error dynamics are stable (discrete-time).
     */
    [[nodiscard]] constexpr bool is_stable() const {
        for (size_t i = 0; i < NM; ++i) {
            if (e[i].abs() >= T{1}) {
                return false;
            }
        }
        return true;
    }
};

/**
 * @brief Design a reduced-order (Gopinath) observer by pole placement (matrix form).
 *
 * Estimates the NX−1 unmeasured states from a single measured output, reading
 * the measured state directly from y (no filtering of the measurement). Places
 * the NX−1 error poles eig(A_bb − L A_ab) with @ref luenberger on the
 * reduced subsystem. Single measured output only (NY == 1).
 *
 * @param A             State matrix (NX × NX)
 * @param B             Input matrix (NX × NU)
 * @param C             Output matrix (1 × NX); one dominant entry (a state measurement)
 * @param desired_poles NX−1 desired error poles (complex; conjugate pairs allowed)
 * @return ReducedLuenbergerResult with runtime matrices
 *
 * @see "On the Synthesis of Minimal-Order Observers" (Gopinath, 1971)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
[[nodiscard]] constexpr ReducedLuenbergerResult<NX, NU, T> reduced_luenberger(
    const Matrix<NX, NX, T>&                A,
    const Matrix<NX, NU, T>&                B,
    const Matrix<NY, NX, T>&                C,
    const ColVec<NX - 1, damp::complex<T>>& desired_poles
) {
    static_assert(NX >= 2, "Reduced-order observer requires at least two states.");
    static_assert(NY == 1, "Reduced-order observer is single-output (NY == 1).");
    constexpr size_t NM = NX - 1;

    ReducedLuenbergerResult<NX, NU, T> result{};
    result.e = desired_poles;

    // Pick the measured-state pivot (largest |C| entry) and build the transform
    // Tf = [C; eⱼ for j ≠ m] so the first transformed state equals the output y.
    size_t m = 0;
    T      best = damp::abs(C(0, 0));
    for (size_t j = 1; j < NX; ++j) {
        const T v = damp::abs(C(0, j));
        if (v > best) {
            best = v;
            m = j;
        }
    }

    Matrix<NX, NX, T> Tf = Matrix<NX, NX, T>::zeros();
    for (size_t j = 0; j < NX; ++j) {
        Tf(0, j) = C(0, j);
    }
    size_t r = 1;
    for (size_t j = 0; j < NX; ++j) {
        if (j == m) {
            continue;
        }
        Tf(r, j) = T{1};
        ++r;
    }

    // Tinv is a stored deliverable (state reconstruction); form it via LU solve
    // Tf · Tinv = I rather than an unguarded inverse (fail closed if singular).
    const auto Tf_inv_opt = mat::lu_solve(Tf, Matrix<NX, NX, T>::identity());
    if (!Tf_inv_opt) {
        return result; // C not full rank / Tf singular: success stays false
    }
    const Matrix<NX, NX, T> Tf_inv = *Tf_inv_opt;

    // Transform to coordinates where the first state is the measurement.
    const Matrix<NX, NX, T> A_bar = Tf * A * Tf_inv;
    const Matrix<NX, NU, T> B_bar = Tf * B;

    const Matrix<1, 1, T>   A_aa = A_bar.template block<1, 1>(0, 0);
    const Matrix<1, NM, T>  A_ab = A_bar.template block<1, NM>(0, 1);
    const Matrix<NM, 1, T>  A_ba = A_bar.template block<NM, 1>(1, 0);
    const Matrix<NM, NM, T> A_bb = A_bar.template block<NM, NM>(1, 1);
    const Matrix<1, NU, T>  B_a = B_bar.template block<1, NU>(0, 0);
    const Matrix<NM, NU, T> B_b = B_bar.template block<NM, NU>(1, 0);

    // Place eig(A_bb − L A_ab) via the dual robust placement on the reduced system.
    const auto sub = luenberger<NM, 1, T>(A_bb, A_ab, desired_poles);
    if (!sub.success) {
        return result; // unobservable reduced subsystem
    }
    const Matrix<NM, 1, T> L = sub.L;

    const Matrix<NM, NM, T> F = A_bb - L * A_ab;
    result.F = F;
    result.L = L;
    result.Gy = Matrix<NM, 1, T>(F * L + (A_ba - L * A_aa));
    result.Gu = Matrix<NM, NU, T>(B_b - L * B_a);
    result.Tinv = Tf_inv;
    result.success = true;
    return result;
}

/**
 * @brief Design a reduced-order observer from real desired poles.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
[[nodiscard]] constexpr ReducedLuenbergerResult<NX, NU, T> reduced_luenberger(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NY, NX, T>& C,
    const ColVec<NX - 1, T>& desired_poles
) {
    ColVec<NX - 1, damp::complex<T>> poles{};
    for (size_t i = 0; i < NX - 1; ++i) {
        poles(i, 0) = damp::complex<T>{desired_poles[i], T{0}};
    }
    return reduced_luenberger<NX, NU, NY, T>(A, B, C, poles);
}

/**
 * @brief Design a reduced-order observer from a StateSpace system (complex poles).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr ReducedLuenbergerResult<NX, NU, T> reduced_luenberger(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const ColVec<NX - 1, damp::complex<T>>&  desired_poles
) {
    return reduced_luenberger<NX, NU, NY, T>(sys.A, sys.B, sys.C, desired_poles);
}

/**
 * @brief Design a reduced-order observer from a StateSpace system (real poles).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr ReducedLuenbergerResult<NX, NU, T> reduced_luenberger(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const ColVec<NX - 1, T>&                 desired_poles
) {
    return reduced_luenberger<NX, NU, NY, T>(sys.A, sys.B, sys.C, desired_poles);
}

} // namespace design

/**
 * @ingroup discrete_estimators
 * @brief Luenberger state observer (runtime).
 *
 * Reconstructs the plant state from inputs and outputs using a precomputed
 * gain L. Call step(y, u) once per sample; the estimate is the one-step-ahead
 * prediction x̂[k+1|k].
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam NY Number of outputs
 * @tparam T  Scalar type (default: float for embedded deployment)
 */
template<size_t NX, size_t NU, size_t NY, typename T = float>
class Luenberger {
public:
    constexpr Luenberger() = default;

    /// Construct from system matrices and a precomputed gain.
    constexpr Luenberger(
        const Matrix<NX, NX, T>& A,
        const Matrix<NX, NU, T>& B,
        const Matrix<NY, NX, T>& C,
        const Matrix<NY, NU, T>& D,
        const Matrix<NX, NY, T>& L
    )
        : A_(A), B_(B), C_(C), D_(D), L_(L), valid_(true) {}

    /// Construct from a StateSpace system and a precomputed gain.
    template<size_t NW, size_t NV>
    constexpr Luenberger(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Matrix<NX, NY, T>& L)
        : A_(sys.A), B_(sys.B), C_(sys.C), D_(sys.D), L_(L), valid_(true) {}

    /// Construct from a StateSpace system and an observer design result.
    /// Failed designs install L = 0 and valid_ = false (open-loop predictor only).
    template<size_t NW, size_t NV>
    constexpr Luenberger(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const design::LuenbergerResult<NX, NY, T>& result)
        : A_(sys.A),
          B_(sys.B),
          C_(sys.C),
          D_(sys.D),
          L_(result.success ? result.L : Matrix<NX, NY, T>{}),
          valid_(result.success) {}

    /**
     * @brief One observer recursion (predictor form):
     *
     *     x̂ ← A x̂ + B u + L (y − C x̂ − D u)
     *
     * The estimate is the one-step-ahead prediction x̂[k+1|k]; the estimation
     * error then evolves as e ← (A − L C) e, matching the placed poles.
     * If constructed from a failed design, L = 0 (open-loop prediction only).
     *
     * @param y Measured output at the current step
     * @param u Input applied at the current step
     */
    constexpr void step(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        const ColVec<NY, T> innovation = ColVec<NY, T>(y - (C_ * x_ + D_ * u));
        // L_ is zero when !valid_; still run open-loop A x + B u (no NaN path).
        x_ = ColVec<NX, T>(A_ * x_ + B_ * u + L_ * innovation);
    }

    /// true if constructed with a successful design (or from explicit L)
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /// step() under the estimator vocabulary (estimators say update, controllers say control).
    constexpr void update(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) { step(y, u); }

    /**
     * @brief Fused per-tick estimate: one observer recursion, returning the state.
     *
     * The StateEstimator concept shape — same as step()/update() plus the
     * state readout in one call.
     */
    constexpr const ColVec<NX, T>& estimate(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        step(y, u);
        return x_;
    }

    /// Current state estimate.
    [[nodiscard]] constexpr const ColVec<NX, T>& state() const { return x_; }

    /// Reset the estimate (default zero).
    constexpr void reset(const ColVec<NX, T>& x0 = ColVec<NX, T>{}) { x_ = x0; }

    /**
     * @brief Overwrite the state estimate (constraint enforcement).
     *
     * Use after step() to clamp the estimate to a physically meaningful range,
     * wrap an angle, or zero a non-physical state before the next step()
     * propagates from it. Unlike reset(), this is for per-tick constraint
     * projection, not re-initialization.
     */
    constexpr void set_state(const ColVec<NX, T>& x_new) { x_ = x_new; }
    constexpr void set_state(size_t i, T value) { x_[i] = value; }

private:
    Matrix<NX, NX, T> A_{};
    Matrix<NX, NU, T> B_{};
    Matrix<NY, NX, T> C_{};
    Matrix<NY, NU, T> D_{};
    Matrix<NX, NY, T> L_{};
    ColVec<NX, T>     x_{};
    bool              valid_{true};
};

/**
 * @ingroup discrete_estimators
 * @brief Reduced-order (Gopinath) state observer (runtime).
 *
 * Estimates the NX−1 unmeasured states from a single measured output and
 * reconstructs the full state, reading the measured state directly from y.
 * Use when the measurement is essentially noise-free (e.g. encoder counts):
 * it does no filtering of the measured channel. For noisy measurements prefer
 * the full-order Luenberger or the Kalman filter.
 *
 * @tparam NX Number of plant states
 * @tparam NU Number of inputs
 * @tparam T  Scalar type (default: float for embedded deployment)
 */
template<size_t NX, size_t NU, typename T = float>
class ReducedLuenberger {
    static constexpr size_t NM = NX - 1;

public:
    constexpr ReducedLuenberger() = default;

    constexpr explicit ReducedLuenberger(const design::ReducedLuenbergerResult<NX, NU, T>& result)
        : F_(result.F), L_(result.L), Gy_(result.Gy), Gu_(result.Gu), Tinv_(result.Tinv) {}

    /**
     * @brief One observer recursion.
     *
     * Forms the unmeasured-state estimate x̂_b = z + L y, reconstructs the full
     * state x̂ = Tinv [y; x̂_b], then advances the internal state z. The reduced
     * estimation error evolves as e ← (A_bb − L A_ab) e.
     *
     * @param y Measured output (single channel)
     * @param u Input applied at the current step
     */
    constexpr void step(const ColVec<1, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        const ColVec<NM, T> x_b = ColVec<NM, T>(z_ + L_ * y);

        ColVec<NX, T> x_bar{};
        x_bar(0, 0) = y(0, 0);
        x_bar.template block<NM, 1>(1, 0) = x_b;
        x_ = ColVec<NX, T>(Tinv_ * x_bar);

        z_ = ColVec<NM, T>(F_ * z_ + Gy_ * y + Gu_ * u);
    }

    /// Scalar-measurement convenience overload.
    constexpr void step(T y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        step(ColVec<1, T>{y}, u);
    }

    /// Current full-state estimate.
    [[nodiscard]] constexpr const ColVec<NX, T>& state() const { return x_; }

    /// Internal recursion state (estimate of the NX−1 unmeasured states in
    /// transformed coordinates), before the +L·y reconstruction.
    [[nodiscard]] constexpr const ColVec<NM, T>& internal_state() const { return z_; }

    /// Reset the internal and reconstructed states to zero.
    constexpr void reset() {
        z_ = ColVec<NM, T>{};
        x_ = ColVec<NX, T>{};
    }

    /**
     * @brief Overwrite the internal recursion state z (constraint enforcement).
     *
     * The reduced observer reconstructs the full estimate x̂ = Tinv·[y; z + L·y]
     * fresh on every step(), so writing the full state directly would be
     * discarded next tick — the persistent state is the internal z. To clamp an
     * unmeasured state, map your constraint into z. The measured channel is read
     * straight from y and is never filtered, so it needs no clamping here.
     */
    constexpr void set_internal_state(const ColVec<NM, T>& z_new) { z_ = z_new; }
    constexpr void set_internal_state(size_t i, T value) { z_[i] = value; }

private:
    Matrix<NM, NM, T> F_{};
    Matrix<NM, 1, T>  L_{};
    Matrix<NM, 1, T>  Gy_{};
    Matrix<NM, NU, T> Gu_{};
    Matrix<NX, NX, T> Tinv_{};
    ColVec<NM, T>     z_{};
    ColVec<NX, T>     x_{};
};

} // namespace damp
