// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file luenberger.hpp
 * @brief The fixed-gain state observer (Luenberger) and its pole-placement designs.
 *
 * Luenberger is the fixed-gain observer: it holds a gain and the plant model, and its gain comes
 * from design::luenberger (pole placement) or design::kalman (the steady-state Kalman gain: the
 * filter DARE). It is the estimator inside LQG, LQGI, OFLQG and OFMPC, with the disturbance anchor
 * the offset-free laws use. For a covariance that updates (R or C changing per reading) use
 * KalmanFilter (kalman.hpp).
 *
 * Both design results carry MATLAB®'s two gains: the predictor gain L and the innovation gain M
 * (MATLAB®'s Mx), with L = A M. The observer runs M in the current form:
 *
 *     x̂[k|k−1] = A x̂[k−1|k−1] + B u[k−1]                   (predict)
 *     x̂[k|k]   = x̂[k|k−1] + M (y[k] − C x̂[k|k−1] − D u)      (update)
 *
 * so the error after an update evolves as e[k|k] = (I − M C) A e[k−1|k−1].
 *
 * design::luenberger places the eigenvalues of the predictor error dynamics,
 * e[k+1] = (A − L C) e[k], by the dual of state-feedback pole placement:
 * eig(A − LC) = eig(Aᵀ − Cᵀ Lᵀ), so Lᵀ is a state-feedback gain for the pair (Aᵀ, Cᵀ) and
 *
 *     L = place(Aᵀ, Cᵀ, p)ᵀ
 *
 * via the robust Kautsky–Nichols–Van Dooren routine (design::place). On a discrete plant
 * A M = L, so (I − M C) A = A⁻¹ (A − L C) A: the current form has the placed poles. On a
 * singular or ill-conditioned A, M is placed directly on (A, C A) instead (design::luenberger).
 * Because place() is multi-input, the observer is multi-output (NY ≤ NX) for free; the
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
 * constexpr auto result = design::luenberger(sys, {0.3, 0.4});
 * static constinit Luenberger observer{result.as<float>()};  // deploy in float; a failed design does not compile
 * // In the loop: auto xhat = observer.step(y, u);  // or observer.predict(u); observer.update(y);
 *
 * // The steady-state Kalman gain instead (a plant with noise inputs G, H):
 * // static constinit Luenberger kf{design::kalman(noisy_sys, {.process = Qn, .measurement = Rn}).as<float>()};
 * @endcode
 *
 * @see "An Introduction to Observers" (Luenberger, 1971), IEEE TAC
 * @see "Digital Control of Dynamic Systems" (Franklin, Powell & Workman, 1998), §8.2 (current estimator)
 * @see kalman.hpp: design::kalman (a gain this observer runs) and KalmanFilter (the recursive filter)
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/design/pole_placement.hpp" // design::place (dual placement)
#include "damp/design/stability.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp" // design::KalmanResult, a gain the observer runs
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/functions.hpp" // mat::one_norm (κ₁(A))
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @struct LuenbergerResult
 * @brief Luenberger observer design result: the plant, the gains L and M, and the achieved error poles.
 *
 * Carries the plant the gain was placed on (the ZOH model when a sample period was
 * given), so `Luenberger observer{res}` needs nothing else. Use .as<float>() to convert
 * for embedded deployment.
 *
 * | Member | What it is | MATLAB® |
 * | --- | --- | --- |
 * | L | predictor gain: e[k+1] = (A − L C) e[k] | place(A', C', p)' |
 * | M | innovation gain, A M = L: e[k\|k] = (I − M C) A e[k−1\|k−1]; zero on a continuous plant | Mx |
 * | e | achieved error poles: eig((I − M C) A), the runtime's, on a discrete plant; eig(A − L C) on a continuous one | |
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam NY Number of outputs (NY ≤ NX)
 * @tparam T  Scalar type
 * @tparam NW Process noise inputs of the plant
 * @tparam NV Measurement noise inputs of the plant
 *
 * @see "An Introduction to Observers" (Luenberger, 1971)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct [[nodiscard]] LuenbergerResult {
    StateSpace<NX, NU, NY, T, NW, NV> sys{};                               ///< Plant the gain was placed on; sys.Ts is the observer's sample time
    Matrix<NX, NY, T>                 L{};                                 ///< Predictor gain (MATLAB®'s L): x̂[k+1] = A x̂ + B u + L (y − C x̂ − D u)
    Matrix<NX, NY, T>                 M{};                                 ///< Innovation gain, A M = L (MATLAB®'s Mx): x̂ ← x̂ + M (y − C x̂ − D u); zero on a continuous plant
    damp::array<damp::complex<T>, NX> e{};                                 ///< Achieved error poles: eig((I − M C) A) when discrete (the runtime's), eig(A − L C) when continuous
    bool                              success{false};                      ///< status == DesignStatus::ok
    DesignStatus                      status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return LuenbergerResult<NX, NU, NY, U, NW, NV>{
            sys.template as<U>(),
            L.template as<U>(),
            M.template as<U>(),
            array_as<damp::complex<U>>(e),
            success,
            status,
        };
    }

    /**
     * @brief Check whether the observer error dynamics are stable
     * @return true if the design succeeded and every error pole is strictly inside the unit circle
     *         (discrete) or in the open left half plane (continuous)
     */
    [[nodiscard]] constexpr bool is_stable() const {
        if (!success) {
            return false; // a refused design's e is zero-filled
        }
        for (size_t i = 0; i < NX; ++i) {
            if (sys.Ts > T{0} ? !(e[i].abs() < T{1}) : !(e[i].real() < T{0})) {
                return false;
            }
        }
        return true;
    }
};

namespace detail {

/**
 * @brief Whether A M = L may be solved for the current-form gain: A invertible with
 *        κ₁(A) = ‖A‖₁ ‖A⁻¹‖₁ ≤ 1/√ε (see luenberger).
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr bool solvable_for_innovation_gain(const Matrix<NX, NX, T>& A) {
    const auto A_inv = mat::solve(A, Matrix<NX, NX, T>::identity());
    if (!A_inv) {
        return false;
    }
    const T kappa = mat::one_norm(A) * mat::one_norm(*A_inv);
    return damp::isfinite(kappa) && kappa * damp::sqrt(eps<T>()) <= T{1};
}

/**
 * @brief The current-form gain of a discrete design: A M = L, or, on a numerically singular A,
 *        M placed on (A, C A) or else A⁺ L, with L = A M; e = eig((I − M C) A) either way (see luenberger).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr LuenbergerResult<NX, NU, NY, T, NW, NV> with_innovation_gain(
    LuenbergerResult<NX, NU, NY, T, NW, NV>  out,
    const damp::array<damp::complex<T>, NX>& poles
) {
    //! eig((I − M C) A) = eig(Aᵀ − (C A)ᵀ Mᵀ), so Mᵀ is a state-feedback gain for (Aᵀ, (C A)ᵀ).
    const Matrix<NX, NX, T>           At = out.sys.A.transpose();
    const Matrix<NX, NY, T>           CAt = Matrix<NY, NX, T>(out.sys.C * out.sys.A).transpose();
    damp::optional<Matrix<NX, NY, T>> M{};
    if (solvable_for_innovation_gain(out.sys.A)) {
        M = mat::solve(out.sys.A, out.L); // (I − M C) A = A⁻¹ (A − L C) A: the placed poles
    }
    if (M) {
        out.M = *M;
    } else {
        //! A numerically singular: place M directly, never dividing by A. A's (near) null space
        //! is (nearly) unobservable through C A: those modes stay at (near) z = 0.
        Matrix<NY, NX, T> Mt{};
        DesignStatus      status = place_gain(At, CAt, poles, Mt);
        if (status != DesignStatus::ok) {
            //! That pair is rank-deficient by construction, and a repeated pole at 0 (two or more
            //! delay registers) reaches a rank test that reads roundoff. L places the poles, and any
            //! M with A M = L has them too: eig((I − M C) A) = eig(A − L C).
            const Matrix<NY, NX, T> from_L = Matrix<NX, NY, T>(mat::pseudo_inverse(out.sys.A) * out.L).transpose();
            if (placement_matches(At, CAt, from_L, poles)) {
                Mt = from_L;
                status = DesignStatus::ok;
            }
        }
        if (status != DesignStatus::ok) {
            out.L = Matrix<NX, NY, T>{}; // a failed design ships no gain
            out.e = damp::array<damp::complex<T>, NX>{};
            return ::damp::detail::set_status(out, status);
        }
        out.M = Mt.transpose();
        out.L = out.sys.A * out.M;
    }
    out.e = stability::closed_loop_poles(At, CAt, Matrix<NY, NX, T>(out.M.transpose())); // the runtime's error poles
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

/** @brief Observer gains on a plant whose domain is already the observer's: L = place(Aᵀ, Cᵀ, p)ᵀ, and M with A M = L when discrete. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr LuenbergerResult<NX, NU, NY, T, NW, NV> luenberger_on(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<damp::complex<T>, NX>& poles
) {
    static_assert(NY <= NX, "An observer cannot place more poles than states (NY <= NX).");
    LuenbergerResult<NX, NU, NY, T, NW, NV> out{};
    out.sys = sys;
    if (!finite_model(sys)) {
        return ::damp::detail::set_status(out, DesignStatus::bad_weights); // placement reads only A and C; the runtime runs B and D
    }
    // eig(A − LC) = eig(Aᵀ − Cᵀ Lᵀ), so Lᵀ is a state-feedback gain for (Aᵀ, Cᵀ).
    Matrix<NY, NX, T>  Lt{};
    const DesignStatus status = place_gain(Matrix<NX, NX, T>(sys.A.transpose()), Matrix<NX, NY, T>(sys.C.transpose()), poles, Lt);
    if (status != DesignStatus::ok) {
        return ::damp::detail::set_status(out, status);
    }
    out.L = Lt.transpose();
    if (sys.Ts > T{0}) {
        return with_innovation_gain(out, poles);
    }
    out.e = stability::closed_loop_poles(Matrix<NX, NX, T>(sys.A.transpose()), Matrix<NX, NY, T>(sys.C.transpose()), Lt);
    return ::damp::detail::set_status(out, DesignStatus::ok); // continuous: no current form, M = 0
}

} // namespace detail

/**
 * @brief Design a Luenberger observer by robust pole placement, in the domain of @p sys.
 *
 * Places the eigenvalues of (A − L C) at @p poles via the dual of design::place:
 * L = place(Aᵀ, Cᵀ, p)ᵀ, the predictor gain. Multi-output (NY ≤ NX). A discrete plant
 * takes z-plane poles; a continuous one takes s-plane poles and gives a continuous gain
 * (not a runtime: the Luenberger runtime refuses it under constant evaluation). Real poles
 * convert implicitly: `design::luenberger(sys, {0.3, 0.4})`.
 *
 * On a discrete plant the result also carries the innovation gain M the Luenberger runtime
 * runs in the current form, from A M = L: then (I − M C) A = A⁻¹ (A − L C) A, so the
 * current-form error e[k|k] = (I − M C) A e[k−1|k−1] has the placed poles. The solve is taken
 * only when κ₁(A) = ‖A‖₁ ‖A⁻¹‖₁ ≤ 1/√ε (A⁻¹ from mat::solve(A, I)): its error in M is about
 * κ₁(A)·ε relative, amplified along the directions A nearly erases. Past that bound M keeps
 * fewer than half its digits, and the error puts large components into x̂[k|k], which the next
 * predict erases but u = −K x̂[k|k] reads. A singular or ill-conditioned A (a pure delay, a
 * register state, a mode much faster than the sample rate) takes the direct placement instead,
 * which never divides by A: Mᵀ = place(Aᵀ, (C A)ᵀ, p), then L = A M. The modes in A's (near)
 * null space are (nearly) unobservable through C A, and the status is the one that placement
 * reports, unless M = A⁺ L has the requested poles: where that pair's rank deficit decides the
 * placement (a pole at 0 repeated, on two or more delay registers), any M with A M = L has them.
 *
 * A mode the plant forgets within a sample (an eigenvalue λ of A near 0) has a current-form
 * error pole near 0 for any bounded gain: request a pole near 0 for it. A slow pole requested
 * for it either takes a gain of order 1/|λ| or, once the mode is numerically unobservable
 * through C A, is refused with not_placeable. An exactly singular A pins those modes at z = 0.
 *
 * @c e holds the error poles the runtime has, eig((I − M C) A) computed from M (eig(A − L C)
 * on a continuous plant). It equals eig(A − L C) whenever L = A M, and differs from @p poles
 * by the rounding of the placement (and of the solve, about κ₁(A)·ε), which grows on a nearly
 * defective spectrum or a nearly unobservable mode. A continuous plant has no current form:
 * M = 0.
 *
 * Status: not_placeable if the pair is unobservable or the spectrum cannot be placed
 * (a pole repeated more than NY times, or complex poles not in adjacent conjugate
 * pairs), singular if C (or C A, on a numerically singular A) is rank-deficient.
 *
 * @note Compare with MATLAB®'s L = place(A', C', p)' and the current estimator's Mx = A \ L.
 * @see "Digital Control of Dynamic Systems" (Franklin, Powell & Workman, 1998), §8.2 (current estimator)
 *
 * @param sys   Plant (A, C used for the gain; the whole model is kept for the runtime)
 * @param poles Desired error poles, adjacent conjugate pairs
 * @return LuenbergerResult with sys, L, M and the achieved poles e
 *
 * @see design::place for the underlying robust placement and its pole requirements
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr LuenbergerResult<NX, NU, NY, T, NW, NV> luenberger(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<damp::complex<T>, NX>& poles
) {
    return detail::luenberger_on(sys, poles);
}

/**
 * @brief Design a discrete Luenberger observer for a continuous plant sampled at @p Ts.
 *
 * ZOH-discretizes @p sys, maps each s-plane pole to z = e^{s·Ts}, and places on the
 * sampled model, which the result carries. A discrete @p sys, or a Ts that is not
 * positive and finite, is refused with DesignStatus::invalid_sample_time.
 *
 * @param sys   Continuous plant (sys.Ts == 0)
 * @param poles Desired continuous error poles (s-plane), adjacent conjugate pairs
 * @param Ts    [s] Observer sample period
 * @return LuenbergerResult with the sampled sys, the discrete L and M, and the z-plane poles e
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr LuenbergerResult<NX, NU, NY, T, NW, NV> luenberger(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<damp::complex<T>, NX>& poles,
    T                                        Ts
) {
    StateSpace<NX, NU, NY, T, NW, NV> sampled{};
    const DesignStatus                status = detail::sample_plant(sys, Ts, sampled);
    if (status != DesignStatus::ok) {
        return ::damp::detail::refused<LuenbergerResult<NX, NU, NY, T, NW, NV>>(status);
    }
    return detail::luenberger_on(sampled, detail::poles_s_to_z(poles, Ts));
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
struct [[nodiscard]] ReducedLuenbergerResult {
    static constexpr size_t NM = NX - 1; ///< Number of estimated (unmeasured) states

    Matrix<NM, NM, T>                 F{};                                 ///< Internal dynamics (Abb − L Aab); error poles
    Matrix<NM, 1, T>                  L{};                                 ///< Reduced gain (x̂_b = z + L y)
    Matrix<NM, 1, T>                  Gy{};                                ///< Internal-state update gain on y
    Matrix<NM, NU, T>                 Gu{};                                ///< Internal-state update gain on u
    Matrix<NX, NX, T>                 Tinv{};                              ///< State reconstruction: x = Tinv [y; x̂_b]
    damp::array<damp::complex<T>, NM> e{};                                 ///< Achieved error poles, eig(F)
    bool                              success{false};                      ///< status == DesignStatus::ok
    DesignStatus                      status{DesignStatus::not_converged}; ///< ok, or why the design failed
    T                                 Ts{0};                               ///< Sample time of the observer: 0 continuous, > 0 discrete

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return ReducedLuenbergerResult<NX, NU, U>{
            F.template as<U>(),
            L.template as<U>(),
            Gy.template as<U>(),
            Gu.template as<U>(),
            Tinv.template as<U>(),
            array_as<damp::complex<U>>(e),
            success,
            status,
            static_cast<U>(Ts),
        };
    }

    /**
     * @brief Check whether the reduced error dynamics are stable
     * @return true if the design succeeded and every error pole is strictly inside the unit circle
     *         (discrete) or in the open left half plane (continuous)
     */
    [[nodiscard]] constexpr bool is_stable() const {
        if (!success) {
            return false; // a refused design's e is zero-filled
        }
        for (size_t i = 0; i < NM; ++i) {
            if (Ts > T{0} ? !(e[i].abs() < T{1}) : !(e[i].real() < T{0})) {
                return false;
            }
        }
        return true;
    }
};

namespace detail {

/** @brief Reduced-order observer on a plant whose domain is already the observer's. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr ReducedLuenbergerResult<NX, NU, T> reduced_luenberger_on(
    const StateSpace<NX, NU, NY, T, NW, NV>&     sys,
    const damp::array<damp::complex<T>, NX - 1>& poles
) {
    static_assert(NX >= 2, "Reduced-order observer requires at least two states.");
    static_assert(NY == 1, "Reduced-order observer is single-output (NY == 1).");
    constexpr size_t NM = NX - 1;
    using Result = ReducedLuenbergerResult<NX, NU, T>;
    const auto& A = sys.A;
    const auto& B = sys.B;
    const auto& C = sys.C;

    Result result{};
    result.Ts = sys.Ts;
    if (!finite_model(sys)) {
        return ::damp::detail::set_status(result, DesignStatus::bad_weights); // placement does not read B; Gu does
    }
    if (!mat::all_of(sys.D, [](T v) { return v == T{0}; })) {
        return ::damp::detail::set_status(result, DesignStatus::bad_weights); // the measured state is read as y: no feedthrough
    }

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
        return ::damp::detail::set_status(result, DesignStatus::singular); // C not full rank / Tf singular
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
    const Matrix<NM, NM, T> A_bb_t = A_bb.transpose();
    const Matrix<NM, 1, T>  A_ab_t = A_ab.transpose();
    Matrix<1, NM, T>        Lt{};
    const DesignStatus      status = place_gain(A_bb_t, A_ab_t, poles, Lt);
    if (status != DesignStatus::ok) {
        return ::damp::detail::set_status(result, status); // unobservable reduced subsystem
    }
    const Matrix<NM, 1, T> L = Lt.transpose();

    const Matrix<NM, NM, T> F = A_bb - L * A_ab;
    result.F = F;
    result.L = L;
    result.Gy = Matrix<NM, 1, T>(F * L + (A_ba - L * A_aa));
    result.Gu = Matrix<NM, NU, T>(B_b - L * B_a);
    result.Tinv = Tf_inv;
    result.e = stability::closed_loop_poles(A_bb_t, A_ab_t, Lt); // eig(F)
    return ::damp::detail::set_status(result, DesignStatus::ok);
}

} // namespace detail

/**
 * @brief Design a reduced-order (Gopinath) observer by pole placement, in the domain of @p sys.
 *
 * Estimates the NX−1 unmeasured states from a single measured output, reading
 * the measured state directly from y (no filtering of the measurement). Places
 * the NX−1 error poles eig(A_bb − L A_ab) with the dual placement on the reduced
 * subsystem. Single measured output only (NY == 1).
 *
 * Status: bad_weights on a non-finite model entry or a feedthrough D ≠ 0 (the measured
 * state is read as y), singular if C has no usable measured-state pivot, not_placeable if
 * the reduced subsystem is unobservable or the spectrum cannot be placed.
 *
 * @param sys   Plant (D = 0); C (1 × NX) has one dominant entry (a state measurement)
 * @param poles NX−1 desired error poles, adjacent conjugate pairs
 * @return ReducedLuenbergerResult with runtime matrices, achieved poles and Ts = sys.Ts
 *
 * @see "On the Synthesis of Minimal-Order Observers" (Gopinath, 1971)
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr ReducedLuenbergerResult<NX, NU, T> reduced_luenberger(
    const StateSpace<NX, NU, NY, T, NW, NV>&     sys,
    const damp::array<damp::complex<T>, NX - 1>& poles
) {
    return detail::reduced_luenberger_on(sys, poles);
}

/**
 * @brief Design a discrete reduced-order observer for a continuous plant sampled at @p Ts.
 *
 * Same domain rule as luenberger(sys, poles, Ts): ZOH model, z = e^{s·Ts}; a discrete
 * @p sys or a Ts that is not positive and finite is refused with
 * DesignStatus::invalid_sample_time.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr ReducedLuenbergerResult<NX, NU, T> reduced_luenberger(
    const StateSpace<NX, NU, NY, T, NW, NV>&     sys,
    const damp::array<damp::complex<T>, NX - 1>& poles,
    T                                            Ts
) {
    StateSpace<NX, NU, NY, T, NW, NV> sampled{};
    const DesignStatus                status = detail::sample_plant(sys, Ts, sampled);
    if (status != DesignStatus::ok) {
        return ::damp::detail::refused<ReducedLuenbergerResult<NX, NU, T>>(status);
    }
    return detail::reduced_luenberger_on(sampled, detail::poles_s_to_z(poles, Ts));
}

} // namespace design

/**
 * @ingroup discrete_estimators
 * @brief The fixed-gain state observer: the innovation gain M held, the current form, and the
 *        optional disturbance anchor
 *
 * Built from design::kalman (the steady-state Kalman gain: the filter DARE), from
 * design::luenberger (pole placement), or from the model and a gain. Either result carries the
 * innovation gain M (MATLAB®'s Mx) and the predictor gain L = A M (MATLAB®'s L); the observer runs
 * M. It is the estimator inside LQG, LQGI, OFLQG and OFMPC. A fixed Kalman gain costs about
 * 2–13 % of the recursive filter's multiply-accumulates per tick (2.4 % at 10 states, 6 outputs),
 * with no per-tick solve and no failure path; for a covariance that updates (R or C changing per
 * reading) use KalmanFilter.
 *
 * The tick is the current form:
 *
 * - predict(u): @f$ \hat x \leftarrow A \hat x + B u @f$, u the input the plant received over the
 *   last period;
 * - update(y[, u]): @f$ \hat x \leftarrow \hat x + M (y - C \hat x - D u) @f$, the feedthrough on
 *   the held input unless @p u is given;
 * - step(y, u): predict(u), then update(y, u) (StateEstimator).
 *
 * The error after an update evolves as @f$ e[k|k] = (I - M C) A\, e[k-1|k-1] @f$: for a
 * design::luenberger result, the poles it placed for A − L C (A M = L makes the two similar).
 *
 * The model is [x; d]: x the law's model (NXE states) and d the ND disturbance states.
 *
 * - ND = 0: the observer on the absolute state, and the input the plant holds (for the feedthrough).
 * - ND > 0: the observer holds @f$ [\hat x - x_a;\; \hat d - d_a] @f$, deviations from an anchor
 *   @f$ (x_a, u_a, d_a) @f$: measurements enter as @f$ y - C x_a - C_d d_a - D u_a @f$ and inputs
 *   as @f$ u - u_a @f$, so the state it stores, and whatever a law forms from it, stay small
 *   however far the plant sits from zero. @ref anchor moves the anchor (SteadyStateTarget::anchor
 *   moves it to the target): it shifts the deviation by the anchor's change, exactly. Between
 *   anchors the deviation carries the reference point's own motion over the period,
 *   @f$ A x^* + B u^* - x^* @f$, so a reference need not be a trajectory of the model; at a
 *   setpoint that term is zero and the deviation observer is the absolute one in exact arithmetic.
 * - An input delay is in the model discretize(plant, Ts, InputDelay<N>{θ}) returns: its registers
 *   v(k−1) … v(k−N) are states like any other, the first loaded through B.
 *
 * A failed or continuous (Ts = 0) design, or a continuous model, does not compile when
 * constant-evaluated; at run time it installs M = 0 (an open-loop predictor) and valid() is false.
 *
 * @tparam NX Observer states: the law's model, then ND disturbance states
 * @tparam NU Inputs
 * @tparam NY Outputs
 * @tparam T  Scalar type (default: float for embedded deployment)
 * @tparam ND Disturbance states among the NX (0: none)
 *
 * @see "Digital Control of Dynamic Systems" (Franklin, Powell & Workman, 1998), §8.2 (current estimator)
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings, 2003), AIChE Journal 49(2)
 * @see "Optimal State Estimation" (Simon, 2006), §5.3 (steady-state Kalman filter)
 */
template<size_t NX, size_t NU, size_t NY, typename T = float, size_t ND = 0>
class Luenberger {
public:
    static_assert(ND <= NX, "The disturbance states are among the observer's.");
    static constexpr size_t NXE = NX - ND; ///< The law's model: the states before the disturbance

    constexpr Luenberger() = default;

    /**
     * @brief From design::luenberger (any scalar): the innovation gain M, held.
     *
     * A failed or continuous design (sys.Ts = 0, a gain for ẋ) does not compile when
     * constant-evaluated; at run time it installs M = 0 and valid() = false.
     */
    template<typename U, size_t NW, size_t NV>
        requires(ND == 0)
    constexpr Luenberger(const design::LuenbergerResult<NX, NU, NY, U, NW, NV>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : sys_(plant_of(result.sys)), valid_(usable(result.success, result.status, result.sys.Ts)) {
        if (valid_) {
            M_ = result.M.template as<T>();
        }
    }

    /**
     * @brief From design::kalman (any scalar): the steady-state innovation gain M, held.
     *
     * ND = 0 accepts a design that augments a disturbance and runs it as an absolute observer on
     * [x; d]; ND > 0 is that design's anchor. A failed design does not compile when
     * constant-evaluated; at run time it installs M = 0 and valid() = false.
     */
    template<typename U, size_t NW, size_t NV, size_t NDR>
    constexpr Luenberger(const design::KalmanResult<NX, NU, NY, U, NW, NV, NDR>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : sys_(plant_of(result.sys)), valid_(usable(result.success, result.status, result.sys.Ts)) {
        static_assert(ND == 0 || NDR == ND, "The anchor's disturbance states are the design's.");
        if (valid_) {
            M_ = result.M.template as<T>();
        }
    }

    /**
     * @brief From the model and an innovation gain from any source, with the initial estimate x̂₀.
     *
     * @p M is MATLAB®'s Mx, the innovation gain: x̂ ← x̂ + M (y − C x̂ − D u). The predictor gain is
     * L = A M. A continuous model (Ts = 0) does not compile when constant-evaluated; at run time it
     * installs M = 0 and valid() = false.
     */
    template<size_t NW, size_t NV>
    constexpr Luenberger(const StateSpace<NX, NU, NY, T, NW, NV>& model, const Matrix<NX, NY, T>& M, const ColVec<NX, T>& x0 = ColVec<NX, T>{})
        : sys_(plant_of(model)), x_(x0), valid_(usable(true, design::DesignStatus::ok, model.Ts)) {
        if (valid_) {
            M_ = M;
        }
    }

    /** @brief The same observer in another scalar. */
    template<typename U>
    constexpr Luenberger(const Luenberger<NX, NU, NY, U, ND>& other) // NOLINT(google-explicit-constructor): a scalar conversion
        : sys_(other.sys_.template as<T>()),
          M_(other.M_.template as<T>()),
          x_(other.x_.template as<T>()),
          innov_(other.innov_.template as<T>()),
          u_prev_(other.u_prev_.template as<T>()),
          valid_(other.valid_) {
        if constexpr (ND > 0) {
            anchor_.x = other.anchor_.x.template as<T>();
            anchor_.u = other.anchor_.u.template as<T>();
            anchor_.d = other.anchor_.d.template as<T>();
            anchor_.drift = other.anchor_.drift.template as<T>();
        }
    }

    /** @brief Time update with the input the plant received over the last period (absolute units): x̂ ← A x̂ + B u. */
    constexpr void predict(const ColVec<NU, T>& u_applied) {
        u_prev_ = u_applied;
        x_ = (sys_.A * x_) + (sys_.B * deviation(u_applied));
        if constexpr (ND > 0) {
            x_.template segment<NXE>(0) += anchor_.drift; // the reference point's motion over the period
        }
    }

    /** @brief Measurement update with @p y (absolute units); the feedthrough takes the input the plant holds now. */
    constexpr void update(const ColVec<NY, T>& y) { update(y, u_prev_); }

    /** @brief Measurement update with an explicit feedthrough input (absolute units): x̂ ← x̂ + M (y − C x̂ − D u). */
    constexpr void update(const ColVec<NY, T>& y, const ColVec<NU, T>& u) {
        innov_ = deviation(y, sys_.C, sys_.D) - ((sys_.C * x_) + (sys_.D * deviation(u)));
        x_ += M_ * innov_;
    }

    /**
     * @brief Fused tick: predict with @p u, correct with @p y, return the stored state.
     *
     * One call is one complete observer tick (StateEstimator).
     *
     * @param y Measurement taken now
     * @param u Input applied over the previous interval
     * @return The updated state (a deviation from the anchor when ND > 0)
     */
    constexpr const ColVec<NX, T>& step(const ColVec<NY, T>& y, const ColVec<NU, T>& u) {
        predict(u);
        update(y, u);
        return x_;
    }

    /**
     * @brief Move the anchor to (@p x_t, @p u_t) under the current d̂, the reference point
     *        (@p x_ref, @p u_ref) the target is taken about (zero at a setpoint)
     *
     * The deviation shifts by the anchor's change, exactly, and the disturbance deviation returns
     * to zero. The reference point's motion A x* + B u* − x* is carried into the next predict.
     */
    constexpr void anchor(
        const ColVec<NXE, T>& x_t,
        const ColVec<NU, T>&  u_t,
        const ColVec<NXE, T>& x_ref = ColVec<NXE, T>{},
        const ColVec<NU, T>&  u_ref = ColVec<NU, T>{}
    )
        requires(ND > 0)
    {
        const ColVec<ND, T> d_hat = disturbance_estimate();
        x_.template segment<NXE>(0) += ColVec<NXE, T>(anchor_.x - x_t);
        x_.template segment<ND>(NXE) = ColVec<ND, T>{};
        anchor_.x = x_t;
        anchor_.u = u_t;
        anchor_.d = d_hat;
        anchor_.drift = ColVec<NXE, T>((sys_.A.template block<NXE, NXE>(0, 0) * x_ref) + (sys_.B.template block<NXE, NU>(0, 0) * u_ref) - x_ref);
    }

    /** @brief The stored state [x; d] the command decided now acts on: a deviation from the anchor when ND > 0, the absolute state when ND = 0. */
    [[nodiscard]] constexpr ColVec<NX, T> at_command() const { return x_; }

    /** @brief The stored state: at_command(). A deviation from the anchor when ND > 0. */
    [[nodiscard]] constexpr const ColVec<NX, T>& state() const { return x_; }

    /** @brief The estimate of the law's model x (absolute: the anchor plus the deviation when ND > 0). */
    [[nodiscard]] constexpr ColVec<NXE, T> state_estimate() const {
        ColVec<NXE, T> xe = x_.template segment<NXE>(0);
        if constexpr (ND > 0) {
            xe += anchor_.x;
        }
        return xe;
    }

    /** @brief The estimated disturbance d̂ (absolute); no entries when ND = 0. */
    [[nodiscard]] constexpr ColVec<ND, T> disturbance_estimate() const {
        ColVec<ND, T> d{};
        if constexpr (ND > 0) {
            d = x_.template segment<ND>(NXE);
            d += anchor_.d;
        }
        return d;
    }

    /** @brief x_a: the target the deviation is relative to (ND > 0). */
    [[nodiscard]] constexpr const ColVec<NXE, T>& anchor_state() const
        requires(ND > 0)
    {
        return anchor_.x;
    }

    /** @brief u_a: the anchor's input (ND > 0). */
    [[nodiscard]] constexpr const ColVec<NU, T>& anchor_input() const
        requires(ND > 0)
    {
        return anchor_.u;
    }

    /** @brief d_a: the disturbance the anchor holds under (ND > 0). */
    [[nodiscard]] constexpr const ColVec<ND, T>& anchor_disturbance() const
        requires(ND > 0)
    {
        return anchor_.d;
    }

    /** @brief The observer's model (A, B, C, D and Ts). */
    [[nodiscard]] constexpr const StateSpace<NX, NU, NY, T>& model() const { return sys_; }

    /** @brief The innovation gain M (MATLAB®'s Mx); the predictor gain is A M. */
    [[nodiscard]] constexpr const Matrix<NX, NY, T>& gain() const { return M_; }

    /** @brief The last innovation y − C x̂ − D u (in deviations when ND > 0). */
    [[nodiscard]] constexpr const ColVec<NY, T>& innovation() const { return innov_; }

    /** @brief True if built from a successful discrete design, or a discrete model and its gain. */
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /**
     * @brief Overwrite the state estimate (constraint enforcement).
     *
     * Use after a tick to clamp the estimate to a physically meaningful range, wrap an angle, or
     * zero a non-physical state before the next predict propagates from it. Unlike reset(x0),
     * this is for per-tick constraint projection, not re-initialization.
     */
    constexpr void set_state(const ColVec<NX, T>& x_new) { x_ = x_new; }

    /** @brief Overwrite one entry of the estimate. */
    constexpr void set_state(size_t i, T value) { x_[i] = value; }

    /** @brief Replace the innovation gain M. */
    constexpr void set_gain(const Matrix<NX, NY, T>& M_new) { M_ = M_new; }

    /** @brief Forget the applied input; the estimate and its anchor are left in place. */
    constexpr void reset() { u_prev_ = ColVec<NU, T>{}; }

    /** @brief Replace the estimate and clear the innovation and the held input. The model and M stay. */
    constexpr void reset(const ColVec<NX, T>& x0) {
        x_ = x0;
        innov_ = ColVec<NY, T>{};
        u_prev_ = ColVec<NU, T>{};
    }

private:
    template<size_t, size_t, size_t, typename, size_t>
    friend class Luenberger;

    struct Anchor {
        ColVec<NXE, T> x{};     ///< x_a: the steady-state target the estimate is relative to
        ColVec<NU, T>  u{};     ///< u_a: its input
        ColVec<ND, T>  d{};     ///< d_a: the disturbance it holds under
        ColVec<NXE, T> drift{}; ///< A x* + B u* − x*: the reference point's motion over a period (zero at a setpoint)
    };
    struct NoAnchor {};

    StateSpace<NX, NU, NY, T>                                            sys_{};        ///< The observer's model, [x; d]
    Matrix<NX, NY, T>                                                    M_{};          ///< The innovation gain (zero when !valid_)
    ColVec<NX, T>                                                        x_{};          ///< Estimate (a deviation from the anchor when ND > 0)
    ColVec<NY, T>                                                        innov_{};      ///< y − C x̂ − D u
    ColVec<NU, T>                                                        u_prev_{};     ///< The input the plant holds, for the feedthrough
    [[no_unique_address]] std::conditional_t<(ND > 0), Anchor, NoAnchor> anchor_{};     ///< The anchor of the deviations (ND > 0)
    bool                                                                 valid_{false}; ///< Built from a successful discrete design or model

    /** @brief A model without its noise inputs, in the observer's scalar. */
    template<typename U, size_t NW, size_t NV>
    [[nodiscard]] static constexpr StateSpace<NX, NU, NY, T> plant_of(const StateSpace<NX, NU, NY, U, NW, NV>& s) {
        return {.A = s.A.template as<T>(), .B = s.B.template as<T>(), .C = s.C.template as<T>(), .D = s.D.template as<T>(), .Ts = static_cast<T>(s.Ts)};
    }

    /** @brief Whether a gain runs: a successful design on a discrete model. Otherwise it does not compile when constant-evaluated (the error names the status). */
    template<typename U>
    [[nodiscard]] static constexpr bool usable(bool success, design::DesignStatus status, U Ts) {
        const bool ok = success && Ts > U{0};
        if (!ok && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(success ? design::DesignStatus::invalid_sample_time : status); // constant-evaluated (constinit / constexpr): a failed or continuous design does not compile
        }
        return ok;
    }

    /** @brief An input in the observer's units: a deviation from u_a when ND > 0. */
    [[nodiscard]] constexpr ColVec<NU, T> deviation(const ColVec<NU, T>& u) const {
        if constexpr (ND > 0) {
            return ColVec<NU, T>(u - anchor_.u);
        } else {
            return u;
        }
    }

    /** @brief A reading in the observer's units: y − C [x_a; d_a] − D u_a when ND > 0. */
    [[nodiscard]] constexpr ColVec<NY, T> deviation(const ColVec<NY, T>& y, const Matrix<NY, NX, T>& C, const Matrix<NY, NU, T>& D) const {
        if constexpr (ND > 0) {
            ColVec<NX, T> a{};
            a.template segment<NXE>(0) = anchor_.x;
            a.template segment<ND>(NXE) = anchor_.d;
            return ColVec<NY, T>(y - (C * a) - (D * anchor_.u));
        } else {
            return y;
        }
    }
};

/**
 * @ingroup discrete_estimators
 * @brief Reduced-order (Gopinath) state observer (runtime).
 *
 * Estimates the NX−1 unmeasured states from a single measured output and
 * reconstructs the full state, reading the measured state directly from y.
 * Use when the measurement is essentially noise-free (e.g. encoder counts):
 * it does no filtering of the measured channel. For noisy measurements prefer
 * the full-order Luenberger (with a Kalman gain) or the KalmanFilter.
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
    /**
     * @brief From a reduced-order observer design result (any scalar); converts via @c .as\<T\>().
     *
     * A failed or continuous design is ill-formed under constant evaluation, as for Luenberger.
     * At run time it leaves the observer inert (zero matrices, a zero estimate) and valid() false,
     * instead of stepping a continuous F as if it were a discrete transition.
     */
    template<typename U>
    constexpr ReducedLuenberger(const design::ReducedLuenbergerResult<NX, NU, U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : valid_(result.success && result.Ts > U{0}) {
        if (!valid_ && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_sample_time : result.status); // constant-evaluated (constinit / constexpr): a failed or continuous design does not compile
        }
        if (valid_) {
            F_ = result.F.template as<T>();
            L_ = result.L.template as<T>();
            Gy_ = result.Gy.template as<T>();
            Gu_ = result.Gu.template as<T>();
            Tinv_ = result.Tinv.template as<T>();
        }
    }

    /**
     * @brief Fused tick: advance z over the last period with @p u, reconstruct the estimate from
     *        @p y, return it.
     *
     * Advances the internal state z ← F z + G_y y⁻ + G_u u (y⁻ the previous reading), forms the
     * unmeasured-state estimate x̂_b = z + L y, and reconstructs the full state
     * x̂ = Tinv [y; x̂_b]. The reduced estimation error evolves as e ← (A_bb − L A_ab) e.
     * One call is one complete observer tick (StateEstimator).
     *
     * @param y Measured output taken now (single channel)
     * @param u Input applied over the previous interval
     * @return The full-state estimate
     */
    constexpr const ColVec<NX, T>& step(const ColVec<1, T>& y, const ColVec<NU, T>& u) {
        z_ = ColVec<NM, T>(F_ * z_ + Gy_ * y_prev_ + Gu_ * u);
        y_prev_ = y;

        ColVec<NX, T> x_bar{};
        x_bar(0, 0) = y(0, 0);
        x_bar.template block<NM, 1>(1, 0) = ColVec<NM, T>(z_ + L_ * y);
        x_ = ColVec<NX, T>(Tinv_ * x_bar);
        return x_;
    }

    /** @brief Scalar-measurement convenience overload. */
    constexpr const ColVec<NX, T>& step(T y, const ColVec<NU, T>& u) { return step(ColVec<1, T>{y}, u); }

    /** @brief Current full-state estimate. */
    [[nodiscard]] constexpr const ColVec<NX, T>& state() const { return x_; }

    /** @brief True if built from a successful discrete design; false leaves the observer inert. */
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /**
     * @brief Internal recursion state (estimate of the NX−1 unmeasured states in
     * transformed coordinates), before the +L·y reconstruction.
     */
    [[nodiscard]] constexpr const ColVec<NM, T>& internal_state() const { return z_; }

    /** @brief Reset the internal and reconstructed states and the held reading to zero. */
    constexpr void reset() {
        z_ = ColVec<NM, T>{};
        x_ = ColVec<NX, T>{};
        y_prev_ = ColVec<1, T>{};
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
    ColVec<1, T>      y_prev_{}; ///< The previous reading, which drives z over the period that followed it
    bool              valid_{false};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar: Luenberger obs{res.as<float>()}; */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
Luenberger(const design::LuenbergerResult<NX, NU, NY, T, NW, NV>&) -> Luenberger<NX, NU, NY, T>;

/** @brief The steady-state Kalman gain, keeping the result's scalar and disturbance states: Luenberger obs{design::kalman(sys, noise)}; */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND>
Luenberger(const design::KalmanResult<NX, NU, NY, T, NW, NV, ND>&) -> Luenberger<NX, NU, NY, T, ND>;

/** @brief From the model and an innovation gain: Luenberger obs{sys, M}; */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
Luenberger(const StateSpace<NX, NU, NY, T, NW, NV>&, const Matrix<NX, NY, T>&) -> Luenberger<NX, NU, NY, T>;

/** @brief From the model, an innovation gain and the initial estimate: Luenberger obs{sys, M, x0}; */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
Luenberger(const StateSpace<NX, NU, NY, T, NW, NV>&, const Matrix<NX, NY, T>&, const ColVec<NX, T>&) -> Luenberger<NX, NU, NY, T>;

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<size_t NX, size_t NU, typename T>
ReducedLuenberger(const design::ReducedLuenbergerResult<NX, NU, T>&) -> ReducedLuenberger<NX, NU, T>;

} // namespace damp
