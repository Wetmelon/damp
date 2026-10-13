// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file synthesis.hpp
 * @brief Controller synthesis APIs
 *
 * @defgroup synthesis Controller Synthesis APIs
 * @brief High-level design bundles combining synthesis, analysis hooks, and runtime objects
 *
 * Provides thin orchestration helpers that connect existing primitives:
 * - design::discrete_lqg(...)
 * - closed-loop analysis models
 * - runtime LQG controller bundles
 * - optional SISO PR internal model composition
 */

#include <cstddef>

#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {

/**
 * @brief Analysis-oriented models produced from an LQG design
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY>
struct LQGAnalysisModels {
    StateSpace<NX, NU, NY, T, NW, NV> plant{};                      ///< Original plant model
    StateSpace<NX, NU, NY, T, NW, NV> state_feedback_closed_loop{}; ///< A_cl = A - B*K
    Matrix<NX, NX, T>                 observer_error_dynamics{};    ///< Predicted-error map A(I − LC)
    StateSpace<2 * NX, NU, NY, T>     augmented_closed_loop{};      ///< [x; xhat] closed-loop model
};

/**
 * @brief Runtime bundle for LQG control
 */
template<size_t NX, size_t NU, size_t NY, typename T = float, size_t NW = NX, size_t NV = NY>
struct LQGRuntimeBundle {
    design::LQGResult<NX, NU, NY, T, NW, NV> design{};     ///< Design at runtime precision
    LQG<NX, NU, NY, T, NW, NV>               controller{}; ///< Runtime LQG built from @c design

    /**
     * @brief One discrete control tick: update estimator, compute control, predict next
     */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NY, T>& y) {
        controller.update(y);
        const auto u = controller.control();
        controller.predict(u);
        return u;
    }

    /**
     * @brief One discrete control tick tracking a state reference (u = -K*(x̂ - x_ref))
     *
     * Feedback-only 2-DOF servo with no integral action, so it droops without an input
     * feedforward — see LQG::control(const ColVec<NX,T>&).
     */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NY, T>& y, const ColVec<NX, T>& x_ref) {
        controller.update(y);
        const auto u = controller.control(x_ref);
        controller.predict(u);
        return u;
    }

    /**
     * @brief Clear the estimator state estimate (covariance stays steady-state)
     */
    constexpr void reset() {
        controller.kf.set_state(ColVec<NX, T>{});
    }
};

/**
 * @brief Synthesis artifact bundle: design + analysis models + runtime bundle
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY, typename TRuntime = float>
struct LQGArtifacts {
    design::LQGResult<NX, NU, NY, T, NW, NV>       design{};       ///< Design-precision LQG result
    LQGAnalysisModels<NX, NU, NY, T, NW, NV>       models{};       ///< Closed-loop analysis models
    LQGRuntimeBundle<NX, NU, NY, TRuntime, NW, NV> runtime{};      ///< Ready-to-run controller
    bool                                           success{false}; ///< Mirrors design.success
};

/**
 * @brief Runtime bundle for SISO LQG + PR internal model compensation
 */
template<size_t NX, size_t NW = NX, size_t NV = 1, typename T = float>
struct LQGPRRuntimeBundle {
    design::LQGResult<NX, 1, 1, T, NW, NV> lqg_design{}; ///< LQG design at runtime precision
    design::PRResult<T>                    pr_design{};  ///< PR design at runtime precision
    LQG<NX, 1, 1, T, NW, NV>               lqg{};        ///< Runtime LQG
    PRController<T>                        pr{};         ///< Runtime PR internal model

    /**
     * @brief One control tick with PR internal model on tracking error
     */
    [[nodiscard]] constexpr ColVec<1, T> step(T reference, T measurement) {
        const ColVec<1, T> y{measurement};
        lqg.update(y);

        ColVec<1, T> u = lqg.control();
        u(0, 0) += pr.control(reference - measurement);

        lqg.predict(u);
        return u;
    }

    /**
     * @brief Clear the estimator state estimate and the PR internal model
     */
    constexpr void reset() {
        lqg.kf.set_state(ColVec<NX, T>{});
        pr.reset();
    }
};

/**
 * @brief Synthesis artifact bundle for SISO LQG + PR
 */
template<size_t NX, size_t NW = NX, size_t NV = 1, typename T = double, typename TRuntime = float>
struct LQGPRArtifacts {
    LQGArtifacts<NX, 1, 1, T, NW, NV, TRuntime> base{};         ///< LQG design, models, runtime
    LQGPRRuntimeBundle<NX, NW, NV, TRuntime>    runtime_pr{};   ///< LQG + PR runtime
    bool                                        success{false}; ///< LQG and PR designs both succeeded
};

/**
 * @brief Analysis-oriented models produced from an LQI design
 * @tparam NI Integrated outputs (default: every output)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0, size_t NI = NY>
struct LQIAnalysisModels {
    StateSpace<NX, NU, NY, T, NW, NV> plant{};                       ///< Original plant model
    StateSpace<NX + NI, NU, NY, T>    augmented_servo_closed_loop{}; ///< [x; xi] servo closed-loop model
};

/**
 * @brief Analysis-oriented models produced from an LQGI design
 * @tparam NI Integrated outputs (default: every output)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY, size_t NI = NY>
struct LQGIAnalysisModels {
    StateSpace<NX, NU, NY, T, NW, NV> plant{};                       ///< Original plant model
    StateSpace<NX + NI, NU, NY, T>    augmented_servo_closed_loop{}; ///< [x; xi] servo closed-loop model
    Matrix<NX, NX, T>                 observer_error_dynamics{};     ///< Predicted-error map A(I − LC)
};

/**
 * @brief Synthesis artifact bundle: LQI servo design + analysis + ready-to-run controller
 *
 * The runtime is a plain LQI at @p TRuntime precision — it owns the integral state, so
 * the tick is `controller.control(r, y, x)` (or the `x_ref` overload for trajectory tracking).
 * @tparam NI Integrated outputs (default: every output)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0, typename TRuntime = float, size_t NI = NY>
struct LQIArtifacts {
    design::LQIResult<NX, NU, NY, T, NI>         design{};       ///< Design-precision LQI result
    LQIAnalysisModels<NX, NU, NY, T, NW, NV, NI> models{};       ///< Servo closed-loop model
    LQI<NX, NU, NY, TRuntime, NI>                runtime{};      ///< Ready-to-run controller
    bool                                         success{false}; ///< Mirrors design.success
};

/**
 * @brief Synthesis artifact bundle: LQGI servo design + analysis + ready-to-run controller
 *
 * The runtime is a plain @ref LQGI at @p TRuntime precision — estimator, integral state, and
 * applied-input memory all live in it, so the tick is `controller.step(r, y)`, or the
 * saturation-aware `feedback` / `commit` pair when the actuator can clamp.
 * @tparam NI Integrated outputs (default: every output)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY, typename TRuntime = float, size_t NI = NY>
struct LQGIArtifacts {
    design::LQGIResult<NX, NU, NY, T, NW, NV, NI> design{};       ///< Design-precision LQGI result
    LQGIAnalysisModels<NX, NU, NY, T, NW, NV, NI> models{};       ///< Servo / observer analysis models
    LQGI<NX, NU, NY, TRuntime, NW, NV, NI>        runtime{};      ///< Ready-to-run controller
    bool                                          success{false}; ///< Mirrors design.success
};

namespace detail {

/// LQI servo closed loop on [x; xi]: xi integrates r − Ci·y, the external input adds to u.
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateSpace<NX + NI, NU, NY, T> lqi_servo_closed_loop(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NU, NX + NI, T>&            K,
    const Matrix<NI, NY, T>&                 Ci
) {
    Matrix<NX + NI, NX + NI, T> A_aug{};
    A_aug.template block<NX, NX>(0, 0) = sys.A;
    A_aug.template block<NI, NX>(NX, 0) = -(Ci * sys.C);
    A_aug.template block<NI, NI>(NX, NX) = Matrix<NI, NI, T>::identity();
    Matrix<NX + NI, NU, T> B_aug{};
    B_aug.template block<NX, NU>(0, 0) = sys.B;
    Matrix<NY, NX + NI, T> C_aug{};
    C_aug.template block<NY, NX>(0, 0) = sys.C;
    return StateSpace<NX + NI, NU, NY, T>{
        .A = A_aug - (B_aug * K),
        .B = B_aug,
        .C = C_aug,
        .D = sys.D,
        .G = Matrix<NX + NI, 0, T>{},
        .H = Matrix<NY, 0, T>{},
        .Ts = sys.Ts,
    };
}

} // namespace detail

/**
 * @brief Build analysis models from an LQG design
 *
 * The augmented [x; x̂] model assumes a strictly proper plant (D = 0): the
 * estimator update uses y = C x. Its D field is copied for reference only.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY>
[[nodiscard]] constexpr LQGAnalysisModels<NX, NU, NY, T, NW, NV> build_lqg_analysis_models(
    const StateSpace<NX, NU, NY, T, NW, NV>&        sys,
    const design::LQGResult<NX, NU, NY, T, NW, NV>& lqg
) {
    const auto&             K = lqg.lqr.K;
    const auto&             L = lqg.kalman.L;
    const Matrix<NX, NX, T> I = Matrix<NX, NX, T>::identity();
    const Matrix<NX, NX, T> ImLC = I - (L * sys.C);
    const Matrix<NX, NX, T> AmBK = sys.A - (sys.B * K);
    const Matrix<NX, NX, T> A_sf = AmBK;

    Matrix<2 * NX, 2 * NX, T> A_aug{};
    Matrix<2 * NX, NU, T>     B_aug{};
    Matrix<NY, 2 * NX, T>     C_aug{};

    // Current estimator, D = 0 in the joint state (u = −K x̂(k|k), then predict).
    // External input is added to the applied u on both the plant and the predictor.
    A_aug.template block<NX, NX>(0, 0) = sys.A - ((sys.B * (K * L)) * sys.C);
    A_aug.template block<NX, NX>(0, NX) = -(sys.B * (K * ImLC));
    A_aug.template block<NX, NX>(NX, 0) = (AmBK * L) * sys.C;
    A_aug.template block<NX, NX>(NX, NX) = AmBK * ImLC;

    B_aug.template block<NX, NU>(0, 0) = sys.B;
    B_aug.template block<NX, NU>(NX, 0) = sys.B;

    C_aug.template block<NY, NX>(0, 0) = sys.C;

    return LQGAnalysisModels<NX, NU, NY, T, NW, NV>{
        .plant = sys,
        .state_feedback_closed_loop = StateSpace<NX, NU, NY, T, NW, NV>{
            .A = A_sf,
            .B = sys.B,
            .C = sys.C,
            .D = sys.D,
            .G = sys.G,
            .H = sys.H,
            .Ts = sys.Ts,
        },
        .observer_error_dynamics = sys.A * ImLC,
        .augmented_closed_loop = StateSpace<2 * NX, NU, NY, T>{
            .A = A_aug,
            .B = B_aug,
            .C = C_aug,
            .D = sys.D,
            .G = Matrix<2 * NX, 0, T>{},
            .H = Matrix<NY, 0, T>{},
            .Ts = sys.Ts,
        },
    };
}

/**
 * @brief Build analysis models from an LQI servo design (its Ci picks the integrated outputs)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0, size_t NI = NY>
[[nodiscard]] constexpr LQIAnalysisModels<NX, NU, NY, T, NW, NV, NI> build_lqi_analysis_models(
    const StateSpace<NX, NU, NY, T, NW, NV>&    sys,
    const design::LQIResult<NX, NU, NY, T, NI>& lqi
) {
    return LQIAnalysisModels<NX, NU, NY, T, NW, NV, NI>{
        .plant = sys,
        .augmented_servo_closed_loop = detail::lqi_servo_closed_loop(sys, lqi.K, lqi.Ci),
    };
}

/**
 * @brief Build analysis models from an LQGI servo design (its Ci picks the integrated outputs)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY, size_t NI = NY>
[[nodiscard]] constexpr LQGIAnalysisModels<NX, NU, NY, T, NW, NV, NI> build_lqgi_analysis_models(
    const StateSpace<NX, NU, NY, T, NW, NV>&             sys,
    const design::LQGIResult<NX, NU, NY, T, NW, NV, NI>& lqgi
) {
    return LQGIAnalysisModels<NX, NU, NY, T, NW, NV, NI>{
        .plant = sys,
        .augmented_servo_closed_loop = detail::lqi_servo_closed_loop(sys, lqgi.lqi.K, lqgi.lqi.Ci),
        .observer_error_dynamics = sys.A * (Matrix<NX, NX, T>::identity() - (lqgi.kalman.L * sys.C)),
    };
}

/**
 * @brief Synthesize the full LQG artifact bundle in one call
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY, typename TRuntime = float>
[[nodiscard]] constexpr LQGArtifacts<NX, NU, NY, T, NW, NV, TRuntime> lqg_bundle(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, NX, T>&                 Q_lqr,
    const Matrix<NU, NU, T>&                 R_lqr,
    const KalmanNoise<NW, NV, T>&            noise = {},
    const Matrix<NX, NU, T>&                 N = Matrix<NX, NU, T>{}
) {
    const auto design_r = design::discrete_lqg(sys, Q_lqr, R_lqr, noise, N);
    const auto models = build_lqg_analysis_models(sys, design_r);
    const auto runtime_design = design_r.template as<TRuntime>();

    return LQGArtifacts<NX, NU, NY, T, NW, NV, TRuntime>{
        .design = design_r,
        .models = models,
        .runtime = LQGRuntimeBundle<NX, NU, NY, TRuntime, NW, NV>{
            .design = runtime_design,
            .controller = LQG<NX, NU, NY, TRuntime, NW, NV>{runtime_design},
        },
        .success = design_r.success,
    };
}

/**
 * @brief Synthesize a SISO LQG + PR design with internal-model compensation
 */
template<size_t NX, size_t NW = NX, size_t NV = 1, typename T = double, typename TRuntime = float>
[[nodiscard]] constexpr LQGPRArtifacts<NX, NW, NV, T, TRuntime> lqg_pr_bundle(
    const StateSpace<NX, 1, 1, T, NW, NV>& sys,
    const Matrix<NX, NX, T>&               Q_lqr,
    const Matrix<1, 1, T>&                 R_lqr,
    const KalmanNoise<NW, NV, T>&          noise,
    const design::PRResult<T>&             pr,
    const Matrix<NX, 1, T>&                N = Matrix<NX, 1, T>{}
) {
    const auto base = lqg_bundle<NX, 1, 1, T, NW, NV, TRuntime>(sys, Q_lqr, R_lqr, noise, N);
    const auto lqg_design_rt = base.design.template as<TRuntime>();
    const auto pr_design_rt = pr.template as<TRuntime>();

    return LQGPRArtifacts<NX, NW, NV, T, TRuntime>{
        .base = base,
        .runtime_pr = LQGPRRuntimeBundle<NX, NW, NV, TRuntime>{
            .lqg_design = lqg_design_rt,
            .pr_design = pr_design_rt,
            .lqg = LQG<NX, 1, 1, TRuntime, NW, NV>{lqg_design_rt},
            .pr = PRController<TRuntime>{pr_design_rt},
        },
        .success = base.success && pr.success,
    };
}

namespace detail {

/// Wrap an LQI design with its analysis models and the runtime at @p TRuntime.
template<typename TRuntime, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI>
[[nodiscard]] constexpr LQIArtifacts<NX, NU, NY, T, NW, NV, TRuntime, NI> lqi_artifacts(
    const StateSpace<NX, NU, NY, T, NW, NV>&    sys,
    const design::LQIResult<NX, NU, NY, T, NI>& design_r
) {
    return LQIArtifacts<NX, NU, NY, T, NW, NV, TRuntime, NI>{
        .design = design_r,
        .models = build_lqi_analysis_models(sys, design_r),
        .runtime = LQI<NX, NU, NY, TRuntime, NI>{design_r.template as<TRuntime>()},
        .success = design_r.success,
    };
}

/// Wrap an LQGI design with its analysis models and the runtime at @p TRuntime.
template<typename TRuntime, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI>
[[nodiscard]] constexpr LQGIArtifacts<NX, NU, NY, T, NW, NV, TRuntime, NI> lqgi_artifacts(
    const StateSpace<NX, NU, NY, T, NW, NV>&             sys,
    const design::LQGIResult<NX, NU, NY, T, NW, NV, NI>& design_r
) {
    // LQGI's converting constructor seeds the filter from the design (steady-state P, or
    // identity if the DARE did not converge) and takes the LQI gain as-is.
    return LQGIArtifacts<NX, NU, NY, T, NW, NV, TRuntime, NI>{
        .design = design_r,
        .models = build_lqgi_analysis_models(sys, design_r),
        .runtime = LQGI<NX, NU, NY, TRuntime, NW, NV, NI>{design_r.template as<TRuntime>()},
        .success = design_r.success,
    };
}

} // namespace detail

/**
 * @brief Synthesize the full LQI servo artifact bundle in one call
 *
 * Integrates every output. To integrate a subset, add the last argument @ref discrete_lqi
 * takes: an index list, `lqi_bundle(sys, Q_aug, R, {0, 3})`, or an NI×NY selector Ci; Q_aug
 * is then (NX+NI)×(NX+NI) and the runtime takes NI references. On a model with delay states
 * Q_aug may be sized to the plant, (NX_plant+NI)×(NX_plant+NI), as in @ref discrete_lqi.
 */
template<
    size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0, typename TRuntime = float,
    size_t NQA = NX + NY>
[[nodiscard]] constexpr LQIArtifacts<NX, NU, NY, T, NW, NV, TRuntime> lqi_bundle(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q_aug,
    const Matrix<NU, NU, T>&                 R
) {
    return detail::lqi_artifacts<TRuntime>(sys, design::discrete_lqi(sys, Q_aug, R));
}

/// LQI bundle integrating the outputs @p Ci selects (NI×NY): integrator i follows r_i − (Ci·y)_i.
template<
    size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename TRuntime = float, size_t NI,
    size_t NQA = NX + NI>
[[nodiscard]] constexpr LQIArtifacts<NX, NU, NY, T, NW, NV, TRuntime, NI> lqi_bundle(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const Matrix<NI, NY, T>&                 Ci
) {
    return detail::lqi_artifacts<TRuntime>(sys, design::discrete_lqi(sys, Q_aug, R, Ci));
}

/// LQI bundle integrating the listed outputs, e.g. `lqi_bundle(sys, Q_aug, R, {0, 3})`.
/// An index outside [0, NY) reports success = false.
template<
    size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename TRuntime = float, size_t NI,
    size_t NQA = NX + NI>
[[nodiscard]] constexpr LQIArtifacts<NX, NU, NY, T, NW, NV, TRuntime, NI> lqi_bundle(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const size_t (&outputs)[NI]
) {
    return detail::lqi_artifacts<TRuntime>(sys, design::discrete_lqi(sys, Q_aug, R, outputs));
}

/**
 * @brief Synthesize the full LQGI servo artifact bundle in one call
 *
 * Integrates every output. To integrate a subset, add the last argument @ref discrete_lqgi
 * takes: an index list, `lqgi_bundle(sys, Q_aug, R, noise, {0, 3})`, or an NI×NY
 * selector Ci; Q_aug is then (NX+NI)×(NX+NI) and the runtime takes NI references. On a model
 * with delay states Q_aug may be sized to the plant, (NX_plant+NI)×(NX_plant+NI), as in
 * @ref discrete_lqgi.
 */
template<
    size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY, typename TRuntime = float,
    size_t NQA = NX + NY>
[[nodiscard]] constexpr LQGIArtifacts<NX, NU, NY, T, NW, NV, TRuntime> lqgi_bundle(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const KalmanNoise<NW, NV, T>&            noise = {}
) {
    return detail::lqgi_artifacts<TRuntime>(sys, design::discrete_lqgi(sys, Q_aug, R, noise));
}

/// LQGI bundle integrating the outputs @p Ci selects (NI×NY).
template<
    size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename TRuntime = float, size_t NI,
    size_t NQA = NX + NI>
[[nodiscard]] constexpr LQGIArtifacts<NX, NU, NY, T, NW, NV, TRuntime, NI> lqgi_bundle(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const KalmanNoise<NW, NV, T>&            noise,
    const Matrix<NI, NY, T>&                 Ci
) {
    return detail::lqgi_artifacts<TRuntime>(sys, design::discrete_lqgi(sys, Q_aug, R, noise, Ci));
}

/// LQGI bundle integrating the listed outputs, e.g. `lqgi_bundle(sys, Q_aug, R, noise, {0, 3})`.
/// An index outside [0, NY) reports success = false.
template<
    size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename TRuntime = float, size_t NI,
    size_t NQA = NX + NI>
[[nodiscard]] constexpr LQGIArtifacts<NX, NU, NY, T, NW, NV, TRuntime, NI> lqgi_bundle(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const KalmanNoise<NW, NV, T>&            noise,
    const size_t (&outputs)[NI]
) {
    return detail::lqgi_artifacts<TRuntime>(sys, design::discrete_lqgi(sys, Q_aug, R, noise, outputs));
}

} // namespace design
} // namespace damp
