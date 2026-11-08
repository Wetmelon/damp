// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file steady_state_target.hpp
 * @brief The steady-state target a law tracks through: design::steady_state_target → SteadyStateTargetResult → SteadyStateTarget
 *
 * The target holds the controlled outputs @f$ z = C_z x + D_z u + C_{d,z} \hat d @f$ (rows of
 * the measured outputs) at r under the estimated disturbance:
 * @f[
 *   \begin{bmatrix} A-I & B \\ C_z & D_z \end{bmatrix}
 *   \begin{bmatrix} x_t \\ u_t \end{bmatrix}
 *   =
 *   \begin{bmatrix} -B_d\hat d \\ r - C_{d,z}\hat d \end{bmatrix},
 * @f]
 * solved once at design time for unit d and r, so the runtime target is two matrix–vector
 * products (with NZ < NU the minimum-norm solution; with NZ > NU the equilibrium whose outputs
 * are closest to r in the least-squares sense). Shared by the laws that track through a target
 * (OFLQG, OFMPC); a law with a disturbance model re-anchors its estimator (a FixedGainObserver
 * with design::kalman's fixed gain) on it.
 *
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings,
 *      2003), AIChE Journal 49(2)
 */

#include "damp/backend.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/options.hpp"
#include "damp/design/status.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @brief The outputs a law without integrators or a disturbance model holds at r through its
 *        steady-state target (design::lqr, design::lqg): `TrackedOutputs<1>{{0}}`, in that order
 *        (the references r follow it).
 *
 * The target is the model's: x* = Nₓ·r, u* = Nᵤ·r. A model error or a constant disturbance
 * leaves an offset; integral action (lqi, lqgi) or a disturbance model (oflqg) removes it.
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin, Powell & Emami-Naeini), the reference
 *      input of state-space design (N_x, N_u)
 *
 * @tparam NZ Outputs held at r, listed in @ref outputs (0: every measured output)
 */
template<size_t NZ = 0>
struct TrackedOutputs {
    damp::array<size_t, NZ> outputs{}; ///< Indices of the measured outputs held at r, each once (none listed: all)
};

/**
 * @struct SteadyStateTargetResult
 * @brief The steady-state target maps: [x_t; u_t] per unit disturbance and per unit reference
 *
 * @tparam NX States of the filter's model before the disturbance (the plant)
 * @tparam ND Disturbance states (0: the target from r alone)
 * @tparam NZ Controlled outputs (the listed measured outputs, or all of them)
 */
template<size_t NX, size_t NU, size_t ND, size_t NZ, typename T = double>
struct [[nodiscard]] SteadyStateTargetResult {
    Matrix<NX + NU, ND, T> target_d{};                          ///< [x_t; u_t] per unit estimated disturbance
    Matrix<NX + NU, NZ, T> target_r{};                          ///< [x_t; u_t] per unit controlled-output reference
    bool                   success{false};                      ///< status == DesignStatus::ok
    DesignStatus           status{DesignStatus::not_converged}; ///< ok, or why the target is refused

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr SteadyStateTargetResult<NX, NU, ND, NZ, U> as() const {
        return SteadyStateTargetResult<NX, NU, ND, NZ, U>{target_d.template as<U>(), target_r.template as<U>(), success, status};
    }
};

namespace detail {

/** @brief The controlled outputs a spec lists: none for NoDisturbance (every output). */
[[nodiscard]] constexpr damp::array<size_t, 0> listed_outputs(const NoDisturbance&) {
    return {};
}

/** @brief The controlled outputs an InputDisturbance, OutputDisturbance, PlantDisturbance or DisturbanceModel lists. */
template<typename Spec>
[[nodiscard]] constexpr auto listed_outputs(const Spec& spec) {
    return spec.outputs;
}

/** @brief The disturbance channels of NoDisturbance: none. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXE, size_t NWE, size_t NVE>
[[nodiscard]] constexpr DisturbanceChannels<NX, NXE, NY, 0, T> disturbance_channels(
    const NoDisturbance&,
    const StateSpace<NX, NU, NY, T, NW, NV>&,
    const StateSpace<NXE, NU, NY, T, NWE, NVE>&,
    T
) {
    return {};
}

/** @brief The disturbance channels of TrackedOutputs: none. */
template<size_t NZ, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXE, size_t NWE, size_t NVE>
[[nodiscard]] constexpr DisturbanceChannels<NX, NXE, NY, 0, T> disturbance_channels(
    const TrackedOutputs<NZ>&,
    const StateSpace<NX, NU, NY, T, NW, NV>&,
    const StateSpace<NXE, NU, NY, T, NWE, NVE>&,
    T
) {
    return {};
}

/** @brief Tracked outputs: no disturbance states. */
template<size_t NZ, size_t NU, size_t NY, size_t NW>
struct disturbance_traits<TrackedOutputs<NZ>, NU, NY, NW> {
    static constexpr size_t states = 0;  ///< ND
    static constexpr size_t listed = NZ; ///< NZ listed
};

/** @brief The options steady_state_target takes: the model features of the design (a KalmanNoise is ignored). */
using target_options = ::damp::detail::option_set<is_disturbance, is_kalman_noise>;

/** @brief The controlled outputs of a target under @p Spec: those it lists, else every output. */
template<typename Spec, size_t NU, size_t NY, size_t NW>
inline constexpr size_t target_outputs = (disturbance_traits<Spec, NU, NY, NW>::listed == 0) ? NY : disturbance_traits<Spec, NU, NY, NW>::listed;

/**
 * @brief The target maps over the model's equilibria: [x_t; u_t] = s_d·d + N·w, with N a basis of
 *        the kernel of E = [A − I, B], s_d E's minimum-norm solution for −B_d, and w the least
 *        squares of ‖S·[x_t; u_t] + C_{d,z}·d − r‖²_Q, S = [C_z, D_z].
 *
 * Always an equilibrium of the model (when one exists for d) and the unique steady state when
 * there is one. A reference no steady state holds (more outputs than inputs, or a velocity on a
 * position plant) gets the closest one measured by Q, and an output Q does not weight does not
 * pull it. Q enters as W = Σ^{1/2}·Uᵀ (Q = U·Σ·Uᵀ, symmetric positive semidefinite) before the
 * pseudoinverse, so its rank cut sees the conditioning of W·S·N rather than of its square; N
 * keeps the columns outside the kernel at exactly zero, so no rank decision falls on roundoff.
 */
template<size_t NXE, size_t NU, size_t NZ, size_t ND, typename T>
[[nodiscard]] constexpr SteadyStateTargetResult<NXE, NU, ND, NZ, T> equilibrium_target(
    const Matrix<NXE, NXE + NU, T>& E,
    const Matrix<NZ, NXE + NU, T>&  S,
    const Matrix<NXE, ND, T>&       Bd,
    const Matrix<NZ, ND, T>&        Cdz,
    const Matrix<NZ, NZ, T>&        Q
) {
    constexpr size_t  NS = NXE + NU;
    const auto        e = mat::svd(E);
    const size_t      rank = mat::rank_from_svd(e);
    Matrix<NS, NS, T> N{};
    for (size_t c = rank; c < NS; ++c) {
        for (size_t i = 0; i < NS; ++i) {
            N(i, c) = e.singular_V(i, c);
        }
    }
    const auto        q = mat::svd(Q);
    Matrix<NZ, NZ, T> W{};
    for (size_t i = 0; i < NZ; ++i) {
        const T root = damp::sqrt(q.singular_values[i]);
        for (size_t j = 0; j < NZ; ++j) {
            W(i, j) = root * q.singular_U(j, i);
        }
    }
    const Matrix<NS, NZ, T>                     G = N * Matrix<NS, NZ, T>(mat::pseudo_inverse(Matrix<NZ, NS, T>(W * (S * N))) * W);
    SteadyStateTargetResult<NXE, NU, ND, NZ, T> out{};
    out.target_r = G;
    if constexpr (ND > 0) {
        //! s_d = E⁺·(−B_d) = −V·Σ⁺·Uᵀ·B_d over the rank nonzero singular values
        const Matrix<NXE, ND, T> UtB = e.singular_U.transpose() * Bd;
        Matrix<NS, ND, T>        sd{};
        for (size_t k = 0; k < rank; ++k) {
            const T inv = T{1} / e.singular_values[k];
            for (size_t i = 0; i < NS; ++i) {
                for (size_t c = 0; c < ND; ++c) {
                    sd(i, c) -= e.singular_V(i, k) * inv * UtB(k, c);
                }
            }
        }
        out.target_d = sd - (G * Matrix<NZ, ND, T>((S * sd) + Cdz));
    } else {
        static_cast<void>(Bd);
        static_cast<void>(Cdz);
    }
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

/** @brief The output weight of a target's least squares when the caller has none: the identity. */
template<typename Spec, size_t NU, size_t NY, size_t NW, typename T>
using target_weight = Matrix<target_outputs<Spec, NU, NY, NW>, target_outputs<Spec, NU, NY, NW>, T>;

/**
 * @brief The target maps on the filter's model of @p sys under the options (see
 *        steady_state_target), in the plant's time domain (@p Discrete) or sampled at @p Ts.
 *        With NZ > NU, @p Q weighs the controlled outputs' least squares (MPC's cost).
 */
template<bool Discrete, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Spec>
[[nodiscard]] constexpr auto target_design(
    const StateSpace<NX, NU, NY, T, NW, NV>&  sys,
    T                                         Ts,
    const Spec&                               spec,
    const target_weight<Spec, NU, NY, NW, T>& Q = target_weight<Spec, NU, NY, NW, T>::identity()
) {
    constexpr size_t NXE = NX;
    constexpr size_t NXL = NXE;
    constexpr size_t ND = disturbance_traits<Spec, NU, NY, NW>::states;
    constexpr size_t NZ = target_outputs<Spec, NU, NY, NW>;
    using Result = SteadyStateTargetResult<NXE, NU, ND, NZ, T>;
    const auto m = sampled_model(sys, Ts, Discrete);
    if (m.status != DesignStatus::ok) {
        return ::damp::detail::refused<Result>(m.status);
    }
    const auto rows = selected_outputs<NY>(listed_outputs(spec));
    if (!rows) {
        return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // a controlled output out of range or listed twice
    }
    const auto ch = disturbance_channels(spec, sys, m.law, Ts);
    if (!mat::all_of(ch.Bd, [](T v) { return damp::isfinite(v); }) || !mat::all_of(ch.Cd, [](T v) { return damp::isfinite(v); })) {
        return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // a non-finite B_d or C_d
    }

    //! [A − I, B; C_z, D_z]: z the listed rows of the law's model
    Matrix<NXE + NZ, NXE + NU, T> M{};
    M.template block<NXE, NXE>(0, 0) = m.law.A - Matrix<NXE, NXE, T>::identity();
    M.template block<NXE, NU>(0, NXE) = m.law.B;
    for (size_t i = 0; i < NZ; ++i) {
        for (size_t j = 0; j < NXL; ++j) {
            M(NXE + i, j) = m.law.C((*rows)[i], j);
        }
        for (size_t j = 0; j < NU; ++j) {
            M(NXE + i, NXE + j) = m.law.D((*rows)[i], j);
        }
    }
    //! d enters the state update through B_d and the controlled outputs through their C_d rows:
    //! z = C_z x + D_z u + C_{d,z} d is held at r.
    Matrix<NXE + NZ, ND, T> rhs_d{};
    rhs_d.template block<NXE, ND>(0, 0) = -ch.Bd;
    for (size_t i = 0; i < NZ; ++i) {
        for (size_t j = 0; j < ND; ++j) {
            rhs_d(NXE + i, j) = -ch.Cd((*rows)[i], j);
        }
    }
    Matrix<NXE + NZ, NZ, T> rhs_r{};
    rhs_r.template block<NZ, NZ>(NXE, 0) = Matrix<NZ, NZ, T>::identity();
    Result out{};
    if constexpr (NZ == NU) {
        const auto td = mat::solve(M, rhs_d);
        const auto tr = mat::solve(M, rhs_r);
        if (!td || !tr) {
            return ::damp::detail::refused<Result>(DesignStatus::singular); // no unique steady state for these outputs
        }
        out.target_d = td.value();
        out.target_r = tr.value();
    } else if constexpr (NZ > NU) {
        //! NZ > NU: no steady state holds every r. Least squares on the stacked system would trade
        //! the equilibrium rows against the outputs (the law would chase a state it cannot hold),
        //! so the target is the equilibrium closest to r.
        out = equilibrium_target<NXE, NU>(
            Matrix<NXE, NXE + NU, T>(M.template block<NXE, NXE + NU>(0, 0)),
            Matrix<NZ, NXE + NU, T>(M.template block<NZ, NXE + NU>(NXE, 0)),
            ch.Bd,
            Matrix<NZ, ND, T>(rhs_d.template block<NZ, ND>(NXE, 0) * T{-1}),
            Q
        );
    } else {
        const Matrix<NXE + NU, NXE + NZ, T> Mp = mat::pseudo_inverse(M);
        out.target_d = Mp * rhs_d;
        out.target_r = Mp * rhs_r;
        {
            //! NZ < NU: a steady state holds r only if [A − I, B; C_z, D_z] has full row rank. The
            //! minimum-norm solution must then solve the system to roundoff (relative to ‖M‖·‖T‖,
            //! which carries the conditioning); a rank-deficient M leaves an O(1) residual.
            const T tol = static_cast<T>(NXE + NU) * T{64} * eps<T>();
            const T res_d = Matrix<NXE + NZ, ND, T>((M * out.target_d) - rhs_d).norm();
            const T res_r = Matrix<NXE + NZ, NZ, T>((M * out.target_r) - rhs_r).norm();
            if (!(res_d <= tol * ((M.norm() * out.target_d.norm()) + rhs_d.norm())) || !(res_r <= tol * ((M.norm() * out.target_r.norm()) + rhs_r.norm()))) {
                return ::damp::detail::refused<Result>(DesignStatus::singular); // no steady state holds every r
            }
        }
    }
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

} // namespace detail

/**
 * @brief The steady-state target of a law on a discrete plant: the maps from the estimated
 *        disturbance and the reference to [x_t; u_t]
 *
 * Solves, for unit d and r,
 * [A − I, B; C_z, D_z]·[x_t; u_t] = [−B_d d; r − C_{d,z} d] on the filter's model of
 * design::kalman with the same options: C_z and D_z are the rows of the measured outputs the
 * disturbance spec lists (all of them when it lists none, or with NoDisturbance), B_d and C_d
 * its channels, C_{d,z} the same rows of C_d (none with NoDisturbance: the target from r
 * alone). With NZ = NU the system is square and solved; with NZ < NU, the minimum-norm
 * solution; with NZ > NU no steady state holds every r, and the target is the steady state
 * (A − I)·x_t + B·u_t = −B_d·d whose outputs are closest to r in the least-squares sense.
 * Offset-free control needs the steady state to exist, NZ ≤ NU, and the augment detectable
 * (design::kalman checks it).
 *
 * Refused with DesignStatus::singular when NZ = NU and the steady state is not unique, or
 * NZ < NU and no steady state holds every r (the system is rank deficient); with
 * DesignStatus::bad_weights for a controlled output out of range or listed twice, or a
 * non-finite entry of the plant, B_d or C_d.
 *
 * @see "Offset-free reference tracking with model predictive control" (Maeder, Borrelli &
 *      Morari, 2009), Automatica 45(9)
 * @see "Disturbance modeling for offset-free linear model predictive control" (Muske &
 *      Badgwell, 2002), Journal of Process Control 12(5)
 *
 * @param sys   Discrete-time plant (Ts > 0)
 * @param opts  The model features of the design, in any order: a disturbance spec
 *              (a KalmanNoise is accepted and ignored)
 * @return SteadyStateTargetResult with the maps and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto steady_state_target(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Opts&... opts) {
    detail::target_options::check<Opts...>();
    return detail::target_design<true>(
        sys, T{0}, ::damp::detail::option<detail::is_disturbance>(NoDisturbance{}, opts...)
    );
}

/**
 * @brief The steady-state target of a law on a continuous plant sampled at @p Ts (see
 *        steady_state_target(sys, opts…)): on the discrete filter's model of
 *        design::kalman(sys, Ts, noise, opts…).
 *
 * @param sys   Continuous plant (sys.Ts = 0)
 * @param Ts    Sampling time of the law [s]
 * @param opts  The model features of the design, in any order
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto steady_state_target(const StateSpace<NX, NU, NY, T, NW, NV>& sys, damp::type_identity_t<T> Ts, const Opts&... opts) {
    detail::target_options::check<Opts...>();
    return detail::target_design<false>(
        sys, Ts, ::damp::detail::option<detail::is_disturbance>(NoDisturbance{}, opts...)
    );
}

} // namespace design

/** @brief A state and the input at it: a steady state (SteadyStateTarget), or a point of a reference trajectory (OutputInverse). */
template<size_t NX, size_t NU, typename T>
struct TargetPoint {
    ColVec<NX, T> x{}; ///< x_t
    ColVec<NU, T> u{}; ///< u_t
};

/**
 * @brief The steady-state target at run time: (x_t, u_t) from the reference and the estimated
 *        disturbance, two matrix–vector products, and the re-anchoring of a FixedGainObserver on it.
 *
 * Shared by the laws that track through a target (OFLQG, MPC, OFMPC): equilibrium() for a law
 * on the absolute state, anchor() for one on the estimator's deviation (ND > 0).
 *
 * @tparam NX States of the filter's model before the disturbance; ND Disturbance states;
 *            NZ Controlled outputs
 */
template<size_t NX, size_t NU, size_t ND, size_t NZ, typename T = float, typename TargetD = Matrix<NX + NU, ND, T>, typename TargetR = Matrix<NX + NU, NZ, T>>
struct SteadyStateTarget {
    [[no_unique_address]] TargetD target_d{}; ///< [x_t; u_t] per unit d̂
    [[no_unique_address]] TargetR target_r{}; ///< [x_t; u_t] per unit r

    constexpr SteadyStateTarget() = default;

    /** @brief From design::steady_state_target (any scalar). A failed design does not compile when constant-evaluated. */
    template<typename U>
        requires(damp::is_same_v<TargetD, Matrix<NX + NU, ND, T>> && damp::is_same_v<TargetR, Matrix<NX + NU, NZ, T>>)
    constexpr SteadyStateTarget(const design::SteadyStateTargetResult<NX, NU, ND, NZ, U>& result) // NOLINT
        : target_d(result.target_d.template as<T>()), target_r(result.target_r.template as<T>()) {
        if (!result.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /** @brief The steady state (x_t, u_t) that holds z = @p r under the disturbance @p d. */
    [[nodiscard]] constexpr TargetPoint<NX, NU, T> equilibrium(const ColVec<NZ, T>& r, const ColVec<ND, T>& d = ColVec<ND, T>{}) const {
        ColVec<NX + NU, T> s{};
        if constexpr (ND > 0) {
            s = (target_d * d) + (target_r * r);
        } else {
            static_cast<void>(d);
            s = target_r * r;
        }
        TargetPoint<NX, NU, T> p{};
        for (size_t i = 0; i < NX; ++i) {
            p.x(i) = s(i);
        }
        for (size_t i = 0; i < NU; ++i) {
            p.u(i) = s(NX + i);
        }
        return p;
    }

    /**
     * @brief Move @p estimator's anchor to the target for @p r about the reference point
     *        (@p x_ref, @p u_ref) under its d̂: x_a = x_t + x*, u_a = u_t + u*.
     */
    template<typename Estimator>
        requires(ND > 0)
    constexpr void anchor(
        Estimator&           estimator,
        const ColVec<NZ, T>& r,
        const ColVec<NX, T>& x_ref = ColVec<NX, T>{},
        const ColVec<NU, T>& u_ref = ColVec<NU, T>{}
    ) const {
        const auto p = equilibrium(r, estimator.disturbance_estimate());
        estimator.anchor(ColVec<NX, T>(p.x + x_ref), ColVec<NU, T>(p.u + u_ref), x_ref, u_ref);
    }
};

/** @brief Deduce the target from its design, keeping the result's scalar: SteadyStateTarget target{res.as<float>()}. */
template<size_t NX, size_t NU, size_t ND, size_t NZ, typename T>
SteadyStateTarget(const design::SteadyStateTargetResult<NX, NU, ND, NZ, T>&) -> SteadyStateTarget<NX, NU, ND, NZ, T>;

} // namespace damp
