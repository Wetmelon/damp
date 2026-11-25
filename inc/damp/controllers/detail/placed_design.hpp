// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private: a public header that exports it provides its names

/**
 * @file placed_design.hpp
 * @brief The placed halves of a state-feedback design: the gain at Poles on the law's model
 *        (placed_law_design), and the filter's gain at ObserverPoles on the filter's model
 *        (placed_filter_design), each in the shape the LQ and Kalman halves have
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/law_detail.hpp"
#include "damp/controllers/detail/lq_cost.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {
namespace detail {

/**
 * @brief Poles lists the count its gain places: the plant's states plus the integrators. A
 *        compile error otherwise, the two counts its template arguments.
 */
template<size_t kPolesExpected, size_t kPolesGiven>
consteval bool controller_pole_count() {
    static_assert(
        kPolesExpected == kPolesGiven,
        "Poles lists one pole per plant state and per integrator of an IntegralAction (the comparison: that count == the poles listed). "
        "An InputDelay's delay states take theirs at 0, and a disturbance spec's states none."
    );
    return kPolesExpected == kPolesGiven;
}

/**
 * @brief ObserverPoles lists the count its filter places: the plant's states plus the disturbance
 *        states. A compile error otherwise, the two counts its template arguments.
 */
template<size_t kPolesExpected, size_t kPolesGiven>
consteval bool observer_pole_count() {
    static_assert(
        kPolesExpected == kPolesGiven,
        "ObserverPoles lists one pole per plant state and per state of the disturbance spec (the comparison: that count == the poles listed). "
        "The delay states of an InputDelay and an OutputDelay take theirs at 0."
    );
    return kPolesExpected == kPolesGiven;
}

/** @brief The listed poles mapped into the model's domain, z = e^{s·Ts} on a sampled one (@p Sampled), then 0 for the @p NZ states placed deadbeat. */
template<size_t NZ, bool Sampled, size_t NP, typename T>
[[nodiscard]] constexpr damp::array<damp::complex<T>, NP + NZ> with_deadbeat(const damp::array<damp::complex<T>, NP>& listed, T Ts) {
    const auto                             mapped = Sampled ? poles_s_to_z(listed, Ts) : listed;
    damp::array<damp::complex<T>, NP + NZ> out{};
    for (size_t i = 0; i < NP; ++i) {
        out[i] = mapped[i];
    }
    return out;
}

/**
 * @brief The gain at @p placed on the law's model: on the plant, or with the integrators of @p Ci
 *        (NI > 0) and the delay states of @p delay (an InputDelay) as law_design and
 *        delayed_law_design build it, the delay states' poles at 0.
 *
 * The listed poles are the plant's and the integrators': s-plane on the sampled form, mapped to
 * z = e^{s·Ts} on the discrete model as design::place(sys, poles, Ts) maps them; z-plane on a
 * discrete plant; s-plane for the continuous gain of a continuous plant without Ts (analysis, as
 * the LQ design's CARE). The delay states hold the commands still on their way, and the gain
 * empties them (deadbeat): N·NU poles at 0 on top of the list. Placement is design::place's on
 * the law's (A, B): KNV, the Jordan construction and Ackermann's gain for one input, certified
 * against the request. With one input the gain is unique and takes any delay; with several, a
 * pole is placed at most NU times, so N ≥ 2 delay states per input (a deadbeat chain, which is
 * defective) is refused with not_placeable.
 *
 * @return The design's @p Result as law_design fills it for the LQ weights, S zero (no Riccati equation); status
 *         invalid_parameter for a pole that is not finite, sampled_model's for the plant,
 *         design::place's otherwise (not_placeable, singular)
 */
template<typename Result, bool Discrete, size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NP, typename In>
[[nodiscard]] constexpr Result placed_law_design(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const PlacedPoles<NP, T>& placed, const Matrix<NI, NY, T>& Ci, const In& delay) {
    constexpr size_t NRU = input_delay_states<In>;
    constexpr size_t NXL = NX + (NRU * NU); // the law's model: the plant, then the delay states
    constexpr size_t NK = NXL + NI;
    static_assert(NP + (NRU * NU) == NK, "The poles listed are the plant's and the integrators'.");
    if (!finite_poles(placed.poles)) {
        return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter); // before z = e^{s·Ts} maps −∞ to 0
    }
    if constexpr (Discrete && NI == 0 && NRU == 0) {
        if (sys.is_continuous()) {
            const auto gain = place_on(sys, placed.poles); // the continuous gain, for analysis
            Result     out{};
            out.K = gain.K;
            out.e = gain.e;
            return ::damp::detail::set_status(out, gain.status);
        }
    }
    const auto m = sampled_model(sys, Ts, Discrete, delay);
    if (m.status != DesignStatus::ok) {
        return ::damp::detail::refused<Result>(m.status);
    }
    const T                  h = Discrete ? sys.Ts : Ts;
    Result                   out{.Ts = h};
    StateSpace<NK, NU, 1, T> law{.Ts = h};
    if constexpr (NI == 0) {
        law.A = m.law.A;
        law.B = m.law.B;
    } else {
        const auto im = integral_model(m.law.A, m.law.B, m.law.C, m.law.D, Ci, h); // as the LQ design's
        law.A = im.A;
        law.B = im.B;
        out.integral = IntegralDesign<NI, NXL, NU, NY, T>{Ci, Ci * m.law.C, Ci * m.law.D, h};
    }
    const auto gain = place_on(law, with_deadbeat<NRU * NU, !Discrete>(placed.poles, Ts));
    out.K = gain.K;
    out.e = gain.e;
    return ::damp::detail::set_status(out, gain.status);
}

/** @brief The noise inputs of the filter's model as kalman_design gives them: the plant's or its states, then d's or none. */
template<bool Discrete, size_t ND, size_t NX, size_t NU, size_t NY, size_t NW, typename Delay>
inline constexpr size_t placed_noise_inputs = (ND > 0) ? (Discrete ? (NW + ND) : (NX + delay_states<Delay, NU, NY> + ND))
                                                       : (Discrete ? NW : ((output_delay_states<Delay> > 0) ? (NX + delay_states<Delay, NU, NY>) : NX));

/**
 * @brief The filter's model a design on @p sys with @p spec and @p delay estimates on, as
 *        kalman_design builds it: the plant, its delay states, then the disturbance's, with the noise
 *        inputs G the Kalman design gives it (none of the noise is read here).
 */
template<bool Discrete, size_t ND, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Spec, typename Delay>
[[nodiscard]] constexpr auto placed_filter_model(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const Spec& spec, const Delay& delay) {
    constexpr size_t                      NXE = NX + delay_states<Delay, NU, NY>;
    constexpr size_t                      NXK = NXE + ND;
    constexpr size_t                      NWK = placed_noise_inputs<Discrete, ND, NX, NU, NY, NW, Delay>;
    SampledModel<NXK, NU, NY, T, NWK, NV> out{};
    const auto                            m = sampled_model(sys, Ts, Discrete, delay);
    out.status = m.status;
    if (m.status != DesignStatus::ok) {
        return out;
    }
    const T h = Discrete ? sys.Ts : Ts;
    if constexpr (ND == 0) {
        static_cast<void>(spec);
        out.law = StateSpace<NXK, NU, NY, T, NWK, NV>{.A = m.law.A, .B = m.law.B, .C = m.law.C, .D = m.law.D, .H = sys.H, .Ts = h};
        if constexpr (Discrete) {
            out.law.G = m.law.G;
        } else {
            out.law.G = Matrix<NXK, NWK, T>{};
            for (size_t i = 0; i < NWK; ++i) {
                out.law.G(i, i) = T{1}; // the plant's, then a late reading's: kalman_sampled's
            }
        }
    } else {
        const auto ch = disturbance_channels(spec, sys, m.law, Ts);
        out.law = disturbance_augment<NWK>(m.law, ch.Bd, ch.Cd, h);
        if constexpr (Discrete) {
            out.law.G.template block<NXE, NW>(0, 0) = m.law.G;
            out.law.G.template block<ND, ND>(NXE, NW) = Matrix<ND, ND, T>::identity();
        } else {
            out.law.G = Matrix<NXK, NXK, T>::identity();
        }
    }
    return out;
}

/**
 * @brief The filter at @p placed on the model of @p sys with @p spec and @p delay: the
 *        FixedGainObserver the Kalman design would run, its innovation gain M placed.
 *
 * The runtime's filter is the current estimator, x̂[k|k] = x̂[k|k−1] + M·(y − C·x̂[k|k−1] − D·u),
 * whose error after an update evolves as e[k|k] = (I − M·C)·A·e[k−1|k−1]. The poles are the
 * eigenvalues of (I − M·C)·A, which are those of the predictor error A − L·C, L = A·M, for any A.
 * The listed ones are the plant's and the disturbance's (s-plane on the sampled form, mapped to
 * z = e^{s·Ts}; z-plane on a discrete plant); the delay states take theirs at 0.
 *
 * An InputDelay's delay states hold the inputs already sent, which predict shifts in: M has no
 * rows on them, as the Kalman gain has none (their covariance is zero), and their error empties
 * in N periods. Their rows of A read no other state but d, which an InputDisturbance adds to u on
 * its way through them: v' = v − G_v·d with G_v = Σ_{j<N} Sʲ·A_vd (S the shift) removes it, so
 * in [c; v'] the error dynamics are block triangular, the core c = [x; late readings; d] on top
 * and the shift S, nilpotent, below. M places the core's (I − M_c·C_c)·A_c with
 * design::place_observer's machinery (place on (Aᵀ, Cᵀ), then A M = L, or the direct placement
 * on a singular A), and its rows on the delay states are G_v·M_d, the correction of d they carry.
 * A late reading's delay states (OutputDelay) are in the core, placed at 0 with the rest:
 * with several readings a pole is placed at most NY times, so two or more delay states per
 * reading are refused with not_placeable then.
 *
 * @return The KalmanResult of the Kalman design on the same model, M and L = A·M placed, Q, R
 *         and P zero (no noise); status invalid_parameter for a pole that is not finite, the
 *         model's, then design::place_observer's (not_placeable for an unobservable mode,
 *         singular for a rank-deficient C)
 */
template<bool Discrete, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NP, typename Spec, typename Delay = NoDelay>
[[nodiscard]] constexpr auto placed_filter_design(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const PlacedPoles<NP, T>&                placed,
    const Spec&                              spec,
    const Delay&                             delay = Delay{}
) {
    constexpr size_t ND = disturbance_traits<Spec, NU, NY, NW>::states;
    static_assert(ND <= NY, "Estimating ND disturbances needs at least ND independent readings (ND <= NY).");
    constexpr size_t NVS = input_delay_states<Delay> * NU;  // the input's delay states, after the plant's
    constexpr size_t NWS = output_delay_states<Delay> * NY; // a late reading's, after those
    constexpr size_t NXK = NX + NVS + NWS + ND;
    constexpr size_t NC = NXK - NVS; // the core the gain is placed on: the plant, the late readings, d
    static_assert(NP == NX + ND, "The poles listed are the plant's and the disturbance's.");
    const auto m = placed_filter_model<Discrete, ND>(sys, Ts, spec, delay);
    using Filter = KalmanResult<NXK, NU, NY, T, placed_noise_inputs<Discrete, ND, NX, NU, NY, NW, Delay>, NV, ND>;
    if (!finite_poles(placed.poles)) {
        return ::damp::detail::refused<Filter>(DesignStatus::invalid_parameter);
    }
    if (m.status != DesignStatus::ok) {
        return ::damp::detail::refused<Filter>(m.status);
    }
    const auto&        A = m.law.A;
    const auto&        C = m.law.C;
    const auto         full = [](size_t c) { return (c < NX) ? c : (c + NVS); }; // a core state's index in the model
    Matrix<NVS, NC, T> Gv{};                                                     // v' = v − Gv·c, nonzero on d's columns alone
    if constexpr (NVS > 0) {
        Matrix<NVS, NVS, T> S{};
        Matrix<NVS, NC, T>  term{};
        for (size_t i = 0; i < NVS; ++i) {
            for (size_t j = 0; j < NVS; ++j) {
                S(i, j) = A(NX + i, NX + j);
            }
            for (size_t c = 0; c < NC; ++c) {
                term(i, c) = A(NX + i, full(c)); // zero but on d's columns: the shift reads only u and d
            }
        }
        Gv = term;
        for (size_t j = 1; j < input_delay_states<Delay>; ++j) {
            term = S * term;
            Gv = Gv + term; // Σ Sʲ·A_vc: S is nilpotent of order N
        }
    }
    StateSpace<NC, NU, NY, T> core{.Ts = m.law.Ts};
    for (size_t r = 0; r < NC; ++r) {
        for (size_t c = 0; c < NC; ++c) {
            T a = A(full(r), full(c));
            for (size_t i = 0; i < NVS; ++i) {
                a += A(full(r), NX + i) * Gv(i, c);
            }
            core.A(r, c) = a;
        }
    }
    for (size_t y = 0; y < NY; ++y) {
        for (size_t c = 0; c < NC; ++c) {
            T v = C(y, full(c));
            for (size_t i = 0; i < NVS; ++i) {
                v += C(y, NX + i) * Gv(i, c);
            }
            core.C(y, c) = v;
        }
    }
    const auto placed_core = place_observer_on(core, with_deadbeat<NWS, !Discrete>(placed.poles, Ts));
    if (!placed_core.success) {
        return ::damp::detail::refused<Filter>(placed_core.status);
    }
    Filter                   out{.sys = m.law};
    const Matrix<NVS, NY, T> Mv = Gv * placed_core.M;
    for (size_t y = 0; y < NY; ++y) {
        for (size_t c = 0; c < NC; ++c) {
            out.M(full(c), y) = placed_core.M(c, y);
        }
        for (size_t i = 0; i < NVS; ++i) {
            out.M(NX + i, y) = Mv(i, y);
        }
    }
    out.L = A * out.M;
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

} // namespace detail
} // namespace design
} // namespace damp
