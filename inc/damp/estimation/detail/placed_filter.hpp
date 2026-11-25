// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private: a public header that exports it provides its names

/**
 * @file placed_filter.hpp
 * @brief The placed filter: the FixedGainObserver the Kalman design would run on the same model,
 *        its gain placed at a list of error poles (placed_filter_design), for design::observer and
 *        the filters design::state_feedback and design::mpc place at ObserverPoles
 */

#include "damp/backend.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/options.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp" // KalmanResult; the model's Domain, sampled_model and delays come with it
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {
namespace detail {

/** @brief The scalar of a pole as written: a floating-point number is its own. */
template<typename P>
struct pole_scalar {
    using type = P; ///< The number's type
};

/** @brief A damp::complex: its parts'. */
template<typename U>
struct pole_scalar<damp::complex<U>> {
    using type = U; ///< The parts' type
};

/** @brief A list of poles in the scalar @p T, as the placement takes them. */
template<typename T, size_t N, typename U>
[[nodiscard]] constexpr damp::array<damp::complex<T>, N> pole_list(const damp::complex<U> (&poles)[N]) {
    damp::array<damp::complex<T>, N> out{};
    for (size_t i = 0; i < N; ++i) {
        out[i] = damp::complex<T>(poles[i]);
    }
    return out;
}

/** @brief The poles of a placed gain or filter, in the plant's scalar: the gain's or the filter's, and the design adds 0 for the delay states. */
template<size_t N, typename T>
struct PlacedPoles {
    damp::array<damp::complex<T>, N> poles{}; ///< The poles listed, in the time domain of the plant
};

/**
 * @brief A filter's poles list the count it places: the plant's states plus the disturbance
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

/** @brief The listed poles mapped into the model's domain, z = e^{s·Ts} on a sampled one (@p D sampled), then 0 for the @p NZ states placed deadbeat. */
template<size_t NZ, Domain D, size_t NP, typename T>
[[nodiscard]] constexpr damp::array<damp::complex<T>, NP + NZ> with_deadbeat(const damp::array<damp::complex<T>, NP>& listed, T Ts) {
    const auto                             mapped = (D == Domain::sampled) ? poles_s_to_z(listed, Ts) : listed;
    damp::array<damp::complex<T>, NP + NZ> out{};
    for (size_t i = 0; i < NP; ++i) {
        out[i] = mapped[i];
    }
    return out;
}

/** @brief The noise inputs of the filter's model as kalman_design gives them: the plant's or its states, then d's or none. */
template<Domain D, size_t ND, size_t NX, size_t NU, size_t NY, size_t NW, typename Delay>
inline constexpr size_t placed_noise_inputs = (ND > 0) ? ((D == Domain::plant) ? (NW + ND) : (NX + delay_states<Delay, NU, NY> + ND))
                                                       : ((D == Domain::plant) ? NW : ((output_delay_states<Delay> > 0) ? (NX + delay_states<Delay, NU, NY>) : NX));

/**
 * @brief The filter's model a design on @p sys with @p spec and @p delay estimates on, as
 *        kalman_design builds it: the plant, its delay states, then the disturbance's, with the noise
 *        inputs G the Kalman design gives it (none of the noise is read here). On readings that
 *        arrive late, a spec's d reaches each as it was when the reading was taken (on_late_readings).
 */
template<Domain D, size_t ND, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Spec, typename Delay>
[[nodiscard]] constexpr auto placed_filter_model(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const Spec& spec, const Delay& delay) {
    constexpr size_t                      NXE = NX + delay_states<Delay, NU, NY>;
    constexpr size_t                      NXK = NXE + ND;
    constexpr size_t                      NWK = placed_noise_inputs<D, ND, NX, NU, NY, NW, Delay>;
    SampledModel<NXK, NU, NY, T, NWK, NV> out{};
    const auto                            m = sampled_model(sys, Ts, D, delay);
    out.status = m.status;
    if (m.status != DesignStatus::ok) {
        return out;
    }
    const T h = (D == Domain::plant) ? sys.Ts : Ts;
    if constexpr (ND == 0) {
        static_cast<void>(spec);
        out.law = StateSpace<NXK, NU, NY, T, NWK, NV>{.A = m.law.A, .B = m.law.B, .C = m.law.C, .D = m.law.D, .H = sys.H, .Ts = h};
        if constexpr (D == Domain::plant) {
            out.law.G = m.law.G;
        } else {
            out.law.G = Matrix<NXK, NWK, T>{};
            for (size_t i = 0; i < NWK; ++i) {
                out.law.G(i, i) = T{1}; // the plant's, then a late reading's: kalman_sampled's
            }
        }
    } else {
        const auto ch = on_late_readings(spec, disturbance_channels(spec, sys, m.law, Ts), sys, Ts, delay); // kalman_design's channels
        out.law = disturbance_augment<NWK>(m.law, ch.Bd, ch.Cd, h);
        if constexpr (D == Domain::plant) {
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
template<Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NP, typename Spec, typename Delay = NoDelay>
[[nodiscard]] constexpr auto placed_filter_design(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const PlacedPoles<NP, T>&                placed,
    const Spec&                              spec,
    const Delay&                             delay = Delay{}
) {
    constexpr size_t ND = disturbance_states<Spec, NU, NY, NW>;
    static_assert(ND <= NY, "Estimating ND disturbances needs at least ND independent readings (ND <= NY).");
    constexpr size_t NVS = input_delay_states<Delay> * NU;  // the input's delay states, after the plant's
    constexpr size_t NWS = output_delay_states<Delay> * NY; // a late reading's, after those
    constexpr size_t NXK = NX + NVS + NWS + ND;
    constexpr size_t NC = NXK - NVS; // the core the gain is placed on: the plant, the late readings, d
    static_assert(NP == NX + ND, "The poles listed are the plant's and the disturbance's.");
    const auto m = placed_filter_model<D, ND>(sys, Ts, spec, delay);
    using Filter = KalmanResult<NXK, NU, NY, T, placed_noise_inputs<D, ND, NX, NU, NY, NW, Delay>, NV, ND>;
    if (!finite_poles(placed.poles)) {
        return refused<Filter>(DesignStatus::invalid_parameter);
    }
    if (m.status != DesignStatus::ok) {
        return refused<Filter>(m.status);
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
    const auto placed_core = place_observer_on(core, with_deadbeat<NWS, D>(placed.poles, Ts));
    if (!placed_core.success) {
        return refused<Filter>(placed_core.status);
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
    return set_status(out, DesignStatus::ok);
}

} // namespace detail
} // namespace design
} // namespace damp
