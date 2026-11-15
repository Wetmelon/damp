// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lq_cost.hpp
 * @brief The discrete LQ problem behind a gain: the integrator-augmented model, and the continuous
 *        cost over one sample on it, with integrators (sampled_integral_cost) and with an input
 *        that arrives late (delayed_cost)
 */

#include "damp/backend.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/riccati.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {
namespace detail {

/** @brief The discrete model of [x; xi] under the regulator (r = 0): xi[k+1] = xi[k] − h·Ci·(C·x[k] + Dy·u[k]). */
template<size_t NX, size_t NU, size_t NI, typename T>
struct IntegralModel {
    Matrix<NX + NI, NX + NI, T> A{}; ///< [A 0; −h·Ci·C I]
    Matrix<NX + NI, NU, T>      B{}; ///< [B; −h·Ci·Dy]
};

/** @brief Augment a discrete model with the integrators of the selected outputs (see IntegralModel). */
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T>
[[nodiscard]] constexpr IntegralModel<NX, NU, NI, T> integral_model(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NY, NX, T>& C,
    const Matrix<NY, NU, T>& Dy,
    const Matrix<NI, NY, T>& Ci,
    T                        h
) {
    IntegralModel<NX, NU, NI, T> m{};
    m.A.template block<NX, NX>(0, 0) = A;
    m.A.template block<NI, NX>(NX, 0) = -(Ci * C) * h;
    m.A.template block<NI, NI>(NX, NX) = Matrix<NI, NI, T>::identity(); // integrator memory
    m.B.template block<NX, NU>(0, 0) = B;
    m.B.template block<NI, NU>(NX, 0) = -(Ci * Dy) * h;
    return m;
}

/**
 * @brief The continuous cost over one sample on [x; xi], for the runtime's integrator
 *
 * Between samples the plant runs continuously on the held input and the integrator ramps on the
 * error sampled at the start, ẋi = −Ci·(Cm·x_k + Dm·u). The Van Loan model is [x; xi; x_k], x_k
 * held over the period; folding x_k onto x gives the exact discrete weights of
 * ∫ ([x; xi]ᵀQ[x; xi] + 2[x; xi]ᵀNu + uᵀRu) dt on [x; xi] and u.
 */
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T>
[[nodiscard]] constexpr damp::optional<LQRCost<NX + NI, NU, T>> sampled_integral_cost(
    const Matrix<NX, NX, T>&       A,
    const Matrix<NX, NU, T>&       B,
    const Matrix<NY, NX, T>&       Cm,
    const Matrix<NY, NU, T>&       Dm,
    const Matrix<NI, NY, T>&       Ci,
    const LQRCost<NX + NI, NU, T>& weights,
    T                              Ts
) {
    constexpr size_t NZ = NX + NI; // [x; xi], the weighted part
    constexpr size_t NS = NZ + NX; // [x; xi; x_k]
    const auto       CiC = Ci * Cm;

    Matrix<NS, NS, T> As{};
    Matrix<NS, NU, T> Bs{};
    As.template block<NX, NX>(0, 0) = A;
    for (size_t i = 0; i < NI; ++i) { // the integrator sees the sampled output
        for (size_t j = 0; j < NX; ++j) {
            As(NX + i, NZ + j) = -CiC(i, j);
        }
    }
    Bs.template block<NX, NU>(0, 0) = B;
    Bs.template block<NI, NU>(NX, 0) = -(Ci * Dm);
    Matrix<NS, NS, T> Qs{};
    Matrix<NS, NU, T> Ns{};
    Qs.template block<NZ, NZ>(0, 0) = weights.Q;
    Ns.template block<NZ, NU>(0, 0) = weights.N;
    const auto cost_opt = discretize_lqr_cost(As, Bs, Qs, weights.R, Ts, Ns);
    if (!cost_opt) {
        return damp::nullopt;
    }
    const auto& cost = *cost_opt;

    //! [x; xi; x_k] = P·[x; xi]: x_k is x at the sample
    Matrix<NS, NZ, T> P{};
    P.template block<NZ, NZ>(0, 0) = Matrix<NZ, NZ, T>::identity();
    P.template block<NX, NX>(NZ, 0) = Matrix<NX, NX, T>::identity();
    const Matrix<NZ, NS, T> Pt = P.t();
    const Matrix<NZ, NZ, T> Qf = Pt * cost.Q * P;
    const Matrix<NZ, NZ, T> Qft = Qf.t();
    return LQRCost<NZ, NU, T>{
        .Q = Matrix<NZ, NZ, T>((Qf + Qft) * static_cast<T>(0.5)),
        .R = cost.R,
        .N = Matrix<NZ, NU, T>(Pt * cost.N),
    };
}

/** @brief The weights on [x; xi]: Q and Q_integral on the diagonal, N on the plant's rows. */
template<size_t NX, size_t NU, typename T, size_t NI>
[[nodiscard]] constexpr LQRCost<NX + NI, NU, T> stacked_weights(const LQRConfig<NX, NU, NI, T>& w) {
    LQRCost<NX + NI, NU, T> out{.R = w.R};
    out.Q.template block<NX, NX>(0, 0) = w.Q;
    out.Q.template block<NI, NI>(NX, NX) = w.Q_integral;
    out.N.template block<NX, NU>(0, 0) = w.N;
    return out;
}

/**
 * @brief The continuous cost over one period of a plant whose input arrives late, on the law's
 *        state [z; ξ], z = [x; v(k−1); …], and the command u: exact under ZOH.
 *
 * Within the period input j holds v_j(k − d_j − 1) until its switch at τ_j, then v_j(k − d_j)
 * (discretize with an InputDelay), and the integrators ramp on the error sampled at the start,
 * ξ̇ = −Ci·y_k. The switches cut the period into stretches over which every input holds one
 * value. Over each, Van Loan gives the exact cost of [x; ξ] from the stretch's start under the
 * held values, and the exponential carries the start to the next stretch. The start and the held
 * values are linear maps of [z; ξ; u], so the stretches sum to a quadratic form on it.
 *
 * The input cost is booked on the command when it is sent, Ts·uᵀRu, not over the periods it
 * reaches the plant in. A held input acts for Ts in all, so over any trajectory the two
 * bookings differ only by the cost of the inputs already in the registers at the start, which
 * no command changes: the same gain, with R on the command where the DARE needs it. Q and the
 * cross term N act on the plant as it runs, x(t) and the input it holds.
 */
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T, size_t NW, size_t NV, size_t NZ>
[[nodiscard]] constexpr damp::optional<LQRCost<NZ + NI, NU, T>> delayed_cost(
    const StateSpace<NX, NU, NY, T, NW, NV>&                 sys,   // continuous
    const StateSpace<NZ, NU, NY, T, NW, NV>&                 model, // discretize(sys, Ts, InputDelay<N>{θ})
    const damp::array<::damp::detail::RegisterDelay<T>, NU>& inputs,
    const Matrix<NI, NY, T>&                                 Ci,
    const LQRConfig<NX, NU, NI, T>&                          w,
    T                                                        Ts
) {
    constexpr size_t NL = NZ + NI;      // the law's state [z; ξ]
    constexpr size_t NQ = NL + NU;      // and the command
    constexpr size_t NE = NX + NI;      // what runs over the period: [x; ξ]
    constexpr size_t NA = NE + NU + NI; // and what a stretch holds: [v; e]

    //! d/dt [x; ξ; v; e] = [A·x + B·v; −e; 0; 0], weighted by [Q 0 N 0; 0 Q_integral 0 0; Nᵀ 0 0 0; 0 0 0 0]
    Matrix<NA, NA, T> F{};
    Matrix<NA, NA, T> S{};
    F.template block<NX, NX>(0, 0) = sys.A;
    F.template block<NX, NU>(0, NE) = sys.B;
    S.template block<NX, NX>(0, 0) = w.Q;
    S.template block<NX, NU>(0, NE) = w.N;
    S.template block<NU, NX>(NE, 0) = w.N.t();
    for (size_t i = 0; i < NI; ++i) {
        F(NX + i, NE + NU + i) = T{-1};
        for (size_t c = 0; c < NI; ++c) {
            S(NX + i, NX + c) = w.Q_integral(i, c);
        }
    }

    //! The stretches: 0, each switch once in time order, then Ts.
    damp::array<T, NU + 2> edge{};
    size_t                 cuts = 1;
    for (size_t j = 0; j < NU; ++j) {
        const T tau = inputs[j].fraction;
        bool    seen = !(tau > T{0});
        for (size_t c = 1; c < cuts; ++c) {
            seen = seen || (edge[c] == tau);
        }
        if (!seen) {
            size_t c = cuts++;
            for (; c > 1 && edge[c - 1] > tau; --c) {
                edge[c] = edge[c - 1];
            }
            edge[c] = tau;
        }
    }
    edge[cuts] = Ts;

    //! [x; ξ] at the start of the period as a map of [z; ξ; u], and the error held over it, e = Ci·(C·z + D·u)
    Matrix<NE, NQ, T> start{};
    for (size_t i = 0; i < NX; ++i) {
        start(i, i) = T{1};
    }
    for (size_t i = 0; i < NI; ++i) {
        start(NX + i, NZ + i) = T{1};
    }
    const Matrix<NI, NZ, T> CiC = Ci * model.C;
    const Matrix<NI, NU, T> CiD = Ci * model.D;

    Matrix<NQ, NQ, T> total{};
    for (size_t s = 0; s < cuts; ++s) {
        const T           a = edge[s];
        const T           h = edge[s + 1] - a;
        Matrix<NA, NQ, T> M{}; // [x; ξ; v; e] at the stretch's start
        M.template block<NE, NQ>(0, 0) = start;
        for (size_t j = 0; j < NU; ++j) {
            const bool   before = inputs[j].fraction > T{0} && a < inputs[j].fraction;
            const size_t age = inputs[j].whole + (before ? 1 : 0);
            M(NE + j, (age == 0) ? NL + j : NX + ((age - 1) * NU) + j) = T{1};
        }
        for (size_t i = 0; i < NI; ++i) {
            for (size_t c = 0; c < NZ; ++c) {
                M(NE + NU + i, c) = CiC(i, c);
            }
            for (size_t j = 0; j < NU; ++j) {
                M(NE + NU + i, NL + j) = CiD(i, j);
            }
        }
        const auto W = ::damp::detail::van_loan_integral(F, S, h);
        const auto E = mat::expm(Matrix<NA, NA, T>(F * h));
        if (!W || !E) {
            return damp::nullopt;
        }
        const Matrix<NQ, NA, T> Mt = M.t();
        total = total + ((Mt * *W) * M);
        start = E->template block<NE, NA>(0, 0) * M;
    }
    total.template block<NU, NU>(NL, NL) = total.template block<NU, NU>(NL, NL) + (w.R * Ts);
    total = (total + total.t()) * static_cast<T>(0.5);
    if (!mat::all_of(total, [](T v) { return damp::isfinite(v); })) {
        return damp::nullopt;
    }
    return LQRCost<NL, NU, T>{
        .Q = total.template block<NL, NL>(0, 0),
        .R = total.template block<NU, NU>(NL, NL),
        .N = total.template block<NL, NU>(0, NL),
    };
}

} // namespace detail
} // namespace design
} // namespace damp
