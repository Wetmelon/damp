// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file dynamics.hpp
 * @brief Evaluating a Dynamics / Measurement model: its Jacobians and one RK4 step over a sample time
 *
 * The shared machinery behind design::linearize, ExtendedKalmanFilter and
 * UnscentedKalmanFilter; the model concepts are in concepts.hpp. A Jacobian comes from the
 * model's own jacobian(x, u) when it has one, else from central differences. The step over a
 * sample time is one classical RK4 step with the input held.
 */


#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/simulation/integrator.hpp"

namespace damp {
namespace detail {

/** @brief Sizes of a model y = g(x, u) whose result and arguments are ColVecs of one scalar; empty otherwise. */
template<typename R, typename X, typename U>
struct model_sizes {};

template<size_t NOUT, size_t NX, size_t NU, typename T>
struct model_sizes<ColVec<NOUT, T>, ColVec<NX, T>, ColVec<NU, T>> {
    static constexpr size_t states = NX;    ///< Rows of x
    static constexpr size_t inputs = NU;    ///< Rows of u
    static constexpr size_t outputs = NOUT; ///< Rows of the result: ẋ for Dynamics, y for Measurement
    using scalar = T;                       ///< Scalar of x, u and the result
};

/** @brief Sizes read from the type of a const call operator g(x, u); empty for any other member. */
template<typename Sig>
struct call_operator_sizes {};

template<typename C, typename R, typename X, typename U>
struct call_operator_sizes<R (C::*)(X, U) const> : model_sizes<damp::remove_cvref_t<R>, damp::remove_cvref_t<X>, damp::remove_cvref_t<U>> {};

template<typename C, typename R, typename X, typename U>
struct call_operator_sizes<R (C::*)(X, U) const noexcept> : model_sizes<damp::remove_cvref_t<R>, damp::remove_cvref_t<X>, damp::remove_cvref_t<U>> {};

/**
 * @brief Sizes of a model type from its call operator: NX, NU and the result's rows
 *
 * The model is a struct (or a lambda) with one non-template operator()(x, u) const taking and
 * returning ColVecs, so the filters read their sizes from the model type and a member reads
 * `ExtendedKalmanFilter<Pendulum, float>`. Empty for an overloaded or template operator().
 */
template<typename F>
struct model_shape {};

template<typename F>
    requires requires { &F::operator(); }
struct model_shape<F> : call_operator_sizes<decltype(&F::operator())> {};

/** @brief A model type whose sizes model_shape can read. */
template<typename F>
concept shaped_model = requires {
    model_shape<F>::states;
    model_shape<F>::inputs;
    model_shape<F>::outputs;
};

/** @brief Default relative central-difference step ∛eps. */
template<typename T>
[[nodiscard]] constexpr T default_jacobian_step() {
    // Central differences: the optimal step ~ ε^(1/3) balances the O(h²) truncation error
    // against round-off, unlike the ε^(1/2) of one-sided differences.
    return damp::cbrt(eps<T>());
}

/**
 * @brief Central-difference Jacobian ∂g/∂p (NOUT × NP) at p: 2·NP evaluations of g.
 *
 * Coordinate i moves by hᵢ = rel·max(|pᵢ|, 1), so the perturbation stays resolvable at large
 * operating points. The divisor is the realized step (pᵢ + hᵢ) − (pᵢ − hᵢ), not 2hᵢ, so the
 * rounding of pᵢ ± hᵢ does not bias the quotient.
 */
template<size_t NOUT, size_t NP, typename T, typename G>
[[nodiscard]] constexpr Matrix<NOUT, NP, T> central_difference(const G& g, const ColVec<NP, T>& p, T rel) {
    Matrix<NOUT, NP, T> J{};
    for (size_t i = 0; i < NP; ++i) {
        const T       mag = damp::abs(p[i]);
        const T       h = rel * ((mag > T{1}) ? mag : T{1});
        ColVec<NP, T> plus = p;
        ColVec<NP, T> minus = p;
        plus[i] += h;
        minus[i] -= h;
        const T               inv_span = T{1} / (plus[i] - minus[i]);
        const ColVec<NOUT, T> g_plus = g(plus);
        const ColVec<NOUT, T> g_minus = g(minus);
        for (size_t r = 0; r < NOUT; ++r) {
            J(r, i) = (g_plus[r] - g_minus[r]) * inv_span;
        }
    }
    return J;
}

/** @brief A model with a jacobian(x, u) member of any return type. */
template<typename M, size_t NX, size_t NU, typename T>
concept declares_jacobian = requires(const M& m, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    m.jacobian(x, u);
};

/** @brief A = ∂f/∂x at (x, u): the model's jacobian(x, u), else 2·NX evaluations of f. */
template<size_t NX, size_t NU, typename T, typename F>
    requires Dynamics<F, NX, NU, T>
[[nodiscard]] constexpr Matrix<NX, NX, T> dynamics_state_jacobian(const F& f, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T rel) {
    if constexpr (AnalyticDynamics<F, NX, NU, T>) {
        const DynamicsJacobian<NX, NU, T> J = f.jacobian(x, u);
        return J.A;
    } else {
        static_assert(!declares_jacobian<F, NX, NU, T>, "a Dynamics model's jacobian(x, u) must return DynamicsJacobian<NX, NU, T>");
        const auto f_of_x = [&](const ColVec<NX, T>& xx) -> ColVec<NX, T> { return f(xx, u); };
        return central_difference<NX>(f_of_x, x, rel);
    }
}

/** @brief A = ∂f/∂x and B = ∂f/∂u at (x, u): the model's jacobian(x, u), else central differences. */
template<size_t NX, size_t NU, typename T, typename F>
    requires Dynamics<F, NX, NU, T>
[[nodiscard]] constexpr DynamicsJacobian<NX, NU, T> dynamics_jacobian(const F& f, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T rel) {
    if constexpr (AnalyticDynamics<F, NX, NU, T>) {
        return f.jacobian(x, u);
    } else {
        const auto f_of_u = [&](const ColVec<NU, T>& uu) -> ColVec<NX, T> { return f(x, uu); };
        return {.A = dynamics_state_jacobian<NX, NU, T>(f, x, u, rel), .B = central_difference<NX>(f_of_u, u, rel)};
    }
}

/** @brief C = ∂h/∂x at (x, u): the model's jacobian(x, u), else 2·NX evaluations of h. */
template<size_t NX, size_t NU, size_t NY, typename T, typename H>
    requires Measurement<H, NX, NU, NY, T>
[[nodiscard]] constexpr Matrix<NY, NX, T> measurement_state_jacobian(const H& h, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T rel) {
    if constexpr (AnalyticMeasurement<H, NX, NU, NY, T>) {
        const MeasurementJacobian<NX, NU, NY, T> J = h.jacobian(x, u);
        return J.C;
    } else {
        static_assert(!declares_jacobian<H, NX, NU, T>, "a Measurement model's jacobian(x, u) must return MeasurementJacobian<NX, NU, NY, T>");
        const auto h_of_x = [&](const ColVec<NX, T>& xx) -> ColVec<NY, T> { return h(xx, u); };
        return central_difference<NY>(h_of_x, x, rel);
    }
}

/** @brief C = ∂h/∂x and D = ∂h/∂u at (x, u): the model's jacobian(x, u), else central differences. */
template<size_t NX, size_t NU, size_t NY, typename T, typename H>
    requires Measurement<H, NX, NU, NY, T>
[[nodiscard]] constexpr MeasurementJacobian<NX, NU, NY, T> measurement_jacobian(const H& h, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T rel) {
    if constexpr (AnalyticMeasurement<H, NX, NU, NY, T>) {
        return h.jacobian(x, u);
    } else {
        const auto h_of_u = [&](const ColVec<NU, T>& uu) -> ColVec<NY, T> { return h(x, uu); };
        return {.C = measurement_state_jacobian<NX, NU, NY, T>(h, x, u, rel), .D = central_difference<NY>(h_of_u, u, rel)};
    }
}

/** @brief x(t + Ts) of ẋ = f(x, u) with u held: one classical RK4 step, four evaluations of f. */
template<size_t NX, size_t NU, typename T, typename F>
    requires Dynamics<F, NX, NU, T>
[[nodiscard]] constexpr ColVec<NX, T> rk4_step(const F& f, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T Ts) {
    const auto rhs = [&](T /*t*/, const ColVec<NX, T>& xx) -> ColVec<NX, T> { return f(xx, u); };
    return sim::RK4<NX, T>{}.evolve(rhs, x, T{0}, Ts).x;
}

/**
 * @brief Φ = I + M + M²/2 + M³/6 + M⁴/24 with M = A·Ts: the Jacobian of rk4_step with A held over the step.
 *
 * The RK4 stability function R(M), evaluated in Horner form with three matrix products. It is
 * the exact Jacobian of one RK4 step of a linear model.
 *
 * @see "Solving Ordinary Differential Equations II" (Hairer & Wanner, 1996), §IV.2
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr Matrix<NX, NX, T> rk4_transition(const Matrix<NX, NX, T>& A, T Ts) {
    const Matrix<NX, NX, T> I = Matrix<NX, NX, T>::identity();
    const Matrix<NX, NX, T> M = A * Ts;
    Matrix<NX, NX, T>       R = I + (M * (T{1} / T{4}));
    R = I + ((M * R) * (T{1} / T{3}));
    R = I + ((M * R) * T{0.5});
    return I + (M * R);
}

/**
 * @brief A nonlinear filter's sample time: refused under constant evaluation unless finite and
 *        positive; at run time a bad one becomes 0, so predict holds the state.
 */
template<typename T>
[[nodiscard]] constexpr T filter_sample_time(T Ts) {
    if (damp::finite_positive(Ts)) {
        return Ts;
    }
    if (damp::is_constant_evaluated()) {
        refuse_failed_design(); // constinit / constexpr: a filter without a valid Ts does not compile
    }
    return T{0};
}

} // namespace detail
} // namespace damp
