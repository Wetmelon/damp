// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/// G2.5: design in double, `.as<float>()`, tick in float. A successful cart-pole
/// (4-state) and pendulum (2-state) LQR stays stabilizing after the conversion,
/// with an honest scaled DARE residual — both `static_assert` and runtime.

#include "damp/backend.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/stability.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "examples/control/cart_pole/cart_pole_controller.hpp"
#include "examples/control/pendulum/pendulum_controller.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

/// DARE residual
///   AᵀSA − S − (AᵀSB + N) (R + BᵀSB)⁻¹ (BᵀSA + Nᵀ) + Q
/// via mat::solve, not an explicit inverse. N defaults to zero.
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr optional<Matrix<NX, NX, T>> dare_residual(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NX, T>& S,
    const Matrix<NX, NU, T>& N = {}
) {
    const auto solved = mat::solve(R + (B.t() * S * B), (B.t() * S * A) + N.t());
    if (!solved) {
        return nullopt;
    }
    return (A.t() * S * A) - S - (((A.t() * S * B) + N) * solved.value()) + Q;
}

/// Scale for an honest residual bound: max(1, ‖S‖∞, ‖Q‖∞, ‖A‖∞² ‖S‖∞).
/// Bound is n · ε · scale so it tracks size and magnitude, not units.
template<size_t NX, typename T>
[[nodiscard]] constexpr T dare_residual_bound(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NX, NX, T>& S
) {
    const T nS = mat::infinity_norm(S);
    const T nQ = mat::infinity_norm(Q);
    const T nA = mat::infinity_norm(A);
    const T scale = damp::max(T{1}, damp::max(nS, damp::max(nQ, (nA * nA) * nS)));
    return static_cast<T>(NX) * eps<T>() * scale;
}

template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr bool dare_residual_ok(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NX, T>& S,
    const Matrix<NX, NU, T>& N = {}
) {
    const auto resid = dare_residual(A, B, Q, R, S, N);
    if (!resid) {
        return false;
    }
    return mat::infinity_norm(*resid) <= dare_residual_bound(A, Q, S);
}

template<size_t NX, typename T>
[[nodiscard]] constexpr bool strictly_schur(const Matrix<NX, NX, T>& Acl) {
    const auto ev = mat::compute_eigenvalues(Acl);
    if (!ev.success) {
        return false;
    }
    for (size_t i = 0; i < NX; ++i) {
        if (!(ev.values[i].abs() < T{1})) {
            return false;
        }
    }
    return true;
}

// ---- cart-pole 4-state (example header; do not copy A, B) ----

namespace cart = examples_cart_pole;

static_assert(cart::lqr_d.success && cart::lqr_d.is_stable());
constexpr auto cart_lqr_f = cart::lqr_d.as<float>();
static_assert(cart_lqr_f.success && cart_lqr_f.is_stable());

constexpr auto cart_zoh_opt = discretize(cart::sys, cart::Ts, DiscretizationMethod::ZOH);
static_assert(cart_zoh_opt.has_value());
constexpr auto cart_zoh = *cart_zoh_opt;
// Van Loan (Qd, Rd, Nd) is the DARE that S actually solves — not the continuous Q, R.
constexpr auto cart_cost = design::discretize_lqr_cost(cart::sys.A, cart::sys.B, cart::Q, Matrix<1, 1>(cart::R), cart::Ts).value();
constexpr auto cart_Ad_f = cart_zoh.A.as<float>();
constexpr auto cart_Bd_f = cart_zoh.B.as<float>();
constexpr auto cart_Qd_f = cart_cost.Q.as<float>();
constexpr auto cart_Rd_f = cart_cost.R.as<float>();
constexpr auto cart_Nd_f = cart_cost.N.as<float>();
static_assert(dare_residual_ok(cart_Ad_f, cart_Bd_f, cart_Qd_f, cart_Rd_f, cart_lqr_f.S, cart_Nd_f));

// ---- pendulum 2-state ----

namespace pend = examples_pendulum;

static_assert(pend::lqr_d.success && pend::lqr_d.is_stable());
constexpr auto pend_lqr_f = pend::lqr_d.as<float>();
static_assert(pend_lqr_f.success && pend_lqr_f.is_stable());

constexpr auto pend_zoh_opt = discretize(pend::sys, pend::Ts, DiscretizationMethod::ZOH);
static_assert(pend_zoh_opt.has_value());
constexpr auto pend_zoh = *pend_zoh_opt;
constexpr auto pend_cost = design::discretize_lqr_cost(pend::sys.A, pend::sys.B, pend::Q, Matrix<1, 1>(pend::R), pend::Ts).value();
constexpr auto pend_Ad_f = pend_zoh.A.as<float>();
constexpr auto pend_Bd_f = pend_zoh.B.as<float>();
constexpr auto pend_Qd_f = pend_cost.Q.as<float>();
constexpr auto pend_Rd_f = pend_cost.R.as<float>();
constexpr auto pend_Nd_f = pend_cost.N.as<float>();
static_assert(dare_residual_ok(pend_Ad_f, pend_Bd_f, pend_Qd_f, pend_Rd_f, pend_lqr_f.S, pend_Nd_f));

// ---- Kalman 2-state (existing golden plant; larger than the 1×1 as<float> round-trip) ----

constexpr auto kf_d = design::kalman(
    StateSpace<2, 1, 1, double, 2, 1>{
        Matrix<2, 2>{{0.9, 0.1}, {0.0, 0.8}},
        Matrix<2, 1>{{0.1}, {0.2}},
        Matrix<1, 2>{{1.0, 0.0}},
        Matrix<1, 1>::zeros(),
        Matrix<2, 2>::identity(),
        Matrix<1, 1>::identity(),
        1.0
    },
    {.process = Matrix<2, 2>{{0.1, 0.0}, {0.0, 0.1}}, .measurement = 0.5}
);
static_assert(kf_d.success);
constexpr auto kf_f = kf_d.as<float>();
static_assert(kf_f.success);
// Filter DARE is the control DARE on the dual (Aᵀ, Cᵀ, Q_eff, R_eff).
constexpr auto kf_At_f = kf_f.sys.A.transpose();
constexpr auto kf_Ct_f = kf_f.sys.C.transpose();
constexpr auto kf_Qeff_f = kf_f.sys.G * kf_f.Q * kf_f.sys.G.transpose();
constexpr auto kf_Reff_f = kf_f.sys.H * kf_f.R * kf_f.sys.H.transpose();
static_assert(dare_residual_ok(kf_At_f, kf_Ct_f, kf_Qeff_f, kf_Reff_f, kf_f.P));

} // namespace

TEST_SUITE("Float deploy") {
    TEST_CASE("cart-pole LQR as<float> stays stabilizing with a scaled DARE residual") {
        CHECK(cart::lqr_d.success);
        CHECK(cart::lqr_d.is_stable());
        CHECK(cart_lqr_f.success);
        CHECK(cart_lqr_f.is_stable());

        const auto poles = stability::closed_loop_poles(cart_Ad_f, cart_Bd_f, cart_lqr_f.K);
        CHECK(strictly_schur(cart_Ad_f - (cart_Bd_f * cart_lqr_f.K)));
        for (size_t i = 0; i < 4; ++i) {
            CHECK(poles[i].abs() < 1.0f);
        }

        const auto resid = dare_residual(cart_Ad_f, cart_Bd_f, cart_Qd_f, cart_Rd_f, cart_lqr_f.S, cart_Nd_f);
        REQUIRE(resid.has_value());
        const float bound = dare_residual_bound(cart_Ad_f, cart_Qd_f, cart_lqr_f.S);
        CHECK(mat::infinity_norm(*resid) <= bound);
    }

    TEST_CASE("pendulum LQR as<float> stays stabilizing with a scaled DARE residual") {
        CHECK(pend::lqr_d.success);
        CHECK(pend::lqr_d.is_stable());
        CHECK(pend_lqr_f.success);
        CHECK(pend_lqr_f.is_stable());

        const auto poles = stability::closed_loop_poles(pend_Ad_f, pend_Bd_f, pend_lqr_f.K);
        CHECK(strictly_schur(pend_Ad_f - (pend_Bd_f * pend_lqr_f.K)));
        for (size_t i = 0; i < 2; ++i) {
            CHECK(poles[i].abs() < 1.0f);
        }

        const auto resid = dare_residual(pend_Ad_f, pend_Bd_f, pend_Qd_f, pend_Rd_f, pend_lqr_f.S, pend_Nd_f);
        REQUIRE(resid.has_value());
        const float bound = dare_residual_bound(pend_Ad_f, pend_Qd_f, pend_lqr_f.S);
        CHECK(mat::infinity_norm(*resid) <= bound);
    }

    TEST_CASE("2-state Kalman as<float> keeps Schur-stable error dynamics") {
        CHECK(kf_d.success);
        CHECK(kf_f.success);

        const auto Acl = kf_f.sys.A - (kf_f.L * kf_f.sys.C);
        CHECK(strictly_schur(Acl));
        // Dual closed-loop poles of (Aᵀ − Cᵀ Lᵀ) = eig(A − LC).
        const auto poles = stability::closed_loop_poles(kf_At_f, kf_Ct_f, kf_f.L.transpose());
        for (size_t i = 0; i < 2; ++i) {
            CHECK(poles[i].abs() < 1.0f);
        }

        const auto resid = dare_residual(kf_At_f, kf_Ct_f, kf_Qeff_f, kf_Reff_f, kf_f.P);
        REQUIRE(resid.has_value());
        const float bound = dare_residual_bound(kf_At_f, kf_Qeff_f, kf_f.P);
        CHECK(mat::infinity_norm(*resid) <= bound);
    }
}
