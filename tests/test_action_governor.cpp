// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>

#include "damp/controllers/action_governor.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/toolbox/bounds.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// =============================================================================
// Box projection
// =============================================================================

TEST_CASE("project_box is identity when already safe") {
    const ColVec<2> u{0.5, -0.5};
    const ColVec<2> umin{-1.0, -1.0};
    const ColVec<2> umax{1.0, 1.0};

    const auto out = project_box(u, umin, umax);
    CHECK(out(0) == doctest::Approx(0.5));
    CHECK(out(1) == doctest::Approx(-0.5));
}

TEST_CASE("project_box clamps violated channels") {
    const ColVec<3> u{-2.0, 0.5, 9.0};
    const ColVec<3> umin{-1.0, 0.0, -5.0};
    const ColVec<3> umax{1.0, 2.0, 5.0};

    const auto out = project_box(u, umin, umax);
    CHECK(out(0) == doctest::Approx(-1.0));
    CHECK(out(1) == doctest::Approx(0.5));
    CHECK(out(2) == doctest::Approx(5.0));
}

TEST_CASE("project_box with Bounds matches ColVec form") {
    const Bounds<2> b{{-1.0, -2.0}, {1.0, 2.0}};
    const ColVec<2> u{1.5, -3.0};
    const auto      out = project_box(u, b);
    CHECK(out(0) == doctest::Approx(1.0));
    CHECK(out(1) == doctest::Approx(-2.0));
}

TEST_CASE("BoxCommandFilter filter and operator() match project_box") {
    const Bounds<2, float>     b{{-1.f, -2.f}, {1.f, 2.f}};
    BoxCommandFilter<2, float> filt{b};

    const ColVec<2, float> u_des{1.5f, -3.f};
    const auto             u1 = filt.filter(u_des);
    const auto             u2 = filt(u_des);
    CHECK(u1(0) == doctest::Approx(1.f));
    CHECK(u1(1) == doctest::Approx(-2.f));
    CHECK(u2(0) == u1(0));
    CHECK(u2(1) == u1(1));
    CHECK(filt.is_safe(u1));
    CHECK_FALSE(filt.is_safe(u_des));
}

TEST_CASE("BoxCommandFilter is constexpr") {
    constexpr Bounds<2>                   b{{-1.0, -2.0}, {1.0, 2.0}};
    constexpr BoxCommandFilter<2, double> filt{b};
    constexpr ColVec<2>                   u_des{1.5, -3.0};
    constexpr auto                        u = filt.filter(u_des);
    static_assert(u(0) == 1.0);
    static_assert(u(1) == -2.0);
    static_assert(filt.is_safe(u));
    static_assert(!filt.is_safe(u_des));

    constexpr auto u_id = project_box(ColVec<2>{0.0, 0.0}, b);
    static_assert(u_id(0) == 0.0);
    static_assert(u_id(1) == 0.0);
}

TEST_CASE("BoxCommandFilter::as converts bounds") {
    const BoxCommandFilter<1, double> d{Bounds<1, double>{-2.0, 2.0}};
    const auto                        f = d.as<float>();
    const auto                        u = f.filter(ColVec<1, float>{3.f});
    CHECK(u(0) == doctest::Approx(2.f));
}

// =============================================================================
// Affine / QP projection
// =============================================================================

TEST_CASE("project_affine is identity when already safe") {
    // Half-space u0 + u1 ≤ 2; desired (0.5, 0.5) is interior.
    const Matrix<1, 2> A{{1.0, 1.0}};
    const ColVec<1>    b{2.0};
    const ColVec<2>    u_des{0.5, 0.5};

    const auto res = design::project_affine(u_des, A, b);
    REQUIRE(res.success);
    CHECK_FALSE(res.modified);
    CHECK(res.u(0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(res.u(1) == doctest::Approx(0.5).epsilon(1e-9));
}

TEST_CASE("project_affine projects onto a half-space boundary") {
    // min ‖u − (2,2)‖² s.t. u0 + u1 ≤ 1 → nearest point on the line is (0.5, 0.5).
    const Matrix<1, 2> A{{1.0, 1.0}};
    const ColVec<1>    b{1.0};
    const ColVec<2>    u_des{2.0, 2.0};

    const auto res = design::project_affine(u_des, A, b);
    REQUIRE(res.success);
    CHECK(res.modified);
    CHECK(res.u(0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(res.u(1) == doctest::Approx(0.5).epsilon(1e-9));

    // Feasibility: A u ≤ b
    CHECK(res.u(0) + res.u(1) == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("project_affine projects onto two half-spaces") {
    // u0 ≤ 1, u1 ≤ 1; desired (2, 3) → (1, 1)
    const Matrix<2, 2> A{{1.0, 0.0}, {0.0, 1.0}};
    const ColVec<2>    b{1.0, 1.0};
    const ColVec<2>    u_des{2.0, 3.0};

    const auto res = design::project_affine(u_des, A, b);
    REQUIRE(res.success);
    CHECK(res.modified);
    CHECK(res.u(0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(res.u(1) == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("project_affine reports infeasible") {
    // u ≤ -1 and -u ≤ -1 (u ≥ 1) cannot both hold.
    const Matrix<2, 1> A{{1.0}, {-1.0}};
    const ColVec<2>    b{-1.0, -1.0};
    const ColVec<1>    u_des{0.0};

    const auto res = design::project_affine(u_des, A, b);
    CHECK_FALSE(res.success);
    CHECK(res.status == design::QPStatus::Infeasible);
    CHECK(res.u(0) == doctest::Approx(0.0)); // fallback = u_des
}

TEST_CASE("ActionGovernor operator is fail-closed") {
    // Feasible first: u ≤ 1; then infeasible pair so operator keeps last safe u.
    ActionGovernor<1, 1, double> ok_filt{Matrix<1, 1>{{1.0}}, ColVec<1>{1.0}};
    const auto                   u_ok = ok_filt(ColVec<1>{0.25});
    CHECK(u_ok(0) == doctest::Approx(0.25));

    ActionGovernor<1, 2, double> bad{Matrix<2, 1>{{1.0}, {-1.0}}, ColVec<2>{-1.0, -1.0}};
    // No prior success → zero (not the wild u_des)
    const auto u0 = bad(ColVec<1>{5.0});
    CHECK(u0(0) == doctest::Approx(0.0));
    CHECK_FALSE(bad.filter(ColVec<1>{5.0}).success);
}

TEST_CASE("project_affine is constexpr") {
    constexpr Matrix<1, 2> A{{1.0, 1.0}};
    constexpr ColVec<1>    b{1.0};
    constexpr ColVec<2>    u_des{2.0, 2.0};
    constexpr auto         res = design::project_affine(u_des, A, b);
    static_assert(res.success);
    static_assert(res.modified);
    static_assert(res.u(0) - 0.5 < 1e-9 && res.u(0) - 0.5 > -1e-9);
    static_assert(res.u(1) - 0.5 < 1e-9 && res.u(1) - 0.5 > -1e-9);
}

TEST_CASE("ActionGovernor runtime matches design::project_affine") {
    const Matrix<1, 2, float>   A{{1.f, 1.f}};
    const ColVec<1, float>      b{1.f};
    ActionGovernor<2, 1, float> filt{A, b};

    const ColVec<2, float> u_des{2.f, 2.f};
    const auto             res = filt.filter(u_des);
    REQUIRE(res.success);
    CHECK(res.u(0) == doctest::Approx(0.5f).epsilon(1e-5));
    CHECK(res.u(1) == doctest::Approx(0.5f).epsilon(1e-5));

    const auto u_op = filt(u_des);
    CHECK(u_op(0) == doctest::Approx(res.u(0)));
    CHECK(u_op(1) == doctest::Approx(res.u(1)));
}

TEST_CASE("ActionGovernor tick-varying b") {
    // Fixed normal A = [1, 0]; b changes: first u0 ≤ 2 (safe), then u0 ≤ 0.5
    ActionGovernor<2, 1, double> filt{Matrix<1, 2>{{1.0, 0.0}}, ColVec<1>{2.0}};
    const ColVec<2>              u_des{1.0, 0.0};

    auto res = filt.filter(u_des);
    REQUIRE(res.success);
    CHECK_FALSE(res.modified);

    res = filt.filter(u_des, ColVec<1>{0.5});
    REQUIRE(res.success);
    CHECK(res.modified);
    CHECK(res.u(0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(res.u(1) == doctest::Approx(0.0).epsilon(1e-9));
}

// =============================================================================
// CBF row helper
// =============================================================================

TEST_CASE("cbf_relative_degree_1 builds the affine row") {
    // ḣ = Lf + Lg·u ≥ −α h  →  (−Lg) u ≤ Lf + α h
    // Lg = (1, 2), Lf = 0.5, h = 1, α = 2 → A = (−1, −2), b = 0.5 + 2 = 2.5
    const auto cbf = design::cbf_relative_degree_1(ColVec<2>{1.0, 2.0}, 0.5, 1.0, 2.0);
    CHECK(cbf.A(0, 0) == doctest::Approx(-1.0));
    CHECK(cbf.A(0, 1) == doctest::Approx(-2.0));
    CHECK(cbf.b == doctest::Approx(2.5));
}

TEST_CASE("cbf_relative_degree_1 enforces barrier via project_affine") {
    // Scalar input: ẋ = u, h = 1 − x (safe for x ≤ 1). At x = 0.9, h = 0.1,
    // Lf = 0 (no drift), Lg = −1 (∂h/∂x · g = −1 · 1).
    // CBF: −α h ≤ ḣ = −u  ⇒  u ≤ α h. With α = 1, u ≤ 0.1.
    // Desired u_des = 1.0 must project to 0.1.
    constexpr double x = 0.9;
    constexpr double h = 1.0 - x;
    constexpr double Lf = 0.0;
    constexpr double Lg = -1.0; // ∂h/∂x * g
    constexpr double alpha = 1.0;

    const auto cbf = design::cbf_relative_degree_1(ColVec<1>{Lg}, Lf, h, alpha);
    const auto res = design::project_affine(ColVec<1>{1.0}, cbf.A, ColVec<1>{cbf.b});
    REQUIRE(res.success);
    CHECK(res.u(0) == doctest::Approx(alpha * h).epsilon(1e-9));

    // Already safe: u_des = 0.05 < 0.1
    const auto res_safe = design::project_affine(ColVec<1>{0.05}, cbf.A, ColVec<1>{cbf.b});
    REQUIRE(res_safe.success);
    CHECK_FALSE(res_safe.modified);
    CHECK(res_safe.u(0) == doctest::Approx(0.05).epsilon(1e-9));
}

TEST_CASE("cbf_relative_degree_1 is constexpr") {
    constexpr auto cbf = design::cbf_relative_degree_1(ColVec<1>{-1.0}, 0.0, 0.1, 1.0);
    static_assert(cbf.A(0, 0) == 1.0); // −Lg = −(−1) = 1
    static_assert(cbf.b - 0.1 < 1e-12 && cbf.b - 0.1 > -1e-12);

    constexpr auto res = design::project_affine(ColVec<1>{1.0}, cbf.A, ColVec<1>{cbf.b});
    static_assert(res.success);
    static_assert(res.u(0) - 0.1 < 1e-9 && res.u(0) - 0.1 > -1e-9);
}
