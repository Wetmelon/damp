// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// Classic CTMS "Inverted Pendulum: State-Space Methods" plant (continuous-time).
// https://ctms.engin.umich.edu/CTMS/  (Control: State-Space section)
//   M=0.5 m=0.2 b=0.1 I=0.006 g=9.8 l=0.3,  p = I(M+m)+Mml^2 = 0.0132
struct CtmsPendulum {
    Matrix<4, 4> A{
        {0.0, 1.0, 0.0, 0.0},
        {0.0, -0.1818181818, 2.6727272727, 0.0},
        {0.0, 0.0, 0.0, 1.0},
        {0.0, -0.4545454545, 31.1818181818, 0.0},
    };
    ColVec<4> B{0.0, 1.8181818182, 0.0, 4.5454545455};
};

} // namespace

TEST_CASE("discrete_lqr scalar gain matches hand solution") {
    //! A=1,B=1,Q=1,R=1: DARE gives S²−S−1=0 → S=φ=1.618…, K=S/(1+S)=0.618…
    const auto res = design::discrete_lqr(
        Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}}
    );

    REQUIRE(res.success);
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    CHECK(res.S(0, 0) == doctest::Approx(phi));
    CHECK(res.K(0, 0) == doctest::Approx(phi / (1.0 + phi))); // 0.6180…
    CHECK(res.is_stable());                                   // pole A−BK = 0.382
    CHECK(res.e[0].abs() == doctest::Approx(1.0 - phi / (1.0 + phi)));
}

TEST_CASE("LQR runtime regulation and tracking laws") {
    const design::LQRResult<2, 1> res = design::discrete_lqr(
        Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}}, Matrix<2, 1>{{0.0}, {0.1}},
        Matrix<2, 2>{{1.0, 0.0}, {0.0, 1.0}}, Matrix<1, 1>{{1.0}}
    );
    REQUIRE(res.success);

    const LQR<2, 1, double> ctrl{res}; // implicit from result

    SUBCASE("u = -K x") {
        const ColVec<2> x{{1.0, 2.0}};
        const auto      u = ctrl.control(x);
        CHECK(u[0] == doctest::Approx(-(res.K(0, 0) * 1.0 + res.K(0, 1) * 2.0)));
    }

    SUBCASE("tracking error is zero at the reference") {
        const ColVec<2> x{{3.0, -1.0}};
        const auto      u = ctrl.control(x, x); // x_ref == x (reference-first order)
        CHECK(u[0] == doctest::Approx(0.0));
    }
}

TEST_CASE("LQRResult::as<float>() round-trips every field") {
    const auto res_d = design::discrete_lqr(
        Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}}
    );
    const auto res_f = res_d.as<float>();

    CHECK(res_f.success);
    CHECK(res_f.K(0, 0) == doctest::Approx(static_cast<float>(res_d.K(0, 0))));
    CHECK(res_f.S(0, 0) == doctest::Approx(static_cast<float>(res_d.S(0, 0))));
    CHECK(res_f.e[0].abs() == doctest::Approx(static_cast<float>(res_d.e[0].abs())));
}

TEST_CASE("discretize_lqr_cost matches closed-form for A=0 integrator") {
    //! ẋ = u (A=0,B=1), x(τ)=x0+uτ over a sample h gives exact weights:
    //!   Qd = h,  Nd = h²/2,  Rd = h + h³/3
    const double h = 0.1;
    const auto   cost = design::discretize_lqr_cost(
        Matrix<1, 1>{{0.0}}, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}}, h
    );

    CHECK(cost.Q(0, 0) == doctest::Approx(h));
    CHECK(cost.N(0, 0) == doctest::Approx(h * h / 2.0));
    CHECK(cost.R(0, 0) == doctest::Approx(h + h * h * h / 3.0));
}

TEST_CASE("discrete_lqr_from_continuous discretizes cost (differs from naive)") {
    const Matrix<2, 2> A{{0.0, 1.0}, {0.0, 0.0}}; // double integrator
    const Matrix<2, 1> B{{0.0}, {1.0}};
    const Matrix<2, 2> Q{{1.0, 0.0}, {0.0, 1.0}};
    const Matrix<1, 1> R{{1.0}};
    const double       Ts = 0.2;

    const auto proper = design::discrete_lqr_from_continuous(A, B, Q, R, Ts);
    REQUIRE(proper.success);
    CHECK(proper.is_stable());

    //! Golden gain from an independent scipy implementation (Van Loan cost
    //! discretization + dlqr), see AGENTS.md "golden reference data".
    CHECK(proper.K(0, 0) == doctest::Approx(0.843695465061778));
    CHECK(proper.K(0, 1) == doctest::Approx(1.548939304133434));

    //! And it must differ from the naive path (ZOH dynamics, continuous cost) —
    //! otherwise the cost discretization is a silent no-op.
    StateSpace<2, 1, 2, double, 2, 2> sys_c{A, B, Matrix<2, 2>::identity()};
    const auto                        sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);
    const auto                        naive = design::discrete_lqr(sys_d.A, sys_d.B, Q, R);
    REQUIRE(naive.success);
    CHECK(std::abs(proper.K(0, 0) - naive.K(0, 0)) > 1e-4);
}

TEST_CASE("continuous_lqr reproduces the CTMS inverted-pendulum gain") {
    const CtmsPendulum sys;

    // CTMS weighting: Q = C'C with the cart-position and pendulum-angle weights
    // bumped to 5000 and 100; R = 1.
    const Matrix<4, 4> Q{
        {5000.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 100.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
    };
    const Matrix<1, 1> R{{1.0}};

    const auto res = design::continuous_lqr(sys.A, sys.B, Q, R);
    REQUIRE(res.success);

    // MATLAB®'s K = lqr(A,B,Q,R) for this problem: [-70.7107 -37.8345 105.5298 20.9238].
    CHECK(res.K(0, 0) == doctest::Approx(-70.7107).epsilon(1e-3));
    CHECK(res.K(0, 1) == doctest::Approx(-37.8345).epsilon(1e-3));
    CHECK(res.K(0, 2) == doctest::Approx(105.5298).epsilon(1e-3));
    CHECK(res.K(0, 3) == doctest::Approx(20.9238).epsilon(1e-3));

    // Continuous-time stability: every closed-loop pole has negative real part.
    CHECK_FALSE(res.discrete);
    CHECK(res.is_stable());
    for (size_t i = 0; i < 4; ++i) {
        CHECK(res.e[i].real() < 0.0);
    }
}

TEST_CASE("matlab::lqr alias returns the same gain as continuous_lqr") {
    const CtmsPendulum sys;
    const Matrix<4, 4> Q{
        {5000.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 100.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
    };
    const Matrix<1, 1> R{{1.0}};

    const auto K = matlab::lqr(sys.A, sys.B, Q, R);
    const auto K2 = design::continuous_lqr(sys.A, sys.B, Q, R).K;
    for (size_t j = 0; j < 4; ++j) {
        CHECK(K(0, j) == doctest::Approx(K2(0, j)));
    }
}

TEST_CASE("design::dlqr and design::lqrd match descriptive names") {
    const Matrix<1, 1> A{{1.0}};
    const Matrix<1, 1> B{{1.0}};
    const Matrix<1, 1> Q{{1.0}};
    const Matrix<1, 1> R{{1.0}};

    const auto via_long = design::discrete_lqr(A, B, Q, R);
    const auto via_dlqr = design::dlqr(A, B, Q, R);
    REQUIRE(via_long.success);
    REQUIRE(via_dlqr.success);
    CHECK(via_dlqr.K(0, 0) == doctest::Approx(via_long.K(0, 0)));
    CHECK(via_dlqr.S(0, 0) == doctest::Approx(via_long.S(0, 0)));

    const Matrix<2, 2> Ac{{0.0, 1.0}, {0.0, 0.0}};
    const Matrix<2, 1> Bc{{0.0}, {1.0}};
    const Matrix<2, 2> Qc{{1.0, 0.0}, {0.0, 1.0}};
    const Matrix<1, 1> Rc{{1.0}};
    const double       Ts = 0.2;

    const auto from_cont = design::discrete_lqr_from_continuous(Ac, Bc, Qc, Rc, Ts);
    const auto via_lqrd = design::lqrd(Ac, Bc, Qc, Rc, Ts);
    REQUIRE(from_cont.success);
    REQUIRE(via_lqrd.success);
    CHECK(via_lqrd.K(0, 0) == doctest::Approx(from_cont.K(0, 0)));
    CHECK(via_lqrd.K(0, 1) == doctest::Approx(from_cont.K(0, 1)));

    StateSpace<2, 1, 2> sys{.A = Ac, .B = Bc, .C = Matrix<2, 2>::identity()};
    const auto          via_lqrd_sys = design::lqrd(sys, Qc, Rc, Ts);
    REQUIRE(via_lqrd_sys.success);
    CHECK(via_lqrd_sys.K(0, 0) == doctest::Approx(from_cont.K(0, 0)));
    CHECK(via_lqrd_sys.K(0, 1) == doctest::Approx(from_cont.K(0, 1)));
}

TEST_SUITE("Design: Compile-Time LQR with Cross-Term N") {
    TEST_CASE("design::discrete_lqr with cross-term N at compile time") {
        // Use constexpr to ensure compile-time evaluation
        constexpr auto result_no_N = design::discrete_lqr(
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}}
        );

        constexpr auto result_with_N = design::discrete_lqr(
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{0.5}}
        );

        static_assert(result_no_N.success);
        static_assert(result_with_N.success);

        // Gains should differ
        static_assert(result_no_N.K(0, 0) != result_with_N.K(0, 0));

        // Runtime checks
        CHECK(result_no_N.success);
        CHECK(result_with_N.success);
        CHECK(result_no_N.K(0, 0) != result_with_N.K(0, 0));
    }
}

TEST_CASE("LQRResult::as compile-time float round-trip") {
    constexpr auto lqr_d = design::discrete_lqr(
        Matrix<1, 1>{{1.0}},
        Matrix<1, 1>{{1.0}},
        Matrix<1, 1>{{1.0}},
        Matrix<1, 1>{{1.0}}
    );
    constexpr auto lqr_f = lqr_d.as<float>();
    static_assert(lqr_f.success);
    static_assert(lqr_f.K(0, 0) != 0.0f);
    CHECK(doctest::Approx(static_cast<double>(lqr_f.K(0, 0))).epsilon(1e-6) == lqr_d.K(0, 0));
}

TEST_CASE("design::discrete_lqr matches scipy/control golden data") {
    constexpr Matrix<2, 2> A{{0.9, 0.1}, {0.0, 0.8}};
    constexpr Matrix<2, 1> B{{0.1}, {0.2}};
    constexpr Matrix<2, 2> Q{{1.0, 0.0}, {0.0, 1.0}};
    constexpr Matrix<1, 1> R{{1.0}};

    constexpr auto result = design::discrete_lqr(A, B, Q, R);
    static_assert(result.success);

    CHECK(doctest::Approx(result.K(0, 0)).epsilon(1e-6) == 0.4156588386070948);
    CHECK(doctest::Approx(result.K(0, 1)).epsilon(1e-6) == 0.435575939822669);
    CHECK(doctest::Approx(result.S(0, 0)).epsilon(1e-6) == 4.201439333611463);
    CHECK(doctest::Approx(result.S(0, 1)).epsilon(1e-6) == 0.5954889098912259);
    CHECK(doctest::Approx(result.S(1, 0)).epsilon(1e-6) == 0.5954889098912259);
    CHECK(doctest::Approx(result.S(1, 1)).epsilon(1e-6) == 2.543807455993601);
    CHECK(doctest::Approx(result.e[0].real()).epsilon(1e-6) == 0.8102357330184564);
    CHECK(doctest::Approx(result.e[1].real()).epsilon(1e-6) == 0.7610831951563004);
}

TEST_CASE("design::discrete_lqr_from_continuous converges with expected gain scale") {
    constexpr Matrix<2, 2> A_c{{-1.0, 1.0}, {0.0, -2.0}};
    constexpr Matrix<2, 1> B_c{{1.0}, {0.5}};
    constexpr Matrix<2, 2> Q{{1.0, 0.0}, {0.0, 1.0}};
    constexpr Matrix<1, 1> R{{1.0}};
    constexpr double       Ts = 0.1;

    constexpr auto result = design::discrete_lqr_from_continuous(A_c, B_c, Q, R, Ts);
    static_assert(result.success);
    CHECK(doctest::Approx(result.K(0, 0)).epsilon(1e-3) == 0.414);
    CHECK(doctest::Approx(result.K(0, 1)).epsilon(1e-3) == 0.231);
    CHECK(result.S(0, 0) > 0.0);
    CHECK(result.S(1, 1) > 0.0);
}
