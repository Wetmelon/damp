// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/math/complex.hpp"
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

    auto ctrl = LQR{res}; // from the gain design

    SUBCASE("u = -K x") {
        const ColVec<2> x{{1.0, 2.0}};
        ctrl.update(x);
        const auto u = ctrl.control();
        CHECK(u[0] == doctest::Approx(-(res.K(0, 0) * 1.0 + res.K(0, 1) * 2.0)));
    }

    SUBCASE("tracking error is zero at the reference") {
        const ColVec<2> x{{3.0, -1.0}};
        ctrl.update(x);
        const auto u = ctrl.control(x, {}); // x_ref == x
        CHECK(u[0] == doctest::Approx(0.0));
    }
}

namespace {
constexpr auto ctad_res = design::discrete_lqr(
    Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}}, Matrix<2, 1>{{0.0}, {0.1}},
    Matrix<2, 2>{{1.0, 0.0}, {0.0, 1.0}}, Matrix<1, 1>{{1.0}}
);
constinit auto ctad_law = LQR{ctad_res.as<float>()}; // sizes and scalar deduced from the result

// The scalar comes from the result: a double design does not narrow to the float default.
static_assert(std::is_same_v<decltype(LQR{ctad_res}), LQR<2, 1, double>>);
static_assert(std::is_same_v<decltype(LQR{ctad_res.as<float>()}), LQR<2, 1, float>>);
static_assert(std::is_same_v<decltype(ctad_law), LQR<2, 1, float>>);
// From a gain matrix (explicit constructor, direct-initialization).
static_assert(std::is_same_v<decltype(LQR{ctad_res.K}), LQR<2, 1, double>>);
static_assert(std::is_same_v<decltype(LQR{ctad_res.K.as<float>()}), LQR<2, 1, float>>);
// Explicit arguments still convert.
static_assert(std::is_same_v<decltype(LQR<2, 1, float>{ctad_res}), LQR<2, 1, float>>);
} // namespace

TEST_CASE("LQR deduces sizes and scalar from its design") {
    CHECK(ctad_law.K(0, 0) == doctest::Approx(static_cast<float>(ctad_res.K(0, 0))));
    CHECK(ctad_law.K(0, 1) == doctest::Approx(static_cast<float>(ctad_res.K(0, 1))));
}

namespace {

// Discrete mass-spring-damper at Ts = 0.05 s, both states measured: no integrator, so
// u = −K·x settles off the reference without the target's feedforward.
constexpr StateSpace<2, 1, 2, double> spring_plant() {
    constexpr double Ts = 0.05;
    return StateSpace<2, 1, 2, double>{
        .A = Matrix<2, 2>{{1.0, Ts}, {-4.0 * Ts, 1.0 - (0.5 * Ts)}},
        .B = Matrix<2, 1>{{0.0}, {Ts}},
        .C = Matrix<2, 2>::identity(),
        .Ts = Ts,
    };
}

} // namespace

TEST_CASE("design::lqr with TrackedOutputs holds the tracked output at r") {
    static constexpr auto res = design::lqr(spring_plant(), {.Q = Matrix<2, 2>::identity(), .R = 0.1}, design::TrackedOutputs<1>{{0}});
    static_assert(res.success);
    static constinit LQR law{res.as<float>()};
    static_assert(std::is_same_v<decltype(law), LQR<2, 1, float, 2, 1>>);

    // The target: x* = Nₓ·r = [r, 0], u* = Nᵤ·r = 4r (the spring's force at x = r).
    CHECK(res.target.target_r(0, 0) == doctest::Approx(1.0));
    CHECK(res.target.target_r(1, 0) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(res.target.target_r(2, 0) == doctest::Approx(4.0));

    auto            sys = spring_plant();
    LQR             tracking{res};
    LQR             regulator{design::lqr(sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1})};
    ColVec<2>       xt{};
    ColVec<2>       xr{};
    const ColVec<1> r{0.3};
    ColVec<1>       u{};
    for (int k = 0; k < 400; ++k) {
        u = tracking.step(r, xt);
        xt = (sys.A * xt) + (sys.B * u);
        // The regulator with the same setpoint as a state, x_ref = [r, 0]: the u* the target adds is missing.
        xr = (sys.A * xr) + (sys.B * regulator.step(ColVec<2>{r[0], 0.0}, {}, xr));
    }
    CHECK(xt[0] == doctest::Approx(0.3).epsilon(1e-9));
    CHECK(u[0] == doctest::Approx(1.2).epsilon(1e-9));
    CHECK(std::abs(xr[0] - 0.3) > 0.01);
}

TEST_CASE("design::lqr with TrackedOutputs refuses an output it cannot hold") {
    // The velocity of a spring at rest is zero: no steady state holds it at r ≠ 0.
    const auto res = design::lqr(spring_plant(), {.Q = Matrix<2, 2>::identity(), .R = 0.1}, design::TrackedOutputs<1>{{1}});
    CHECK_FALSE(res.success);
    CHECK(res.status == design::DesignStatus::singular);
    // An output out of range.
    const auto bad = design::lqr(spring_plant(), {.Q = Matrix<2, 2>::identity(), .R = 0.1}, design::TrackedOutputs<1>{{2}});
    CHECK(bad.status == design::DesignStatus::bad_weights);
}

TEST_CASE("design::lqr with TrackedOutputs refuses an output no steady state holds, with more inputs than outputs") {
    // Two parallel input columns on a double integrator: the velocity has no steady state at r ≠ 0
    // with either input, so two redundant inputs are refused like one.
    const StateSpace<2, 2, 2> sys{
        .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
        .B = Matrix<2, 2>{{0.005, 0.01}, {0.1, 0.2}},
        .C = Matrix<2, 2>::identity(),
        .Ts = 0.1,
    };
    const auto weights = design::LQRConfig<2, 2>{.Q = Matrix<2, 2>::identity(), .R = Matrix<2, 2>::identity()};
    const auto res = design::lqr(sys, weights, design::TrackedOutputs<1>{{1}});
    CHECK_FALSE(res.success);
    CHECK(res.status == design::DesignStatus::singular);
    CHECK_FALSE(res.is_stable());
    // The position does have one (NZ < NU, the minimum-norm target), and the law holds it.
    const auto pos = design::lqr(sys, weights, design::TrackedOutputs<1>{{0}});
    REQUIRE(pos.success);
    auto      law = LQR{pos};
    ColVec<2> x{};
    for (int k = 0; k < 2000; ++k) {
        x = (sys.A * x) + (sys.B * law.step(ColVec<1>{1.0}, x));
    }
    CHECK(x[0] == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("LQRResult::to_ss closed around the plant has the design's poles") {
    // The plant reports its state, the law's input; unity positive feedback closes u = −Kx.
    const StateSpace<2, 1, 2> sys{
        .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
        .B = Matrix<2, 1>{{0.005}, {0.1}},
        .C = Matrix<2, 2>::identity(),
        .Ts = 0.1,
    };
    const auto res = design::lqr(sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1});
    REQUIRE(res.success);

    const auto ss = res.to_ss();
    CHECK(ss.D(0, 1) == doctest::Approx(-res.K(0, 1)));
    CHECK(ss.Ts == doctest::Approx(sys.Ts));

    StateSpace<0, 2, 2> positive{.Ts = sys.Ts};
    positive.D = -Matrix<2, 2>::identity();
    const auto loop = series(ss, sys);
    REQUIRE(loop.has_value());
    const auto cl = feedback(*loop, positive);
    REQUIRE(cl.has_value());

    const auto poles = mat::compute_eigenvalues(cl->A).values;
    for (size_t i = 0; i < 2; ++i) {
        double nearest = std::numeric_limits<double>::max();
        for (size_t j = 0; j < 2; ++j) {
            nearest = damp::min(nearest, (poles[j] - res.e[i]).abs());
        }
        CHECK(nearest < 1e-9);
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
    )
                            .value();

    CHECK(cost.Q(0, 0) == doctest::Approx(h));
    CHECK(cost.N(0, 0) == doctest::Approx(h * h / 2.0));
    CHECK(cost.R(0, 0) == doctest::Approx(h + h * h * h / 3.0));
}

TEST_CASE_TEMPLATE("discretize_lqr_cost and lqr(sys, Ts, w) keep a fast stable mode", T, float, double) {
    //! A = [[−1, 1], [0, a]] sampled at Ts = 1. One exponential of the Van Loan block holds e^{|a|·Ts}
    //! and its product cancels: in double a = −40 gave Nd₂ = −0.0154 and a wrong-sign gain, a = −200
    //! Rd = −6.6e63; in float a = −20 gave Nd₂ 19× too large and a = −200 NaN, all as values.
    //! Reference: SciPy quad_vec of ∫₀^Ts e^{Āᵀt} Q̄ e^{Āt} dt, then solve_discrete_are with s = Nd.
    struct Case {
        double a, Rd, Nd0, Nd1, K0, K1;
    };
    const Case cases[] = {
        {-20.0, 1.0026832864045, 0.0089770065728728, 0.00171934119308932, 0.0149975986964182, 0.00203212486167899},
        {-40.0, 1.00070039834416, 0.00473320365573618, 0.00043366876165428, 0.00770189356221218, 0.00050952719860189},
        {-200.0, 1.00002896485745, 0.000988203383289507, 1.7465533639371e-05, 0.00157238224993514, 2.04006799641705e-05},
    };
    const double tol = std::is_same_v<T, float> ? 1e-4 : 1e-9;
    for (const Case& c : cases) {
        CAPTURE(c.a);
        const StateSpace<2, 1, 1, T>        sys{.A = Matrix<2, 2, T>{{T{-1}, T{1}}, {T{0}, static_cast<T>(c.a)}}, .B = Matrix<2, 1, T>{{T{0}}, {T{1}}}, .C = Matrix<1, 2, T>{{T{1}, T{0}}}};
        const design::LQRConfig<2, 1, 0, T> w{.Q = Matrix<2, 2, T>::identity(), .R = Matrix<1, 1, T>{{T{1}}}};
        const auto                          cost = design::discretize_lqr_cost(sys.A, sys.B, w.Q, w.R, T{1});
        REQUIRE(cost.has_value());
        CHECK(static_cast<double>(cost->R(0, 0)) == doctest::Approx(c.Rd).epsilon(tol));
        CHECK(static_cast<double>(cost->N(0, 0)) == doctest::Approx(c.Nd0).epsilon(tol));
        CHECK(static_cast<double>(cost->N(1, 0)) == doctest::Approx(c.Nd1).epsilon(10.0 * tol));
        const auto res = design::lqr(sys, T{1}, w);
        REQUIRE(res.success);
        CHECK(static_cast<double>(res.K(0, 0)) == doctest::Approx(c.K0).epsilon(tol));
        CHECK(static_cast<double>(res.K(0, 1)) == doctest::Approx(c.K1).epsilon(10.0 * tol));
    }
}

TEST_CASE("design::lqr on a continuous plant with Ts discretizes the cost (differs from naive)") {
    const Matrix<2, 2> A{{0.0, 1.0}, {0.0, 0.0}}; // double integrator
    const Matrix<2, 1> B{{0.0}, {1.0}};
    const Matrix<2, 2> Q{{1.0, 0.0}, {0.0, 1.0}};
    const Matrix<1, 1> R{{1.0}};
    const double       Ts = 0.2;

    const auto proper = design::lqr(StateSpace{.A = A, .B = B, .C = Matrix<1, 2>{{1.0, 0.0}}}, Ts, {.Q = Q, .R = R});
    REQUIRE(proper.success);
    CHECK(proper.Ts == Ts);
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
    CHECK(res.Ts == 0.0);
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

    const auto via_matlab = matlab::lqr(sys.A, sys.B, Q, R);
    const auto via_design = design::continuous_lqr(sys.A, sys.B, Q, R);
    REQUIRE(via_matlab.success);
    CHECK(via_matlab.status == via_design.status);
    for (size_t j = 0; j < 4; ++j) {
        CHECK(via_matlab.K(0, j) == doctest::Approx(via_design.K(0, j)));
    }
}

TEST_CASE("matlab::dlqr and matlab::lqrd forward to the core names") {
    const Matrix<1, 1> A{{1.0}};
    const Matrix<1, 1> B{{1.0}};
    const Matrix<1, 1> Q{{1.0}};
    const Matrix<1, 1> R{{1.0}};

    const auto via_long = design::discrete_lqr(A, B, Q, R);
    const auto via_dlqr = matlab::dlqr(A, B, Q, R);
    REQUIRE(via_long.success);
    REQUIRE(via_dlqr.success);
    CHECK(via_dlqr.K(0, 0) == doctest::Approx(via_long.K(0, 0)));
    CHECK(via_dlqr.S(0, 0) == doctest::Approx(via_long.S(0, 0)));

    const Matrix<2, 2> Ac{{0.0, 1.0}, {0.0, 0.0}};
    const Matrix<2, 1> Bc{{0.0}, {1.0}};
    const Matrix<2, 2> Qc{{1.0, 0.0}, {0.0, 1.0}};
    const Matrix<1, 1> Rc{{1.0}};
    const double       Ts = 0.2;

    StateSpace<2, 1, 2> sys{.A = Ac, .B = Bc, .C = Matrix<2, 2>::identity()};
    const auto          from_cont = design::lqr(sys, Ts, {.Q = Qc, .R = Rc});
    const auto          via_lqrd = matlab::lqrd(Ac, Bc, Qc, Rc, Ts);
    REQUIRE(from_cont.success);
    REQUIRE(via_lqrd.success);
    CHECK(via_lqrd.Ts == Ts);
    CHECK(via_lqrd.K(0, 0) == doctest::Approx(from_cont.K(0, 0)));
    CHECK(via_lqrd.K(0, 1) == doctest::Approx(from_cont.K(0, 1)));

    const auto via_lqrd_sys = matlab::lqrd(sys, Qc, Rc, Ts);
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

TEST_CASE("design::lqr sampled-data design converges with expected gain scale") {
    constexpr Matrix<2, 2> A_c{{-1.0, 1.0}, {0.0, -2.0}};
    constexpr Matrix<2, 1> B_c{{1.0}, {0.5}};
    constexpr Matrix<2, 2> Q{{1.0, 0.0}, {0.0, 1.0}};
    constexpr Matrix<1, 1> R{{1.0}};
    constexpr double       Ts = 0.1;

    constexpr auto result = design::lqr(StateSpace{.A = A_c, .B = B_c, .C = Matrix<1, 2>{{1.0, 0.0}}}, Ts, {.Q = Q, .R = R});
    static_assert(result.success);
    CHECK(doctest::Approx(result.K(0, 0)).epsilon(1e-3) == 0.414);
    CHECK(doctest::Approx(result.K(0, 1)).epsilon(1e-3) == 0.231);
    CHECK(result.S(0, 0) > 0.0);
    CHECK(result.S(1, 1) > 0.0);
}

// ---------------------------------------------------------------------------
// design::lqr: one name, dispatching on the plant's time domain
// ---------------------------------------------------------------------------

namespace {
constexpr Matrix<2, 2> kDiA{{0.0, 1.0}, {0.0, 0.0}}; // double integrator
constexpr Matrix<2, 1> kDiB{{0.0}, {1.0}};
constexpr Matrix<1, 2> kDiC{{1.0, 0.0}};
constexpr StateSpace   kDiPlant{.A = kDiA, .B = kDiB, .C = kDiC};
constexpr auto         kDiDiscrete = *discretize(kDiPlant, 0.1);

constexpr design::LQRConfig<2, 1> kDiW{.Q = Matrix<2, 2>::identity(), .R = 1.0};

// Row 2 at compile time: continuous plant + Ts is the sampled-data design, and the law deploys from it.
constexpr auto kSampled = design::lqr(kDiPlant, 0.2, {.Q = kDiW.Q, .R = kDiW.R});
static_assert(kSampled.success && kSampled.status == design::DesignStatus::ok && kSampled.Ts == 0.2);
constinit auto kSampledLaw = LQR{kSampled.as<float>()};
static_assert(std::is_same_v<decltype(kSampledLaw), LQR<2, 1, float>>);

// Row 4 at compile time: a discrete plant plus Ts is refused, even at the plant's own Ts.
constexpr auto kTwice = design::lqr(kDiDiscrete, 0.1, kDiW);
static_assert(!kTwice.success && kTwice.status == design::DesignStatus::invalid_sample_time);

// Row 1 designs fine but is a continuous gain (Ts = 0), not a sampled law: it does not deploy.
constexpr auto kContinuous = design::lqr(kDiPlant, kDiW);
static_assert(kContinuous.success && kContinuous.Ts == 0.0);
// constinit auto kContinuousLaw = LQR{kContinuous.as<float>()}; // ill-formed: refuse_failed_design under constant evaluation
} // namespace

TEST_CASE("LQR from a continuous result: refused when constant-evaluated, zero gain at run time") {
    auto law = LQR{kContinuous}; // a continuous gain on a discrete tick: loaded as a failed design
    CHECK(law.K(0, 0) == 0.0);
    CHECK(law.K(0, 1) == 0.0);
}

TEST_SUITE("design::lqr dispatch table") {
    TEST_CASE("sys.Ts = 0, no Ts: CARE, same as continuous_lqr") {
        const auto res = design::lqr(kDiPlant, kDiW);
        const auto ref = design::continuous_lqr(kDiA, kDiB, kDiW.Q, kDiW.R);
        REQUIRE(res.success);
        CHECK(res.status == design::DesignStatus::ok);
        CHECK(res.Ts == 0.0);
        CHECK(res.K(0, 0) == doctest::Approx(ref.K(0, 0)));
        CHECK(res.K(0, 1) == doctest::Approx(ref.K(0, 1)));
        CHECK(res.K(0, 0) == doctest::Approx(1.0)); // MATLAB® lqr: K = [1 √3]
        CHECK(res.K(0, 1) == doctest::Approx(std::sqrt(3.0)));
    }

    TEST_CASE("sys.Ts = 0, Ts > 0: ZOH + Van Loan cost, then DARE") {
        //! The numbers the removed discrete_lqr_from_continuous gave (scipy golden data).
        CHECK(kSampled.K(0, 0) == doctest::Approx(0.843695465061778));
        CHECK(kSampled.K(0, 1) == doctest::Approx(1.548939304133434));
        CHECK(kSampled.is_stable());
        CHECK(kSampledLaw.K(0, 0) == doctest::Approx(0.843695465061778f));
    }

    TEST_CASE("sys.Ts > 0, no Ts: DARE on the plant as given, same as discrete_lqr") {
        const auto res = design::lqr(kDiDiscrete, kDiW);
        const auto ref = design::discrete_lqr(kDiDiscrete.A, kDiDiscrete.B, kDiW.Q, kDiW.R);
        REQUIRE(res.success);
        CHECK(res.Ts == 0.1);
        CHECK(res.K(0, 0) == doctest::Approx(ref.K(0, 0)));
        CHECK(res.K(0, 1) == doctest::Approx(ref.K(0, 1)));
        CHECK(res.is_stable());
    }

    TEST_CASE("sys.Ts > 0 and a Ts: refused, even when equal") {
        for (const double Ts : {0.1, 0.2}) {
            const auto res = design::lqr(kDiDiscrete, Ts, kDiW);
            CHECK_FALSE(res.success);
            CHECK(res.status == design::DesignStatus::invalid_sample_time);
        }
    }

    TEST_CASE("Ts <= 0 or non-finite: refused") {
        const double inf = std::numeric_limits<double>::infinity();
        for (const double Ts : {0.0, -0.1, inf, -inf}) {
            const auto res = design::lqr(kDiPlant, Ts, kDiW);
            CHECK_FALSE(res.success);
            CHECK(res.status == design::DesignStatus::invalid_sample_time);
        }
        auto bad_plant = kDiPlant;
        bad_plant.Ts = -0.1;
        CHECK(design::lqr(bad_plant, kDiW).status == design::DesignStatus::invalid_sample_time);
        bad_plant.Ts = inf;
        CHECK(design::lqr(bad_plant, kDiW).status == design::DesignStatus::invalid_sample_time);
    }

    TEST_CASE("the sampled design takes discretize's status: a non-finite entry, an overflow") {
        auto bad_plant = kDiPlant;
        bad_plant.A(0, 1) = std::numeric_limits<double>::quiet_NaN();
        CHECK(design::lqr(bad_plant, 0.1, kDiW).status == design::DesignStatus::bad_weights);
        bad_plant.A = Matrix<2, 2>{{1000.0, 1.0}, {0.0, 1000.0}}; // e^(1000·Ts) past the range of double
        CHECK(design::lqr(bad_plant, 10.0, kDiW).status == design::DesignStatus::not_converged);
        CHECK(design::lqr(bad_plant, 1e-3, kDiW).success); // the same plant at a rate it fits
    }

    TEST_CASE("an unstabilizable (A, B) DARE says not_stabilizable") {
        //! The mode at 1.5 is unstable and B cannot reach it.
        const Matrix<2, 2> A{{1.5, 0.0}, {0.0, 0.5}};
        const Matrix<2, 1> B{{0.0}, {1.0}};
        const auto         res = design::discrete_lqr(A, B, Matrix<2, 2>::identity(), Matrix<1, 1>{{1.0}});
        CHECK_FALSE(res.success);
        CHECK(res.status == design::DesignStatus::not_stabilizable);
    }

    TEST_CASE("R not positive definite in CARE says bad_weights") {
        for (const double r : {0.0, -1.0}) {
            const auto res = design::continuous_lqr(kDiA, kDiB, Matrix<2, 2>::identity(), Matrix<1, 1>{{r}});
            CHECK_FALSE(res.success);
            CHECK(res.status == design::DesignStatus::bad_weights);
        }
    }

    TEST_CASE("Q not PSD says bad_weights, CARE and DARE") {
        const Matrix<2, 2> Q{{1.0, 0.0}, {0.0, -1.0}};
        CHECK(design::lqr(kDiPlant, {.Q = Q, .R = kDiW.R}).status == design::DesignStatus::bad_weights);
        CHECK(design::lqr(kDiDiscrete, {.Q = Q, .R = kDiW.R}).status == design::DesignStatus::bad_weights);
    }

    TEST_CASE("success and status agree, and survive .as<float>()") {
        const auto ok = design::lqr(kDiPlant, 0.2, kDiW);
        const auto bad = design::lqr(kDiDiscrete, 0.2, kDiW);
        CHECK(ok.success == (ok.status == design::DesignStatus::ok));
        CHECK(bad.success == (bad.status == design::DesignStatus::ok));
        CHECK(ok.as<float>().status == design::DesignStatus::ok);
        CHECK(ok.as<float>().Ts == 0.2f);
        CHECK(bad.as<float>().status == design::DesignStatus::invalid_sample_time);
        CHECK_FALSE(bad.as<float>().success);
    }

    TEST_CASE("design::lqr fails on a non-positive sample time") {
        const StateSpace<2, 1, 1>     plant{.A = {{0.0, 1.0}, {0.0, 0.0}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        const design::LQRConfig<2, 1> w{.Q = Matrix<2, 2>::identity(), .R = 1.0};
        CHECK_FALSE(design::lqr(plant, 0.0, w).success);
        CHECK(design::lqr(plant, 0.01, w).success);
    }
}
