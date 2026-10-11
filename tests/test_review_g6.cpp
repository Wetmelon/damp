// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Regression tests for the g6 review (simulation + toolbox + matlab aliases).

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numbers>
#include <plotlypp/figure.hpp>
#include <vector>

#include "damp/matlab.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/simulation/cached_zoh.hpp"
#include "damp/simulation/hybrid.hpp"
#include "damp/simulation/integrator.hpp"
#include "damp/simulation/multirate.hpp"
#include "damp/simulation/npy_export.hpp"
#include "damp/simulation/plot_plotly.hpp"
#include "damp/simulation/simulate.hpp"
#include "damp/simulation/solver.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/conditioning.hpp"
#include "damp/toolbox/encoder.hpp"
#include "damp/toolbox/iec61131.hpp"
#include "damp/toolbox/logic.hpp"
#include "damp/toolbox/thermistor.hpp"
#include "damp/toolbox/timing.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::sim;

TEST_SUITE("Review g6: timers and logic") {
    TEST_CASE("TON<double> fires on the tick where n·dt reaches PT despite rounding") {
        // 0.1 summed 8 times is 0.7999999999999999 < 0.8; the timer must still fire on tick 8.
        plc::TON<double> ton;
        ton.PT = 0.8;
        for (int i = 0; i < 7; ++i) {
            CHECK_FALSE(ton(true, 0.1));
        }
        CHECK(ton(true, 0.1));
        CHECK(ton.ET == 0.8);
    }

    TEST_CASE("OnDelayTimer<double> fires after exactly ten 0.1 s ticks") {
        OnDelayTimer<double> t{1.0};
        for (int i = 0; i < 9; ++i) {
            CHECK_FALSE(t(true, 0.1));
        }
        CHECK(t(true, 0.1));
    }

    TEST_CASE("TP re-arms when IN rises right after a pulse that ended with IN low") {
        plc::TP<double> tp;
        tp.PT = 0.25;
        CHECK(tp(true, 0.125));           // start, ET = 0.125
        CHECK(tp(false, 0.125) == false); // ET = 0.25 -> pulse ends, IN already low
        CHECK(tp(true, 0.125));           // new rising edge must start a new pulse
    }

    TEST_CASE("TOF elapsed time stays 0 at cold start while IN was never high") {
        plc::TOF<double> tof;
        tof.PT = 0.5;
        for (int i = 0; i < 10; ++i) {
            CHECK_FALSE(tof(false, 0.1));
        }
        CHECK(tof.ET == 0.0);
    }

    TEST_CASE("CTD reset leaves Q consistent with CV = 0") {
        plc::CTD<uint32_t> ctd;
        ctd.PV = 3;
        ctd(false, true); // load
        ctd.reset();
        CHECK(ctd.CV == 0);
        CHECK(ctd.Q);
    }

    TEST_CASE("BLINK keeps its average period when dt does not divide the phase times") {
        plc::BLINK<double> blink{0.25, 0.25}; // 2 Hz
        RisingEdge         edge;
        int                n_edges = 0;
        for (int i = 0; i < 1000; ++i) { // 100 s at dt = 0.1
            if (edge(blink(true, 0.1))) {
                ++n_edges;
            }
        }
        CHECK(n_edges >= 199);
        CHECK(n_edges <= 201);
    }

    TEST_CASE("Counter saturates instead of overflowing") {
        Counter<int> c{std::numeric_limits<int>::max()};
        CHECK(c(true, false) == std::numeric_limits<int>::max());
        Counter<int> d{std::numeric_limits<int>::min()};
        CHECK(d(false, true) == std::numeric_limits<int>::min());
    }

    TEST_CASE("Timeout<uint32_t> expires when dt would wrap the accumulator") {
        Timeout<uint32_t> to{4'000'000'000u};
        CHECK_FALSE(to.tick(3'000'000'000u));
        CHECK(to.tick(2'000'000'000u));
    }
}

TEST_SUITE("Review g6: toolbox") {
    TEST_CASE("QuadratureDecoder X1 does not drift while dithering on an A edge") {
        QuadratureDecoder dec{QuadMode::X1};
        dec.update(false, true); // state 01
        for (int i = 0; i < 10; ++i) {
            dec.update(true, true);  // 01 -> 11 (A rises, forward)
            dec.update(false, true); // 11 -> 01 (A falls, reverse)
        }
        CHECK(dec.position() == 0);
    }

    TEST_CASE("soft_sign does not collapse to 0 when x² overflows") {
        CHECK(soft_sign(1e20f, 1.0f) == doctest::Approx(1.0f));
        CHECK(soft_sign(-1e200, 1.0) == doctest::Approx(-1.0));
        CHECK(soft_sign(0.0, 0.0) == 0.0);
        CHECK(soft_sign(1e-30, 0.0) == doctest::Approx(1.0)); // eps = 0 is the hard sign
    }

    TEST_CASE("steinhart_hart reports failure for a degenerate fit") {
        const auto bad = design::steinhart_hart<double>({10000.0, 273.15}, {10000.0, 298.15}, {3603.0, 323.15});
        CHECK_FALSE(bad.success);
        const auto good = design::steinhart_hart<double>({32650.0, 273.15}, {10000.0, 298.15}, {3603.0, 323.15});
        CHECK(good.success);
        CHECK(design::beta(10000.0, 298.15, 3950.0).success);
        CHECK_FALSE(design::beta(10000.0, 298.15, 0.0).success);
    }
}

TEST_SUITE("Review g6: integrators") {
    TEST_CASE("BackwardEuler nonlinear step is invariant to the state's units") {
        // ẋ = −λx, hλ = 0.5: implicit Euler gives x/(1 + hλ) whatever the scale of x.
        const auto       f = [](double, const ColVec<1>& x) { return ColVec<1>{-x(0)}; };
        BackwardEuler<1> be;
        for (double scale : {1.0, 1e-12, 1e12}) {
            const auto r = be.evolve(f, ColVec<1>{scale}, 0.0, 0.5);
            CHECK(r.x(0) / scale == doctest::Approx(1.0 / 1.5).epsilon(1e-9));
        }
    }

    TEST_CASE("Trapezoidal nonlinear overload agrees with the LTI overload") {
        const auto         f = [](double, const ColVec<1>& x) { return ColVec<1>{-x(0)}; };
        Trapezoidal<1>     tr;
        const auto         nl = tr.evolve(f, ColVec<1>{1.0}, 0.0, 0.5);
        const Matrix<1, 1> A{{-1.0}};
        const Matrix<1, 1> B{{0.0}};
        const auto         lti = tr.evolve(A, B, ColVec<1>{1.0}, ColVec<1>{0.0}, 0.5);
        CHECK(lti.x(0) == doctest::Approx(0.75 / 1.25));
        CHECK(nl.x(0) == doctest::Approx(lti.x(0)).epsilon(1e-12));
    }
}

TEST_SUITE("Review g6: solvers") {
    const auto decay = [](double, const ColVec<1>& x) { return ColVec<1>{-x(0)}; };

    TEST_CASE("FixedStepSolver reaches tf on a femtosecond time scale") {
        const auto r = fixed_solve<RK4>(decay, ColVec<1>{1.0}, {0.0, 1e-14}, 1e-15);
        CHECK(r.size() == 11);
        CHECK(r.t.back() == 1e-14);
    }

    TEST_CASE("FixedStepSolver takes no sliver step from accumulated time drift") {
        const auto r = fixed_solve<RK4>(decay, ColVec<1>{1.0}, {0.0, 100.0}, 0.1);
        CHECK(r.size() == 1001);
        CHECK(r.t.back() == 100.0);
        CHECK(r.t[3] == 0.1 * 3);
    }

    TEST_CASE("FixedStepSolver rejects a non-positive step instead of hanging") {
        const auto r = fixed_solve<RK4>(decay, ColVec<1>{1.0}, {0.0, 1.0}, 0.0);
        CHECK_FALSE(r.success);
    }

    TEST_CASE("Adaptive initial step uses Hairer/SciPy exponent 1/error_order") {
        const auto r = adaptive_solve<DP45>(decay, ColVec<1>{1.0}, {0.0, 1.0});
        // SciPy select_initial_step for y' = −y, y0 = 1, atol = 1e-6, rtol = 1e-3, order 4.
        const double scale = 1e-6 + 1e-3;
        const double d1 = 1.0 / scale;
        const double h0 = 0.01;
        const double d2 = std::abs(-0.99 + 1.0) / scale / h0;
        const double h1 = std::pow(0.01 / std::max(d1, d2), 1.0 / 5.0);
        REQUIRE(r.size() > 1);
        CHECK(r.t[1] == doctest::Approx(std::min(100 * h0, h1)).epsilon(1e-9));
    }

    TEST_CASE("Adaptive solver shrinks h after a non-finite error estimate") {
        // f is undefined (NaN) beyond x = 3: a too-large trial step must shrink, not grow.
        const auto f = [](double, const ColVec<1>& x) {
            return ColVec<1>{x(0) >= 3.0 ? std::numeric_limits<double>::quiet_NaN() : x(0) * x(0)};
        };
        const auto r = adaptive_solve<DP45>(f, ColVec<1>{1.0}, {0.0, 0.6}, {.first_step = 0.6});
        CHECK(r.success);
        CHECK(r.x.back()(0) == doctest::Approx(1.0 / (1.0 - 0.6)).epsilon(1e-3));
    }

    TEST_CASE("Zero crossing is located on the integrated trajectory, not a chord") {
        // Free fall from 1 m: impact at sqrt(2/g). DP45 is exact on the quadratic, so h grows large.
        const auto         fall = [](double, const ColVec<2>& x) { return ColVec<2>{x(1), -9.81}; };
        DP45<2>            dp;
        AdaptiveStepSolver solver(dp, AdaptiveOptions<2>{.first_step = 0.6});
        solver.add_zero_crossing([](double, const ColVec<2>& x) { return x(0); });
        solver.set_stop_condition([](double, const ColVec<2>& x) { return x(0) <= 0.0; });
        const auto r = solver.solve(fall, ColVec<2>{1.0, 0.0}, {0.0, 2.0});
        CHECK(r.t.back() == doctest::Approx(std::sqrt(2.0 / 9.81)).epsilon(1e-9));
        CHECK(std::abs(r.x.back()(0)) < 1e-9);
    }

    TEST_CASE("Earliest of several zero crossings in one step is the one recorded") {
        const auto         ramp = [](double, const ColVec<1>&) { return ColVec<1>{1.0}; };
        DP45<1>            dp;
        AdaptiveStepSolver solver(dp, AdaptiveOptions<1>{.first_step = 1.0});
        solver.add_zero_crossing([](double, const ColVec<1>& x) { return x(0) - 0.8; });
        solver.add_zero_crossing([](double, const ColVec<1>& x) { return x(0) - 0.3; });
        const auto r = solver.solve(ramp, ColVec<1>{0.0}, {0.0, 1.0});
        bool       hit = false;
        for (double t : r.t) {
            hit = hit || std::abs(t - 0.3) < 1e-9;
        }
        CHECK(hit);
    }

    TEST_CASE("Zero-crossing location is scale-invariant in z") {
        // z in tiny units (1e-15): the old absolute |z| < 1e-12 stop quit after one bisection.
        const auto         ramp = [](double, const ColVec<1>&) { return ColVec<1>{1.0}; };
        DP45<1>            dp;
        AdaptiveStepSolver solver(dp, AdaptiveOptions<1>{.first_step = 1.0});
        solver.add_zero_crossing([](double, const ColVec<1>& x) { return 1e-15 * (x(0) - 0.3); });
        const auto r = solver.solve(ramp, ColVec<1>{0.0}, {0.0, 1.0});
        bool       hit = false;
        for (double t : r.t) {
            hit = hit || std::abs(t - 0.3) < 1e-9;
        }
        CHECK(hit);
    }
}

namespace {
struct ConstController {
    double u{2.0};
    double control(double, double) { return u; }
    void   reset() {}
};
} // namespace

TEST_SUITE("Review g6: closed-loop simulation") {
    TEST_CASE("simulate_sampled ticks the controller once per period, no drift sliver") {
        int        calls = 0;
        const auto plant = [](double, const ColVec<1>& x, const ColVec<1>& u) { return ColVec<1>{-x(0) + u(0)}; };
        const auto output = [](const ColVec<1>& x) { return x; };
        const auto ctrl = [&](double, const ColVec<1>&) {
            ++calls;
            return ColVec<1>{1.0};
        };
        FixedStepSolver solver(RK4<1>{}, 0.05);
        const auto      sim = simulate_sampled<1, 1, 1, double>(plant, output, ctrl, solver, 0.1, ColVec<1>{0.0}, {0.0, 100.0});
        CHECK(calls == 1000);
        CHECK(sim.t.back() == 100.0);
    }

    TEST_CASE("simulate_multirate runs exactly one control step per base tick") {
        std::size_t calls = 0;
        const auto  plant = [](double, const ColVec<1>& x, const ColVec<1>& u) { return ColVec<1>{-x(0) + u(0)}; };
        const auto  output = [](const ColVec<1>& x) { return x; };
        const auto  step = [&](std::size_t, double, const ColVec<1>&) {
            ++calls;
            return ColVec<1>{1.0};
        };
        FixedStepSolver solver(RK4<1>{}, 0.05);
        const auto      sim = simulate_multirate<1, 1, 1, double>(plant, output, step, solver, 0.1, ColVec<1>{0.0}, {0.0, 100.0});
        CHECK(calls == 1000);
        CHECK(sim.t.back() == 100.0);
    }

    TEST_CASE("simulate_lti_siso records y = Cx + Du") {
        const StateSpace<1, 1, 1> sys{
            .A = Matrix<1, 1>{{-1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{0.0}},
            .D = Matrix<1, 1>{{1.0}},
            .Ts = 0.0,
        };
        ConstController c;
        FixedStepSolver solver(RK4<1>{}, 0.01);
        const auto      sim = simulate_lti_siso(sys, c, 1.0, solver, 0.1, ColVec<1>{0.0}, {0.0, 0.5});
        for (std::size_t j = 1; j < sim.t.size(); ++j) {
            CHECK(sim.y[j](0) == doctest::Approx(sim.u[j](0)));
        }
    }

    TEST_CASE("simulate_cached_zoh over a span never steps past tf") {
        const StateSpace<1, 1, 1> sys{.A = Matrix<1, 1>{{-1.0}}, .B = Matrix<1, 1>{{1.0}}, .C = Matrix<1, 1>{{1.0}}, .Ts = 0.0};
        const auto                sim = simulate_cached_zoh(sys, ColVec<1>{0.0}, ColVec<1>{1.0}, {0.0, 1.0}, 0.4);
        CHECK(sim.t.back() <= 1.0);
        CHECK(sim.t.size() == 3); // 0, 0.4, 0.8
        const auto exact = simulate_cached_zoh(sys, ColVec<1>{0.0}, ColVec<1>{1.0}, {0.0, 0.3}, 0.1);
        CHECK(exact.t.size() == 4); // 0.3/0.1 = 2.9999999999999996 still counts three steps
    }
}

TEST_SUITE("Review g6: hybrid") {
    const StateSpace<1, 1, 1> integ{.A = Matrix<1, 1>{{0.0}}, .B = Matrix<1, 1>{{1.0}}, .C = Matrix<1, 1>{{1.0}}, .Ts = 0.0};

    TEST_CASE("A non-advancing schedule is reported, not integrated open-loop to tf") {
        PiecewiseLTI<1, 1, 1, 1> plant{.modes = {integ}};
        const auto               schedule = [](double t, std::size_t, const auto&, const auto&) { return t; };
        const auto               on_event = [](double, std::size_t m, const auto& x, const auto&, const auto& u) {
            return HybridEventAction<1, 1, 1>{.mode = m, .u = u, .x = x};
        };
        const auto tr = simulate_piecewise_lti(plant, 0, ColVec<1>{0.0}, ColVec<1>{1.0}, {0.0, 1.0}, schedule, on_event);
        CHECK(tr.truncated);
    }

    TEST_CASE("FixedStride logs h_fine samples without drift duplicates") {
        PiecewiseLTI<1, 1, 1, 1> plant{.modes = {integ}};
        const auto               schedule = [](double, std::size_t, const auto&, const auto&) { return 1.0; };
        const auto               on_event = [](double, std::size_t m, const auto& x, const auto&, const auto& u) {
            return HybridEventAction<1, 1, 1>{.mode = m, .u = u, .x = x};
        };
        EventSimConfig<> cfg{.log = LogPolicy::FixedStride, .h_fine = 0.1};
        const auto       tr = simulate_piecewise_lti(plant, 0, ColVec<1>{0.0}, ColVec<1>{1.0}, {0.0, 1.0}, schedule, on_event, cfg);
        CHECK(tr.t.size() == 11);
        CHECK(tr.x.back()(0) == doctest::Approx(1.0));
    }

    TEST_CASE("locate_zero_crossing_exact default tolerance is relative to the segment") {
        // ẋ = u = 1e12: x = 1e12 t crosses 3 at t = 3e-12, on a picosecond segment.
        const ColVec<1> x0{0.0};
        const ColVec<1> u{1e12};
        const double    t = locate_zero_crossing_exact(integ, x0, u, 0.0, 1e-11, [](const ColVec<1>& x) { return x(0) - 3.0; });
        CHECK(std::abs(t - 3e-12) <= 1e-9 * 3e-12);
    }

    TEST_CASE("NPY header is padded to a 64-byte boundary") {
        const char* path = "tests/build/test_review_g6.npy";
        REQUIRE(write_npy_f64(path, {3}, {1.0, 2.0, 3.0}));
        FILE* f = std::fopen(path, "rb");
        REQUIRE(f != nullptr);
        unsigned char pre[10]{};
        REQUIRE(std::fread(pre, 1, 10, f) == 10);
        std::fclose(f);
        const std::size_t hdr_len = pre[8] | (pre[9] << 8);
        CHECK((10 + hdr_len) % 64 == 0);
    }
}

TEST_SUITE("Review g6: matlab aliases and plots") {
    TEST_CASE("matlab::reg builds the y→u compensator (separation principle)") {
        const Matrix<2, 2> A{{0.0, 1.0}, {2.0, -1.0}};
        const Matrix<2, 1> B{{0.0}, {1.0}};
        const Matrix<1, 2> C{{1.0, 0.0}};
        const Matrix<1, 1> D{{0.5}};
        const auto         sys = matlab::ss(A, B, C, D);
        const Matrix<1, 2> K{{8.0, 4.0}};    // eig(A − BK) = {−2, −3}
        const Matrix<2, 1> L{{9.0}, {17.0}}; // eig(A − LC) = {−4, −6}
        const auto         rsys = matlab::reg(sys, K, L);
        // Close the loop: plant u = C_r ξ + D_r y, compensator input y.
        Matrix<4, 4> Acl{};
        Acl.template block<2, 2>(0, 0) = A;
        Acl.template block<2, 2>(0, 2) = B * rsys.C;
        Acl.template block<2, 2>(2, 0) = rsys.B * C;
        Acl.template block<2, 2>(2, 2) = rsys.A + rsys.B * D * rsys.C;
        auto                ev = mat::compute_eigenvalues(Acl).values;
        std::vector<double> re;
        for (std::size_t i = 0; i < 4; ++i) {
            CHECK(std::abs(ev[i].imag()) < 1e-9);
            re.push_back(ev[i].real());
        }
        std::sort(re.begin(), re.end());
        CHECK(re[0] == doctest::Approx(-6.0));
        CHECK(re[1] == doctest::Approx(-4.0));
        CHECK(re[2] == doctest::Approx(-3.0));
        CHECK(re[3] == doctest::Approx(-2.0));
    }

    TEST_CASE("matlab::lqi / lqg refuse a continuous plant") {
        const auto sys = matlab::ss(Matrix<1, 1>{{-1.0}}, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}});
        CHECK_FALSE(matlab::lqi(sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{1.0}}).success);
    }

    TEST_CASE("matlab::pidstd default N is MATLAB's (no derivative filter)") {
        const auto r = matlab::pidstd(2.0, 1.0, 0.5);
        CHECK(r.Tf == 0.0);
    }

    TEST_CASE("matlab::lqr is fail-closed") {
        // Unstable, uncontrollable mode: no stabilizing gain exists.
        const Matrix<2, 2> A{{1.0, 0.0}, {0.0, -1.0}};
        const Matrix<2, 1> B{{0.0}, {1.0}};
        CHECK_FALSE(matlab::lqr(A, B, Matrix<2, 2>::identity(), Matrix<1, 1>{{1.0}}).has_value());
    }

    TEST_CASE("Nichols M-circle samples satisfy |G/(1+G)| = M") {
        plotlypp::Figure fig;
        plot::detail::add_nichols_m_circle(fig, 2.0, "M=2");
        const auto& data = fig.json()["data"][0];
        const auto& xs = data["x"];
        const auto& ys = data["y"];
        REQUIRE(xs.size() > 10);
        for (std::size_t i = 0; i < xs.size(); ++i) {
            const double               ph = xs[i].get<double>() * std::numbers::pi / 180.0;
            const double               mag = std::pow(10.0, ys[i].get<double>() / 20.0);
            const std::complex<double> G = std::polar(mag, ph);
            CHECK(std::abs(G / (1.0 + G)) == doctest::Approx(2.0).epsilon(1e-9));
        }
    }
}
