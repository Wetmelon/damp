// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include "damp/analysis/analysis.hpp"
#include "damp/backend.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/simulation/integrator.hpp"
#include "damp/simulation/plot_plotly.hpp"
#include "damp/simulation/solver.hpp"
#include "damp/systems/transfer_function.hpp"
#include "plotlypp/figure.hpp"
#include "plotlypp/layout/layout.hpp"
#include "plotlypp/traces/scatter.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::sim;

// Simple exponential decay: dx/dt = -x, x(0) = 1 => x(t) = exp(-t)
static auto exp_decay = [](double /*t*/, const ColVec<1>& x) -> ColVec<1> {
    return ColVec<1>{-x(0, 0)};
};

TEST_CASE("SolveResult iterator") {
    SolveResult<1> result;
    result.t = {0.0, 0.5, 1.0};
    result.x = {ColVec<1>{1.0}, ColVec<1>{0.6}, ColVec<1>{0.3}};

    size_t count = 0;
    for (const auto& [t, x] : result) {
        CHECK(t == result.t[count]);
        CHECK(x(0, 0) == result.x[count](0, 0));
        ++count;
    }
    CHECK(count == 3);
    CHECK(result.size() == 3);
}

TEST_CASE("FixedStepSolver - RK4 exponential decay") {
    RK4<1>          rk4;
    FixedStepSolver solver(rk4, 0.01);

    ColVec<1> x0{1.0};
    auto      result = solver.solve(exp_decay, x0, {0.0, 1.0});

    CHECK(result.success);
    CHECK(result.t.size() > 1);

    // x(1) = exp(-1) ≈ 0.3678794
    double x_final = result.x.back()(0, 0);
    CHECK(x_final == doctest::Approx(0.36787944117).epsilon(1e-8));
}

TEST_CASE("FixedStepSolver - matches manual RK4 loop") {
    // Manually integrate with RK4 and compare
    RK4<1>    rk4;
    double    h = 0.1;
    double    t = 0.0;
    ColVec<1> x{1.0};

    for (int i = 0; i < 10; ++i) {
        auto r = rk4.evolve(exp_decay, x, t, h);
        x = r.x;
        t += h;
    }
    double manual_result = x(0, 0);

    FixedStepSolver solver(rk4, 0.1);
    auto            result = solver.solve(exp_decay, ColVec<1>{1.0}, {0.0, 1.0});

    CHECK(result.x.back()(0, 0) == doctest::Approx(manual_result).epsilon(1e-14));
}

TEST_CASE("FixedStepSolver - stop condition") {
    RK4<1>          rk4;
    FixedStepSolver solver(rk4, 0.01);

    // Stop when x < 0.5
    solver.set_stop_condition([](double /*t*/, const ColVec<1>& x) {
        return x(0, 0) < 0.5;
    });

    ColVec<1> x0{1.0};
    auto      result = solver.solve(exp_decay, x0, {0.0, 5.0});

    // Should have stopped well before t=5
    CHECK(result.t.back() < 1.0);
    CHECK(result.x.back()(0, 0) < 0.5);
}

TEST_CASE("FixedStepSolver - 2D harmonic oscillator") {
    // dx/dt = [x1, -x0]  ==> x0(t)=cos(t), x1(t)=-sin(t)
    auto harmonic = [](double /*t*/, const ColVec<2>& x) -> ColVec<2> {
        return ColVec<2>{x(1, 0), -x(0, 0)};
    };

    RK4<2>          rk4;
    FixedStepSolver solver(rk4, 0.001);

    ColVec<2> x0{1.0, 0.0};
    auto      result = solver.solve(harmonic, x0, {0.0, 6.283185307});

    // After one full period, should return near initial condition
    CHECK(result.x.back()(0, 0) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(result.x.back()(1, 0) == doctest::Approx(0.0).epsilon(1e-6));

    // Plot x0=cos(t) and x1=-sin(t)
    auto t = plot::to_double_vector(result.t);
    auto x0_hist = plot::extract_channel(result.x, 0);
    auto x1_hist = plot::extract_channel(result.x, 1);

    using namespace plotlypp;
    auto fig = Figure()
                   .addTrace(Scatter().x(t).y(x0_hist).mode({Scatter::Mode::Lines}).name("x0 = cos(t)"))
                   .addTrace(Scatter().x(t).y(x1_hist).mode({Scatter::Mode::Lines}).name("x1 = -sin(t)"))
                   .setLayout(Layout().title([](auto& title) { title.text("Harmonic Oscillator"); }).xaxis(Layout::Xaxis().title([](auto& title) { title.text("Time (s)"); })).yaxis(Layout::Yaxis().title([](auto& title) { title.text("State"); })));
    plot::write_html(fig, "tests/build/harmonic_oscillator.html");
}

TEST_CASE("AdaptiveStepSolver - DP45 exponential decay") {
    DP45<1>            dp45;
    AdaptiveStepSolver solver(dp45, {.atol = 1e-8, .rtol = 1e-8, .first_step = 0.1});

    ColVec<1> x0{1.0};
    auto      result = solver.solve(exp_decay, x0, {0.0, 1.0});

    CHECK(result.success);
    CHECK(result.n_accepted > 0);
    CHECK(result.x.back()(0, 0) == doctest::Approx(0.36787944117).epsilon(1e-6));
}

TEST_CASE("adaptive_solve<DP45> one-liner matches AdaptiveStepSolver") {
    ColVec<1> x0{1.0};
    auto      one = adaptive_solve<DP45>(exp_decay, x0, {0.0, 1.0}, {.atol = 1e-8, .rtol = 1e-8, .first_step = 0.1});
    CHECK(one.success);
    CHECK(one.x.back()(0, 0) == doctest::Approx(0.36787944117).epsilon(1e-6));
}

TEST_CASE("adaptive_solve defaults match SciPy atol/rtol") {
    ColVec<1> x0{1.0};
    auto      result = adaptive_solve<DP45>(exp_decay, x0, {0.0, 1.0});
    CHECK(result.success);
    // SciPy defaults (rtol=1e-3, atol=1e-6) are looser than 1e-8 absolute.
    CHECK(result.x.back()(0, 0) == doctest::Approx(0.36787944117).epsilon(1e-3));
}

TEST_CASE("adaptive_solve multi-scale states use weighted error") {
    // x0 = [1e-6, 1e3]: pure absolute L2 control would be dominated by the large
    // component; weighted atol/rtol should still integrate both channels.
    auto multi = [](double /*t*/, const ColVec<2>& x) -> ColVec<2> {
        return ColVec<2>{-x(0, 0), -x(1, 0)};
    };
    ColVec<2> x0{1e-6, 1e3};
    auto      result = adaptive_solve<DP45>(multi, x0, {0.0, 1.0}, {.atol = 1e-9, .rtol = 1e-6, .first_step = 0.01});
    CHECK(result.success);
    CHECK(result.x.back()(0, 0) == doctest::Approx(1e-6 * std::exp(-1.0)).epsilon(1e-3));
    CHECK(result.x.back()(1, 0) == doctest::Approx(1e3 * std::exp(-1.0)).epsilon(1e-3));
}

TEST_CASE("adaptive_solve Hairer first_step when first_step is zero") {
    // No explicit first_step → select_initial_step from f(t0, x0).
    ColVec<1> x0{1.0};
    auto      result = adaptive_solve<DP45>(exp_decay, x0, {0.0, 1.0}, {.atol = 1e-8, .rtol = 1e-8});
    CHECK(result.success);
    CHECK(result.n_accepted > 0);
    CHECK(result.x.back()(0, 0) == doctest::Approx(0.36787944117).epsilon(1e-6));
}

TEST_CASE("adaptive_solve per-component atol/rtol") {
    // Tight abs tol on the small channel only; large channel uses loose abs + moderate rel.
    auto multi = [](double /*t*/, const ColVec<2>& x) -> ColVec<2> {
        return ColVec<2>{-x(0, 0), -x(1, 0)};
    };
    ColVec<2>          x0{1e-6, 1e3};
    AdaptiveOptions<2> opts{
        .atol = ColVec<2>{1e-12, 1e0},
        .rtol = ColVec<2>{1e-6, 1e-3},
        .first_step = 0.01,
    };
    auto result = adaptive_solve<DP45>(multi, x0, {0.0, 1.0}, opts);
    CHECK(result.success);
    CHECK(result.x.back()(0, 0) == doctest::Approx(1e-6 * std::exp(-1.0)).epsilon(1e-3));
    CHECK(result.x.back()(1, 0) == doctest::Approx(1e3 * std::exp(-1.0)).epsilon(1e-3));
}

TEST_CASE("adaptive_solve Gustafsson PI stays accurate on smooth problem") {
    // PI controller is the default after the first accept; check accuracy + efficiency.
    ColVec<1> x0{1.0};
    auto      result = adaptive_solve<DP45>(exp_decay, x0, {0.0, 2.0}, {.atol = 1e-8, .rtol = 1e-8});
    CHECK(result.success);
    CHECK(result.x.back()(0, 0) == doctest::Approx(std::exp(-2.0)).epsilon(1e-6));
    // Should not thrash: rejections << accepts on a smooth linear problem.
    CHECK(result.n_accepted > result.n_rejected);
}

TEST_CASE("AdaptiveStepSolver fails at min_step when error stays large") {
    // Stiff for explicit DP45: with a large min_step the local error cannot meet tol.
    auto stiff = [](double /*t*/, const ColVec<1>& x) -> ColVec<1> {
        return ColVec<1>{-1000.0 * x(0, 0)};
    };
    DP45<1>            dp45;
    AdaptiveStepSolver solver(dp45, AdaptiveOptions<1>{
                                        .atol = 1e-12,
                                        .rtol = 1e-12,
                                        .first_step = 0.1,
                                        .min_step = 0.05,
                                        .max_step = 0.1,
                                        .fail_on_min_step = true,
                                    });
    auto               result = solver.solve(stiff, ColVec<1>{1.0}, {0.0, 1.0});
    CHECK_FALSE(result.success);
}

TEST_CASE("adaptive_solve<TRBDF2> handles stiff decay") {
    // dx/dt = -1000 x: explicit fixed-step needs h ≪ 0.001; adaptive TR-BDF2
    // should take large stable steps and land near e^{-10} at t = 0.01.
    auto stiff = [](double /*t*/, const ColVec<1>& x) -> ColVec<1> {
        return ColVec<1>{-1000.0 * x(0, 0)};
    };
    ColVec<1> x0{1.0};
    auto      result = adaptive_solve<TRBDF2>(stiff, x0, {0.0, 0.01}, {.atol = 1e-6, .rtol = 1e-6, .first_step = 1e-4});
    CHECK(result.success);
    CHECK(result.x.back()(0, 0) == doctest::Approx(std::exp(-10.0)).epsilon(1e-3));
}

TEST_CASE("fixed_solve<RK4> one-liner matches FixedStepSolver") {
    ColVec<1> x0{1.0};
    auto      one = fixed_solve<RK4>(exp_decay, x0, {0.0, 1.0}, 0.01);
    CHECK(one.success);
    CHECK(one.x.back()(0, 0) == doctest::Approx(0.36787944117).epsilon(1e-8));
}

TEST_CASE("AdaptiveStepSolver - fewer steps than fixed for smooth problem") {
    DP45<1>            dp45;
    AdaptiveStepSolver adaptive(dp45, {.atol = 1e-6, .rtol = 1e-6, .first_step = 0.1});

    RK4<1>          rk4;
    FixedStepSolver fixed_solver(rk4, 0.001);

    ColVec<1> x0{1.0};
    auto      result_adaptive = adaptive.solve(exp_decay, x0, {0.0, 2.0});
    auto      result_fixed = fixed_solver.solve(exp_decay, x0, {0.0, 2.0});

    // Adaptive should use fewer evolve calls for same accuracy
    CHECK(result_adaptive.nfev < result_fixed.nfev);
    CHECK(result_adaptive.n_accepted > 0);

    // Both should converge to the same answer
    CHECK(result_adaptive.x.back()(0, 0) == doctest::Approx(result_fixed.x.back()(0, 0)).epsilon(1e-5));
}

TEST_CASE("AdaptiveStepSolver - zero-crossing detection") {
    // dx/dt = 1 (linear ramp: x(t) = t - 0.5, crosses zero at t=0.5)
    auto ramp = [](double /*t*/, const ColVec<1>& /*x*/) -> ColVec<1> {
        return ColVec<1>{1.0};
    };

    DP45<1>            dp45;
    AdaptiveStepSolver solver(dp45, {.atol = 1e-8, .rtol = 1e-8, .first_step = 0.1});

    // Detect when x crosses zero
    solver.add_zero_crossing([](double /*t*/, const ColVec<1>& x) -> double {
        return x(0, 0);
    });

    ColVec<1> x0{-0.5};
    auto      result = solver.solve(ramp, x0, {0.0, 1.0});

    // Should have a point very close to t=0.5, x=0
    bool found_crossing = false;
    for (const auto& [t, x] : result) {
        if (std::abs(x(0, 0)) < 1e-6 && std::abs(t - 0.5) < 0.01) {
            found_crossing = true;
            break;
        }
    }
    CHECK(found_crossing);
}

TEST_CASE("FixedStepSolver - Van der Pol oscillator") {
    // Van der Pol oscillator: x'' - mu*(1-x^2)*x' + x = 0
    // State form: x0' = x1, x1' = mu*(1 - x0^2)*x1 - x0
    constexpr double mu = 1.0;
    auto             vdp = [](double /*t*/, const ColVec<2>& x) -> ColVec<2> {
        return ColVec<2>{x(1, 0), (mu * (1.0 - (x(0, 0) * x(0, 0))) * x(1, 0)) - x(0, 0)};
    };

    RK4<2>          rk4;
    FixedStepSolver solver(rk4, 0.001);

    ColVec<2> x0{2.0, 0.0};
    auto      result = solver.solve(vdp, x0, {0.0, 30.0});

    CHECK(result.success);

    // The Van der Pol oscillator converges to a limit cycle with amplitude ~2
    // After 30s the trajectory should be on the limit cycle
    double x_max = 0.0;
    // Check last ~10 seconds worth of data for peak amplitude
    size_t start_idx = result.t.size() * 2 / 3;
    for (size_t i = start_idx; i < result.t.size(); ++i) {
        double val = std::abs(result.x[i](0, 0));
        x_max = damp::max(val, x_max);
    }
    // Limit cycle peak amplitude is slightly above 2 for mu=1
    CHECK(x_max == doctest::Approx(2.009).epsilon(0.01));

    // Plot time history and phase portrait as stacked subplots
    auto t = plot::to_double_vector(result.t);
    auto x0_hist = plot::extract_channel(result.x, 0);
    auto x1_hist = plot::extract_channel(result.x, 1);

    using namespace plotlypp;

    auto fig = Figure()
                   // Time history (top)
                   .addTrace(Scatter().x(t).y(x0_hist).mode({Scatter::Mode::Lines}).name("x (position)"))
                   .addTrace(Scatter().x(t).y(x1_hist).mode({Scatter::Mode::Lines}).name("x' (velocity)"))
                   // Phase portrait (bottom)
                   .addTrace(Scatter().x(x0_hist).y(x1_hist).mode({Scatter::Mode::Lines}).name("Trajectory").xaxis("x2").yaxis("y2"))
                   .setLayout(Layout().title([](auto& title) { title.text("Van der Pol Oscillator (mu=1)"); }).grid(Layout::Grid().rows(2).columns(1).pattern(Layout::Grid::Pattern::Independent).roworder(Layout::Grid::Roworder::TopToBottom)).xaxis(Layout::Xaxis().title([](auto& title) { title.text("Time (s)"); })).yaxis(Layout::Yaxis().title([](auto& title) { title.text("State"); })).xaxis(2, Layout::Xaxis().title([](auto& title) { title.text("x"); })).yaxis(2, Layout::Yaxis().title([](auto& title) { title.text("x'"); })));
    plot::write_html(fig, "tests/build/vanderpol.html");
}

TEST_CASE("MATLAB®-style plot wrappers instantiate and render") {
    // Smoke test: each wrapper must compile against the analysis result types
    // and emit HTML carrying its title. 1/(s^2+s+1) exercises poles + a Bode/Nyquist sweep.
    TransferFunction<1, 3, double> tf{{1.0}, {1.0, 1.0, 1.0}};
    const auto                     ss = tf.to_state_space().value();
    const auto                     t = analysis::linspace(0.0, 5.0, 51);
    const std::vector<double>      u(t.size(), 1.0);
    const auto                     omega = analysis::logspace(0.1, 100.0, 50);

    const auto s = analysis::step(ss, t);
    const auto im = analysis::impulse(ss, t);
    const auto ls = analysis::lsim(ss, u, t);
    const auto bd = analysis::bode(ss, omega);
    const auto nq = analysis::nyquist(ss, omega);
    const auto pz = analysis::pzmap(tf);
    const auto nc = analysis::nichols(ss, omega);
    const auto sg = analysis::sigma(ss, omega);
    const auto rl = analysis::rlocus(tf, analysis::linspace(0.0, 10.0, 21));

    auto carries = [](const plotlypp::Figure& f, const std::string& needle) {
        std::ostringstream os;
        f.toHtml(os);
        return os.str().find(needle) != std::string::npos;
    };

    CHECK(carries(plot::stepplot(s), "Step Response"));
    CHECK(carries(plot::impulseplot(im), "Impulse Response"));
    CHECK(carries(plot::lsimplot(ls), "Linear Simulation"));
    CHECK(carries(plot::bodeplot(bd), "Bode Plot"));
    CHECK(carries(plot::bodemag(bd), "Bode Magnitude"));
    CHECK(carries(plot::nyquistplot(nq), "Nyquist Plot"));
    CHECK(carries(plot::pzplot(pz), "Pole-Zero Map"));
    CHECK(carries(plot::nicholsplot(nc), "Nichols Chart"));
    CHECK(carries(plot::sigmaplot(sg), "Singular Values"));
    CHECK(carries(plot::rlocusplot(rl), "Root Locus"));
}
