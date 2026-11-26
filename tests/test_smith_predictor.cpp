// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_smith_predictor.cpp
 * @brief Smith predictor: design, L=0 equivalence to PID, dead-time closed loop,
 *        a fraction of a sample of delay, custom SisoLaw primary, handover
 */

#include <cmath>
#include <cstddef>
#include <limits>

#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/controllers/smith_predictor.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/design/status.hpp"
#include "damp/filters/delay.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/toolbox/bounds.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

/// ZOH first-order plant step: x <- a x + b u, y = x (before update)
struct FirstOrderPlant {
    double a{};
    double b{};
    double x{0.0};

    double step(double u) {
        const double y = x;
        x = (a * x) + (b * u);
        return y;
    }

    void reset() { x = 0.0; }
};

/// First-order + pure delay plant (FOPDT discrete)
template<size_t MaxD>
struct FopdtPlant {
    FirstOrderPlant     fo{};
    Delay<MaxD, double> delay{};
    size_t              d{0};

    void init(double K, double tau, double L, double Ts) {
        fo.a = std::exp(-Ts / tau);
        fo.b = K * (1.0 - fo.a);
        d = static_cast<size_t>((L / Ts) + 0.5);
        delay.init(d);
        reset();
    }

    double step(double u) {
        // Delay-free FO steps on u; measurement is delayed output
        const double y_undelayed = fo.step(u);
        return delay(y_undelayed);
    }

    void reset() {
        fo.reset();
        delay.reset();
    }
};

/// Minimal custom primary to prove Primary is not PID-only: a proportional law on the SISO contract
struct ProportionalPrimary {
    /// Its limits: the command's
    struct Limits {
        Bounds<1, double> input{};
    };

    double            k{1.0};
    Bounds<1, double> lim{};
    double            u_cmd{0.0};

    [[nodiscard]] constexpr double step(double r, double y) {
        u_cmd = lim.saturate(k * (r - y));
        return u_cmd;
    }

    constexpr void back_calculate(double /*u_sat*/, double /*u_ff*/ = 0.0) {}

    constexpr void track(double u, double /*r*/, double /*y*/) { u_cmd = u; }

    [[nodiscard]] constexpr Limits limits() const { return {lim}; }

    [[nodiscard]] constexpr bool set_limits(const Limits& l) {
        if (!l.input.valid()) {
            return false;
        }
        lim = l.input;
        return true;
    }

    [[nodiscard]] constexpr double last_command() const { return u_cmd; }

    [[nodiscard]] constexpr bool valid() const { return true; }

    constexpr void reset() { u_cmd = 0.0; }
};

/// The Smith predictor as it stood before its model ran on the delayed command: the delay-free FO
/// model x and a Delay line of its output. The oracle for the rebuilt law.
template<size_t MaxDelay, typename Primary>
struct LegacySmith {
    Primary                 primary{};
    Delay<MaxDelay, double> delay{};
    double                  a{0.0};
    double                  b{0.0};
    double                  x{0.0};

    LegacySmith(const Primary& c, double a_, double b_, size_t d) : primary(c), a(a_), b(b_) { delay.init(d); }

    double step(double r, double y) {
        const double y_hat = x;
        const double y_delayed = delay(y_hat);
        const double u = primary.step(r, y + y_hat - y_delayed);
        x = (a * x) + (b * u);
        return u;
    }
};

/// Largest |u_new − u_old| over @p N ticks of the rebuilt and the legacy law, each in closed loop
/// with its own copy of the FOPDT plant (K, τ, L at Ts), reference steps at 0 and N/2.
template<size_t MaxD, typename Primary>
double max_command_gap(const Primary& c, double K, double tau, double L, double Ts, size_t N) {
    const auto res = design::smith_predictor(design::FOPDT{.K = K, .tau = tau}, Ts, InputDelay<MaxD>{L});
    REQUIRE(res.success);
    const double                          a = std::exp(-Ts / tau);
    const double                          b = K * (1.0 - a);
    const auto                            d = static_cast<size_t>((L / Ts) + 0.5);
    SmithPredictor<MaxD, double, Primary> rebuilt(c, res);
    LegacySmith<MaxD + 1, Primary>        legacy(c, a, b, d);
    FopdtPlant<MaxD + 1>                  plant_new{};
    FopdtPlant<MaxD + 1>                  plant_old{};
    plant_new.init(K, tau, L, Ts);
    plant_old.init(K, tau, L, Ts);
    double y_new = 0.0;
    double y_old = 0.0;
    double gap = 0.0;
    for (size_t k = 0; k < N; ++k) {
        const double r = (k < N / 2) ? 1.0 : -0.5;
        const double u_new = rebuilt.step(r, y_new);
        const double u_old = legacy.step(r, y_old);
        gap = std::max(gap, std::abs(u_new - u_old));
        y_new = plant_new.step(u_new);
        y_old = plant_old.step(u_old);
    }
    return gap;
}

static_assert(SisoLaw<ProportionalPrimary, double>);
static_assert(SisoLaw<SmithPredictor<8, double>, double>);
static_assert(SisoLaw<SmithPredictor<8, float, PIController<float>>, float>);
static_assert(SisoLaw<SmithPredictor<8, double, ProportionalPrimary>, double>);

} // namespace

TEST_SUITE("Smith Predictor") {

    TEST_CASE("design::smith_predictor packs the FO model, the delay and a PI on the delay-free plant") {
        constexpr double K = 1.0;
        constexpr double tau = 5.0;
        constexpr double L = 2.0;
        constexpr double Ts = 0.1;
        constexpr double tau_c = 2.0;

        constexpr auto sp = design::smith_predictor(design::FOPDT{.K = K, .tau = tau}, Ts, InputDelay<20>{L}, design::SIMC{.tau_c = tau_c});
        static_assert(sp.success);
        static_assert(sp.delay_samples == 20); // 2.0/0.1
        static_assert(sp.b_stale == 0.0);      // whole samples

        const double a_ref = std::exp(-Ts / tau);
        const double b_ref = K * (1.0 - a_ref);
        CHECK(sp.a == doctest::Approx(a_ref).epsilon(1e-12));
        CHECK(sp.b == doctest::Approx(b_ref).epsilon(1e-12));
        CHECK(sp.b_fresh == sp.b);
        CHECK(sp.delay == L);
        CHECK(sp.pid.Ts == doctest::Approx(Ts));
        CHECK(sp.pid.gains.Kp > 0.0);
        CHECK(sp.pid.gains.Ki > 0.0); // PI from SIMC

        // SIMC on delay-free plant: Kp = tau / (K * tau_c). The packed PIDResult carries the
        // continuous gains and Ts (the discrete Ki·Ts is the runtime's).
        const auto cont = design::pid(design::FOPDT{.K = K, .L = 0.0, .tau = tau}, Ts, design::SIMC{.tau_c = tau_c}).gains;
        CHECK(sp.pid.gains.Kp == doctest::Approx(cont.Kp));
        CHECK(sp.pid.gains.Ki == doctest::Approx(cont.Ki));
        CHECK(sp.pid.gains.Kp == doctest::Approx(tau / (K * tau_c)));

        // Without a rule: SIMC's PI with τ_c = τ, as the old τ_c ≤ 0 default.
        constexpr auto by_default = design::smith_predictor(design::FOPDT{.K = K, .tau = tau}, Ts, InputDelay<20>{L});
        CHECK(by_default.pid.gains.Kp == doctest::Approx(tau / (K * tau)));
    }

    TEST_CASE("the design refuses a non-physical plant, a dead time in the model, and a delay past N·Ts") {
        CHECK_FALSE(design::smith_predictor(design::FOPDT{.K = 0.0, .tau = 5.0}, 0.1, InputDelay<16>{1.0}).success);
        CHECK_FALSE(design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 0.0}, 0.1, InputDelay<16>{1.0}).success);
        CHECK(design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 5.0}, 0.0, InputDelay<16>{1.0}).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 5.0}, 0.1, InputDelay<16>{-1.0}).status == design::DesignStatus::invalid_parameter);
        //! The dead time is the InputDelay's: an FOPDT that carries one too is refused, not added.
        CHECK(design::smith_predictor(design::FOPDT{.K = 1.0, .L = 1.0, .tau = 5.0}, 0.1, InputDelay<16>{1.0}).status == design::DesignStatus::invalid_parameter);
        //! N delay states hold N·Ts: θ = 0.8 fits N = 8 at 0.1 s, θ = 0.9 does not.
        CHECK(design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 1.0}, 0.1, InputDelay<8>{0.8}, design::SIMC{.tau_c = 1.0}).success);
        CHECK(design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 1.0}, 0.1, InputDelay<8>{0.9}, design::SIMC{.tau_c = 1.0}).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("as<float> preserves success, the delay and N") {
        constexpr auto sp = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 3.0}, 0.05, InputDelay<32>{1.5}, design::SIMC{.tau_c = 1.0});
        static_assert(sp.success);
        constexpr auto spf = sp.as<float>();
        static_assert(spf.success);
        static_assert(spf.delay_samples == sp.delay_samples);
        static_assert(std::is_same_v<decltype(spf), const design::SmithPredictorResult<float, 32>>);
        CHECK(static_cast<double>(spf.a) == doctest::Approx(sp.a).epsilon(1e-6));
    }

    TEST_CASE("gains you have as the primary: Parallel, with the PID's options") {
        constexpr auto sp = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 4.0}, 0.1, InputDelay<10>{1.0}, design::Parallel{.Kp = 2.0, .Ki = 0.4}, design::InputLimits{.lower = {-1.0}, .upper = {1.0}});
        static_assert(sp.success);
        CHECK(sp.pid.gains.Kp == doctest::Approx(2.0));
        CHECK(sp.pid.gains.Ki == doctest::Approx(0.4)); // continuous, at sp.pid.Ts = 0.1
        CHECK(sp.pid.Ts == doctest::Approx(0.1));
        CHECK(sp.pid.limits.input.upper(0) == 1.0);
        CHECK(sp.delay_samples == 10);
        SmithPredictor smith{sp};
        CHECK(smith.limits().input.upper(0) == 1.0); // the primary's limits are the predictor's
        CHECK(smith.step(10.0, 0.0) == 1.0);
    }

    TEST_CASE("L=0 Smith predictor matches bare PI on FO plant") {
        constexpr double K = 1.0;
        constexpr double tau = 2.0;
        constexpr double Ts = 0.05;
        constexpr double tau_c = 1.0;

        constexpr auto sp = design::smith_predictor(design::FOPDT{.K = K, .tau = tau}, Ts, InputDelay<8>{0.0}, design::SIMC{.tau_c = tau_c});
        static_assert(sp.success);
        static_assert(sp.delay_samples == 0);

        SmithPredictor<8, double> smith(sp);
        PIController<double>      bare(sp.pid);

        FirstOrderPlant plant{.a = sp.a, .b = sp.b, .x = 0.0};
        FirstOrderPlant plant_bare{.a = sp.a, .b = sp.b, .x = 0.0};

        constexpr double r = 1.0;
        constexpr size_t N = 400;
        double           max_diff = 0.0;
        double           y = 0.0;
        double           y_bare = 0.0;
        for (size_t k = 0; k < N; ++k) {
            const double u = smith.step(r, y);
            const double u_b = bare.step(r, y_bare);
            max_diff = std::max(max_diff, std::abs(u - u_b));
            y = plant.step(u);
            y_bare = plant_bare.step(u_b);
        }

        // With L=0, y_pred = y so the primary PI is identical; model path cancels.
        CHECK(max_diff < 1e-9);
        CHECK(std::abs(y - r) < 0.05); // settled near setpoint
    }

    TEST_CASE("with dead time Smith tracks better than aggressive bare PI") {
        // Plant: G(s) = e^{-3s}/(s+1). PI designed for the delay-free plant with
        // tight tau_c; raw feedback of the delayed measurement is poorly damped,
        // while Smith restores the delay-free loop (output still waits ~L).
        constexpr double K = 1.0;
        constexpr double tau = 1.0;
        constexpr double L = 3.0;
        constexpr double Ts = 0.05;
        constexpr double tau_c = 0.5; // aggressive on the FO plant alone

        constexpr auto sp = design::smith_predictor(design::FOPDT{.K = K, .tau = tau}, Ts, InputDelay<128>{L}, design::SIMC{.tau_c = tau_c});
        REQUIRE(sp.success);
        REQUIRE(sp.delay_samples == 60);

        SmithPredictor       smith(sp);
        PIController<double> bare(sp.pid);
        static_assert(std::is_same_v<decltype(smith), SmithPredictor<128, double>>);

        FopdtPlant<128> plant_s{};
        FopdtPlant<128> plant_b{};
        plant_s.init(K, tau, L, Ts);
        plant_b.init(K, tau, L, Ts);

        constexpr double r = 1.0;
        constexpr size_t N = 800; // 40 s
        double           y_s = 0.0;
        double           y_b = 0.0;
        double           iae_s = 0.0;
        double           iae_b = 0.0;
        double           max_y_s = 0.0;
        double           max_y_b = 0.0;
        double           y_s_final = 0.0;
        double           peak_abs_e_b = 0.0;

        for (size_t k = 0; k < N; ++k) {
            const double u_s = smith.step(r, y_s);
            const double u_b = bare.step(r, y_b);
            y_s = plant_s.step(u_s);
            y_b = plant_b.step(u_b);
            iae_s += std::abs(r - y_s) * Ts;
            iae_b += std::abs(r - y_b) * Ts;
            max_y_s = std::max(max_y_s, y_s);
            max_y_b = std::max(max_y_b, y_b);
            peak_abs_e_b = std::max(peak_abs_e_b, std::abs(r - y_b));
            if (k + 100 >= N) {
                y_s_final = y_s;
            }
        }

        // Smith settles near the reference with modest overshoot
        CHECK(std::abs(y_s_final - r) < 0.05);
        CHECK(max_y_s < 1.25);
        // Same gains without Smith: larger overshoot / IAE on the delayed plant
        CHECK(max_y_b > max_y_s);
        CHECK(iae_s < iae_b);
        CHECK(iae_s < 0.75 * iae_b);
        // Bare loop is not well damped for this L/tau with delay-free tuning
        CHECK(peak_abs_e_b > 0.2);
    }

    TEST_CASE("a delay of a fraction of a sample: the model is the plant's, exactly") {
        //! The delay was rounded to whole samples, L/Ts = 2.5 to 3 here, and the model's ŷ_d
        //! disagreed with the plant's output by the half sample, in every transient.
        constexpr double K = 2.0;
        constexpr double tau = 0.4;
        constexpr double Ts = 0.1;
        constexpr double theta = 0.25; // 2.5 samples
        constexpr auto   sp = design::smith_predictor(design::FOPDT{.K = K, .tau = tau}, Ts, InputDelay<4>{theta}, design::SIMC{.tau_c = 0.3});
        REQUIRE(sp.success);
        CHECK(sp.delay_samples == 2);
        CHECK(sp.b_stale > 0.0);
        CHECK(sp.b_fresh + sp.b_stale == doctest::Approx(sp.b).epsilon(1e-15));

        // The plant: K/(τs + 1) with the delay, as discretize samples it.
        const StateSpace<1, 1, 1> plant{.A = {{-1.0 / tau}}, .B = {{K / tau}}, .C = {{1.0}}};
        const auto                delayed = discretize(plant, Ts, InputDelay<4>{theta});
        REQUIRE(delayed.success);
        SmithPredictor<4, double> smith{sp};
        ColVec<5, double>         x{};
        double                    worst = 0.0;
        for (int k = 0; k < 200; ++k) {
            const double y = (delayed->C * x)(0, 0);
            worst = std::max(worst, std::abs(y - smith.delayed_model_output()));
            const double u = smith.step((k < 100) ? 1.0 : -0.5, y);
            x = (delayed->A * x) + (delayed->B * ColVec<1, double>{u});
        }
        CHECK(worst < 1e-12); // ŷ_d is the plant's output: y − ŷ_d = 0
    }

    TEST_CASE("reset clears model and delay so next step is cold-start") {
        constexpr auto             sp = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 2.0}, 0.1, InputDelay<16>{0.5}, design::SIMC{.tau_c = 1.0});
        SmithPredictor<16, double> ctrl(sp);
        REQUIRE(ctrl.valid());

        double          y = 0.0;
        FirstOrderPlant plant{.a = sp.a, .b = sp.b};
        // Drive some state
        for (int i = 0; i < 30; ++i) {
            const double u = ctrl.step(1.0, y);
            y = plant.step(u);
        }
        CHECK(std::abs(ctrl.model_output()) > 1e-6);

        ctrl.reset();
        CHECK(ctrl.model_output() == doctest::Approx(0.0));
        CHECK(ctrl.last_command() == 0.0);
        // First post-reset command with y=0,r=0 should be ~0 (no integral yet)
        CHECK(ctrl.step(0.0, 0.0) == doctest::Approx(0.0));
    }

    TEST_CASE("an invalid design yields zero control; a design of fewer delay states fits a larger runtime") {
        SmithPredictor<8, double> bad(design::smith_predictor(design::FOPDT{.K = 0.0, .tau = 1.0}, 0.1, InputDelay<8>{0.0}));
        CHECK_FALSE(bad.valid());
        CHECK(bad.step(1.0, 0.0) == doctest::Approx(0.0));

        constexpr auto             d_eq = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 1.0}, 0.1, InputDelay<8>{0.8}, design::SIMC{.tau_c = 1.0});
        SmithPredictor<16, double> roomy(d_eq); // N = 8 into 16
        SmithPredictor             exact(d_eq);
        CHECK(roomy.valid());
        for (int k = 0; k < 50; ++k) {
            CHECK(roomy.step(1.0, 0.0) == exact.step(1.0, 0.0));
        }
    }

    TEST_CASE("the law on its two models reproduces the Delay-line law's commands") {
        // The same plants as the cases above. ŷ_d is now the model run on the delayed command,
        // which is the delay-free model delayed for whole samples: the two agree to the rounding
        // of the multiply-adds on commands of order 1–10, 1e-12 absolute.
        const auto pi = [](double K, double tau, double tau_c, double Ts) {
            return PIController<double>(design::pid(design::FOPDT{.K = K, .L = 0.0, .tau = tau}, Ts, design::SIMC{.tau_c = tau_c}));
        };
        const double gap_l0 = max_command_gap<8>(pi(1.0, 2.0, 1.0, 0.05), 1.0, 2.0, 0.0, 0.05, 400);
        const double gap_l3 = max_command_gap<128>(pi(1.0, 1.0, 0.5, 0.05), 1.0, 1.0, 3.0, 0.05, 800);
        const double gap_l05 = max_command_gap<16>(pi(1.0, 2.0, 1.0, 0.1), 1.0, 2.0, 0.5, 0.1, 400);
        const double gap_p = max_command_gap<32>(ProportionalPrimary{.k = 1.5}, 1.0, 2.0, 1.0, 0.05, 400);
        CHECK(gap_l0 < 1e-12); // d = 0: ŷ − ŷ_d is zero in both
        CHECK(gap_l3 < 1e-12);
        CHECK(gap_l05 < 1e-12);
        CHECK(gap_p < 1e-12);
    }

    TEST_CASE("custom SisoLaw primary (not PID-family)") {
        constexpr double K = 1.0;
        constexpr double tau = 2.0;
        constexpr double Ts = 0.05;
        const auto       model = design::smith_predictor(design::FOPDT{.K = K, .tau = tau}, Ts, InputDelay<8>{0.0});
        REQUIRE(model.success);

        // L = 0: y_pred = y, so Smith + pure P matches bare P on the FO plant
        ProportionalPrimary primary{.k = 1.5};
        SmithPredictor      smith(primary, model);
        static_assert(std::is_same_v<decltype(smith), SmithPredictor<8, double, ProportionalPrimary>>);
        REQUIRE(smith.valid());

        ProportionalPrimary bare{.k = 1.5};
        FirstOrderPlant     plant{.a = model.a, .b = model.b};
        FirstOrderPlant     plant_bare{.a = model.a, .b = model.b};

        constexpr double r = 1.0;
        double           y = 0.0;
        double           y_bare = 0.0;
        double           max_diff = 0.0;
        for (size_t k = 0; k < 200; ++k) {
            const double u = smith.step(r, y);
            const double u_b = bare.step(r, y_bare);
            max_diff = std::max(max_diff, std::abs(u - u_b));
            y = plant.step(u);
            y_bare = plant_bare.step(u_b);
        }
        // L=0 residual is identity: Smith+P and bare P must match sample-for-sample.
        // Pure P on a type-0 plant does not eliminate steady-state error — only compare u.
        CHECK(max_diff < 1e-12);
        CHECK(smith.primary().k == doctest::Approx(1.5));
        //! The limits are the primary's.
        REQUIRE(smith.set_limits(ProportionalPrimary::Limits{{.lower = -0.1, .upper = 0.1}}));
        CHECK(smith.primary().lim.upper(0) == 0.1);
        CHECK(smith.step(10.0, 0.0) == 0.1);
    }

    TEST_CASE("full PID primary via template parameter") {
        // Smith SIMC hides the delay from C, so the packed primary is a PI; the point of this
        // case is that Primary can be PIDController (full mode), not only PI.
        constexpr auto sp = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 2.0}, 0.05, InputDelay<8>{0.0}, design::SIMC{.tau_c = 1.0});
        REQUIRE(sp.success);
        SmithPredictor<8, double, PIDController<double>> smith(sp);
        REQUIRE(smith.valid());
        CHECK(std::isfinite(smith.step(1.0, 0.0)));
    }

    TEST_CASE("as<float> converting ctor preserves model coeffs") {
        constexpr auto             sp = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 3.0}, 0.1, InputDelay<16>{0.55}, design::SIMC{.tau_c = 1.0});
        SmithPredictor<16, double> src(sp);
        REQUIRE(src.valid());
        for (int k = 0; k < 10; ++k) {
            (void)src.step(1.0, 0.0);
        }
        SmithPredictor<16, float> dst(src);
        REQUIRE(dst.valid());
        CHECK(static_cast<double>(dst.model_output()) == doctest::Approx(src.model_output()).epsilon(1e-6));
        CHECK(static_cast<double>(dst.model_gains().stale) == doctest::Approx(src.model_gains().stale).epsilon(1e-6));
        CHECK(dst.delay_samples() == src.delay_samples());
        CHECK(static_cast<double>(dst.step(1.0f, 0.0f)) == doctest::Approx(src.step(1.0, 0.0)).epsilon(1e-5));
    }

    TEST_CASE("constexpr constructibility of design result") {
        constexpr auto sp = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 1.0}, 0.05, InputDelay<8>{0.2}, design::SIMC{.tau_c = 0.5});
        static_assert(sp.success);
        static_assert(sp.delay_samples == 4);
        constexpr auto disc = sp.pid;
        static_assert(disc.Ts > 0.0);
        (void)disc;
        // The runtime ticks under constant evaluation too (the commands in flight are a fixed array).
        constexpr double u2 = [] {
            SmithPredictor<8, double> c(design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 1.0}, 0.05, InputDelay<8>{0.2}, design::SIMC{.tau_c = 0.5}));
            (void)c.step(1.0, 0.0);
            return c.step(1.0, 0.0);
        }();
        static_assert(u2 > 0.0);
    }

    TEST_CASE("to_tf recovers SIMC primary; plant_tf is K/(τs+1)") {
        constexpr auto sp = design::smith_predictor(design::FOPDT{.K = 2.0, .tau = 5.0}, 0.1, InputDelay<16>{1.0}, design::SIMC{.tau_c = 2.0});
        REQUIRE(sp.success);
        const auto C = design::pid(design::FOPDT{.K = 2.0, .L = 0.0, .tau = 5.0}, 0.1, design::SIMC{.tau_c = 2.0}).gains;
        const auto tf = sp.to_tf();
        const auto cref = C.to_tf();
        CHECK(tf.num[1] == doctest::Approx(cref.num[1])); // Kp
        CHECK(tf.num[2] == doctest::Approx(cref.num[2])); // Ki
        CHECK(tf.den[1] == doctest::Approx(cref.den[1]));
        const auto Ghat = sp.plant_tf();
        CHECK(Ghat.num[0] == doctest::Approx(2.0));
        CHECK(Ghat.den[0] == doctest::Approx(5.0)); // τ
        CHECK(Ghat.den[1] == doctest::Approx(1.0));
        REQUIRE(sp.to_ss().has_value());
        CHECK(sp.to_ss()->D(0, 0) == doctest::Approx(C.to_ss()->D(0, 0)));
    }

    TEST_CASE("smith designs carry a failed primary's status and refuse a delay that is not finite") {
        // A negative tracking time is refused by the PID design; the Smith design used to report success.
        const auto from_bad_tt = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 4.0}, 0.1, InputDelay<16>{1.0}, design::Parallel{.Kp = 2.0, .Ki = 0.4}, design::AntiWindup{.tracking_time = -1.0});
        CHECK_FALSE(from_bad_tt.success);
        CHECK(from_bad_tt.status == design::DesignStatus::invalid_parameter);
        SmithPredictor<16, double> runtime(from_bad_tt);
        CHECK_FALSE(runtime.valid());

        const double inf = std::numeric_limits<double>::infinity();
        const double nan = std::numeric_limits<double>::quiet_NaN();
        CHECK(design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 4.0}, 0.1, InputDelay<16>{inf}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 4.0}, 0.1, InputDelay<16>{1e300}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::smith_predictor(design::FOPDT{.K = 1.0, .tau = nan}, 0.1, InputDelay<16>{1.0}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 4.0}, nan, InputDelay<16>{1.0}).status == design::DesignStatus::invalid_sample_time);
    }

    TEST_CASE("back_calculate: the model runs on the input that reached the plant, and the primary holds") {
        //! A clamp after the law could not be reported: the model stepped on the unclamped command
        //! and predicted an output the plant never produced, and the primary wound up behind it.
        constexpr double K = 1.0;
        constexpr double tau = 1.0;
        constexpr double L = 1.0;
        constexpr double Ts = 0.05;
        constexpr double rail = 0.8;
        constexpr auto   sp = design::smith_predictor(design::FOPDT{.K = K, .tau = tau}, Ts, InputDelay<64>{L}, design::SIMC{.tau_c = 0.5});
        REQUIRE(sp.success);
        SmithPredictor<64, double> told(sp);
        SmithPredictor<64, double> untold(sp);
        FopdtPlant<64>             plant_t{};
        FopdtPlant<64>             plant_u{};
        plant_t.init(K, tau, L, Ts);
        plant_u.init(K, tau, L, Ts);
        double y_t = 0.0;
        double y_u = 0.0;
        int    cut = 0;
        for (int k = 0; k < 1000; ++k) {
            const double r = (k < 600) ? 2.0 : 0.5; // out of the rail's reach, then within it
            const double asked = told.step(r, y_t);
            const double u_t = damp::min(asked, rail);
            told.back_calculate(u_t);
            if (asked > rail) {
                ++cut;
                CHECK(told.last_command() == asked);                        // the command, not the input reported
                CHECK(told.commands_in_flight()[told.ring_head()] == rail); // which the models run on
            }
            const double u_u = damp::min(untold.step(r, y_u), rail);
            y_t = plant_t.step(u_t);
            y_u = plant_u.step(u_u);
            if (k == 599) {
                CHECK(told.delayed_model_output() == doctest::Approx(rail * K).epsilon(1e-3)); // the plant's own output
                CHECK(untold.delayed_model_output() > 10.0 * rail);                            // an output the plant never had
            }
        }
        CHECK(cut > 100);
        CHECK(y_t == doctest::Approx(0.5).epsilon(0.02)); // back on the reference
        CHECK(y_u > 0.75);                                // still unwinding, on the rail
    }

    TEST_CASE("track: the models run on the input applied, and the next step at the same r and y returns it") {
        //! It had no handover: a predictor switched in started from its last commands in flight and
        //! its primary's last integral.
        constexpr double K = 1.0;
        constexpr double tau = 1.0;
        constexpr double L = 0.5;
        constexpr double Ts = 0.05;
        constexpr auto   sp = design::smith_predictor(design::FOPDT{.K = K, .tau = tau}, Ts, InputDelay<16>{L}, design::SIMC{.tau_c = 0.5});
        REQUIRE(sp.success);
        SmithPredictor smith{sp};
        // The plant K/(τs + 1) with its delay, as discretize samples it, read before each period's input.
        const StateSpace<1, 1, 1> plant{.A = {{-1.0 / tau}}, .B = {{K / tau}}, .C = {{1.0}}};
        const auto                delayed = discretize(plant, Ts, InputDelay<16>{L});
        REQUIRE(delayed.success);
        ColVec<17, double> x{};
        const double       u_op = 0.6; // an operator drives the plant
        for (int k = 0; k < 200; ++k) {
            smith.track(u_op, 1.0, (delayed->C * x)(0, 0));
            x = (delayed->A * x) + (delayed->B * ColVec<1, double>{u_op});
        }
        const double y = (delayed->C * x)(0, 0);
        CHECK(smith.last_command() == u_op);
        CHECK(smith.delayed_model_output() == doctest::Approx(y).epsilon(1e-12)); // the model followed the plant
        smith.track(u_op, 1.0, y);
        CHECK(smith.step(1.0, y) == doctest::Approx(u_op).epsilon(1e-12));
    }

    TEST_CASE("Law of a Smith predictor design: N deduced, a failed design refused where it is declared") {
        static constexpr auto kSmith = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 2.0}, 0.1, InputDelay<12>{1.0}, design::SIMC{.tau_c = 1.0});
        Law<kSmith, float>    law;
        static_assert(std::is_same_v<Law<kSmith, float>::Runtime, SmithPredictor<12, float>>);
        SmithPredictor<12, float> held{kSmith};
        for (int k = 0; k < 30; ++k) {
            CHECK(law.step(1.0f, 0.0f) == held.step(1.0f, 0.0f));
        }
    }
}
