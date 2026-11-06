// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_smith_predictor.cpp
 * @brief Smith predictor: design, L=0 equivalence to PID, dead-time closed loop,
 *        custom SisoController primary
 */

#include <cmath>
#include <cstddef>
#include <limits>

#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/controllers/smith_predictor.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/filters/delay.hpp"
#include "damp/systems/transfer_function.hpp"

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

/// Minimal custom primary to prove Primary is not PID-only
struct ProportionalPrimary {
    double k{1.0};

    [[nodiscard]] constexpr double control(double r, double y) {
        return k * (r - y);
    }

    constexpr void reset() {}
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

    double control(double r, double y) {
        const double y_hat = x;
        const double y_delayed = delay(y_hat);
        const double u = primary.control(r, y + y_hat - y_delayed);
        x = (a * x) + (b * u);
        return u;
    }
};

/// Largest |u_new − u_old| over @p N ticks of the rebuilt and the legacy law, each in closed loop
/// with its own copy of the FOPDT plant (K, τ, L at Ts), reference steps at 0 and N/2.
template<size_t MaxD, typename Primary>
double max_command_gap(const Primary& c, double K, double tau, double L, double Ts, size_t N) {
    const double                          a = std::exp(-Ts / tau);
    const double                          b = K * (1.0 - a);
    const auto                            d = static_cast<size_t>((L / Ts) + 0.5);
    SmithPredictor<MaxD, double, Primary> rebuilt(c, a, b, d);
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
        const double u_new = rebuilt.control(r, y_new);
        const double u_old = legacy.control(r, y_old);
        gap = std::max(gap, std::abs(u_new - u_old));
        y_new = plant_new.step(u_new);
        y_old = plant_old.step(u_old);
    }
    return gap;
}

static_assert(SisoController<ProportionalPrimary, double>);
static_assert(SisoController<SmithPredictor<8, double>, double>);
static_assert(SisoController<SmithPredictor<8, float, PIController<float>>, float>);

} // namespace

TEST_SUITE("Smith Predictor") {

    TEST_CASE("design::smith_predictor_from_fopdt packs FO model and delay") {
        constexpr double K = 1.0;
        constexpr double tau = 5.0;
        constexpr double L = 2.0;
        constexpr double Ts = 0.1;
        constexpr double tau_c = 2.0;

        constexpr auto sp = design::smith_predictor_from_fopdt(K, tau, L, Ts, tau_c);
        static_assert(sp.success);
        static_assert(sp.delay_samples == 20); // round(2.0/0.1)

        const double a_ref = std::exp(-Ts / tau);
        const double b_ref = K * (1.0 - a_ref);
        CHECK(sp.a == doctest::Approx(a_ref).epsilon(1e-12));
        CHECK(sp.b == doctest::Approx(b_ref).epsilon(1e-12));
        CHECK(sp.pid.Ts == doctest::Approx(Ts));
        CHECK(sp.pid.Kp > 0.0);
        CHECK(sp.pid.Ki > 0.0); // PI from SIMC

        // SIMC on delay-free plant: Kp = tau / (K * tau_c)
        const auto cont = design::simc(K, 0.0, tau, tau_c, PIDMode::PI);
        CHECK(sp.pid.Kp == doctest::Approx(cont.Kp));
        CHECK(sp.pid.Ki == doctest::Approx(cont.Ki * Ts));
    }

    TEST_CASE("design rejects non-physical FOPDT") {
        CHECK_FALSE(design::smith_predictor_from_fopdt(0.0, 5.0, 1.0, 0.1).success);
        CHECK_FALSE(design::smith_predictor_from_fopdt(1.0, 0.0, 1.0, 0.1).success);
        CHECK_FALSE(design::smith_predictor_from_fopdt(1.0, 5.0, 1.0, 0.0).success);
        CHECK_FALSE(design::smith_predictor_from_fopdt(1.0, 5.0, -1.0, 0.1).success);
    }

    TEST_CASE("as<float> preserves success and delay samples") {
        constexpr auto sp = design::smith_predictor_from_fopdt(1.0, 3.0, 1.5, 0.05, 1.0);
        static_assert(sp.success);
        constexpr auto spf = sp.as<float>();
        static_assert(spf.success);
        static_assert(spf.delay_samples == sp.delay_samples);
        CHECK(static_cast<double>(spf.a) == doctest::Approx(sp.a).epsilon(1e-6));
    }

    TEST_CASE("smith_predictor_from_pid matches manual discretize path") {
        constexpr auto pid = design::pid(2.0, 0.4, 0.0);
        constexpr auto sp = design::smith_predictor_from_pid(pid, 1.0, 4.0, 1.0, 0.1);
        static_assert(sp.success);
        CHECK(sp.pid.Kp == doctest::Approx(2.0));
        CHECK(sp.pid.Ki == doctest::Approx(0.4 * 0.1));
        CHECK(sp.delay_samples == 10);
    }

    TEST_CASE("L=0 Smith predictor matches bare PI on FO plant") {
        constexpr double K = 1.0;
        constexpr double tau = 2.0;
        constexpr double L = 0.0;
        constexpr double Ts = 0.05;
        constexpr double tau_c = 1.0;

        constexpr auto sp = design::smith_predictor_from_fopdt(K, tau, L, Ts, tau_c);
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
            const double u = smith.control(r, y);
            const double u_b = bare.control(r, y_bare);
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

        constexpr auto sp = design::smith_predictor_from_fopdt(K, tau, L, Ts, tau_c);
        REQUIRE(sp.success);
        REQUIRE(sp.delay_samples == 60);

        SmithPredictor<128, double> smith(sp);
        PIController<double>        bare(sp.pid);

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
            const double u_s = smith.control(r, y_s);
            const double u_b = bare.control(r, y_b);
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

    TEST_CASE("reset clears model and delay so next step is cold-start") {
        constexpr auto             sp = design::smith_predictor_from_fopdt(1.0, 2.0, 0.5, 0.1, 1.0);
        SmithPredictor<16, double> ctrl(sp);
        REQUIRE(ctrl.valid());

        double          y = 0.0;
        FirstOrderPlant plant{.a = sp.a, .b = sp.b};
        // Drive some state
        for (int i = 0; i < 30; ++i) {
            const double u = ctrl.control(1.0, y);
            y = plant.step(u);
        }
        CHECK(std::abs(ctrl.model_output()) > 1e-6);

        ctrl.reset();
        CHECK(ctrl.model_output() == doctest::Approx(0.0));
        // First post-reset command with y=0,r=0 should be ~0 (no integral yet)
        CHECK(ctrl.control(0.0, 0.0) == doctest::Approx(0.0));
    }

    TEST_CASE("invalid design or MaxDelay overflow yields zero control") {
        SmithPredictor<8, double> bad(design::smith_predictor_from_fopdt(0.0, 1.0, 0.0, 0.1));
        CHECK_FALSE(bad.valid());
        CHECK(bad.control(1.0, 0.0) == doctest::Approx(0.0));

        // delay_samples = 50 > MaxDelay=8
        constexpr auto long_delay = design::smith_predictor_from_fopdt(1.0, 1.0, 5.0, 0.1, 1.0);
        static_assert(long_delay.delay_samples == 50);
        SmithPredictor<8, double> too_short(long_delay);
        CHECK_FALSE(too_short.valid());
        CHECK(too_short.control(1.0, 0.0) == doctest::Approx(0.0));

        // MaxDelay is a capacity in whole samples: d = MaxDelay fits, d = MaxDelay + 1 does not.
        constexpr auto d_eq = design::smith_predictor_from_fopdt(1.0, 1.0, 0.8, 0.1, 1.0);
        static_assert(d_eq.delay_samples == 8);
        SmithPredictor<8, double> edge(d_eq);
        CHECK(edge.valid());
        constexpr auto d_over = design::smith_predictor_from_fopdt(1.0, 1.0, 0.9, 0.1, 1.0);
        static_assert(d_over.delay_samples == 9);
        SmithPredictor<8, double> over(d_over);
        CHECK_FALSE(over.valid());
    }

    TEST_CASE("the law on its ring of commands in flight reproduces the Delay-line law's commands") {
        // The same plants as the cases above. ŷ is now the delayed model carried over the d
        // commands in flight instead of a stored delay-free model, so the two agree to the
        // rounding of d multiply-adds on commands of order 1–10: 1e-12 absolute.
        const auto pi = [](double K, double tau, double tau_c, double Ts) {
            return PIController<double>(design::simc(K, 0.0, tau, tau_c, PIDMode::PI).discretize(Ts));
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

    TEST_CASE("custom SisoController primary (not PID-family)") {
        constexpr double K = 1.0;
        constexpr double tau = 2.0;
        constexpr double Ts = 0.05;
        const double     a = std::exp(-Ts / tau);
        const double     b = K * (1.0 - a);

        // L = 0: y_pred = y, so Smith + pure P matches bare P on the FO plant
        ProportionalPrimary                            primary{.k = 1.5};
        SmithPredictor<8, double, ProportionalPrimary> smith(primary, a, b, 0);
        REQUIRE(smith.valid());

        ProportionalPrimary bare{.k = 1.5};
        FirstOrderPlant     plant{.a = a, .b = b};
        FirstOrderPlant     plant_bare{.a = a, .b = b};

        constexpr double r = 1.0;
        double           y = 0.0;
        double           y_bare = 0.0;
        double           max_diff = 0.0;
        for (size_t k = 0; k < 200; ++k) {
            const double u = smith.control(r, y);
            const double u_b = bare.control(r, y_bare);
            max_diff = std::max(max_diff, std::abs(u - u_b));
            y = plant.step(u);
            y_bare = plant_bare.step(u_b);
        }
        // L=0 residual is identity: Smith+P and bare P must match sample-for-sample.
        // Pure P on a type-0 plant does not eliminate steady-state error — only compare u.
        CHECK(max_diff < 1e-12);
        CHECK(smith.primary().k == doctest::Approx(1.5));
        (void)y;
        (void)r;
    }

    TEST_CASE("full PID primary via template parameter") {
        // Smith SIMC hides L from C, so the packed primary is a PI; the point of this
        // case is that Primary can be PIDController (full mode), not only PI.
        constexpr auto sp = design::smith_predictor_from_fopdt(1.0, 2.0, 0.0, 0.05, 1.0);
        REQUIRE(sp.success);
        SmithPredictor<8, double, PIDController<double>> smith(sp);
        REQUIRE(smith.valid());
        CHECK(std::isfinite(smith.control(1.0, 0.0)));
    }

    TEST_CASE("as<float> converting ctor preserves model coeffs") {
        constexpr auto             sp = design::smith_predictor_from_fopdt(1.0, 3.0, 0.5, 0.1, 1.0);
        SmithPredictor<16, double> src(sp);
        REQUIRE(src.valid());
        SmithPredictor<16, float> dst(src);
        REQUIRE(dst.valid());
        CHECK(static_cast<double>(dst.model_output()) == doctest::Approx(0.0));
        CHECK(dst.delay_samples() == src.delay_samples());
    }

    TEST_CASE("constexpr constructibility of design result") {
        constexpr auto sp = design::smith_predictor_from_fopdt(1.0, 1.0, 0.2, 0.05, 0.5);
        static_assert(sp.success);
        static_assert(sp.delay_samples == 4);
        constexpr auto disc = sp.pid;
        static_assert(disc.Ts > 0.0);
        (void)disc;
        // The runtime ticks under constant evaluation too (the commands in flight are a fixed array).
        constexpr double u2 = [] {
            SmithPredictor<8, double> c(design::smith_predictor_from_fopdt(1.0, 1.0, 0.2, 0.05, 0.5));
            (void)c.control(1.0, 0.0);
            return c.control(1.0, 0.0);
        }();
        static_assert(u2 > 0.0);
    }

    TEST_CASE("to_tf recovers SIMC primary; plant_tf is K/(τs+1)") {
        constexpr auto sp = design::smith_predictor_from_fopdt(2.0, 5.0, 1.0, 0.1, 2.0);
        REQUIRE(sp.success);
        const auto C = design::simc(2.0, 0.0, 5.0, 2.0, PIDMode::PI);
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

    TEST_CASE("smith designs carry a failed primary's status and refuse an unrepresentable delay") {
        // A Kbc against Ki's sign is refused by discretize; the Smith design used to report success.
        auto wrong_kbc = design::pid(2.0, 0.4, 0.0);
        wrong_kbc.Kbc = -1.0;
        const auto from_bad_kbc = design::smith_predictor_from_pid(wrong_kbc, 1.0, 4.0, 1.0, 0.1);
        CHECK_FALSE(from_bad_kbc.success);
        CHECK(from_bad_kbc.status == design::DesignStatus::invalid_parameter);
        SmithPredictor<16, double> runtime(from_bad_kbc);
        CHECK_FALSE(runtime.valid());

        const auto failed_tuning = design::simc(0.0, 0.0, 4.0, 1.0); // K = 0: refused
        REQUIRE_FALSE(failed_tuning.success);
        CHECK(design::smith_predictor_from_pid(failed_tuning, 1.0, 4.0, 1.0, 0.1).status == failed_tuning.status);

        // round(L/Ts) has no size_t value for an infinite or huge L.
        const double inf = std::numeric_limits<double>::infinity();
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const auto   pid = design::pid(2.0, 0.4, 0.0);
        CHECK(design::smith_predictor_from_pid(pid, 1.0, 4.0, inf, 0.1).status == design::DesignStatus::invalid_parameter);
        CHECK(design::smith_predictor_from_pid(pid, 1.0, 4.0, 1e300, 0.1).status == design::DesignStatus::invalid_parameter);
        CHECK(design::smith_predictor_from_fopdt(1.0, 4.0, inf, 0.1).status == design::DesignStatus::invalid_parameter);
        CHECK(design::smith_predictor_from_fopdt(1.0, nan, 1.0, 0.1).status == design::DesignStatus::invalid_parameter);
        CHECK(design::smith_predictor_from_fopdt(1.0, 4.0, 1.0, nan).status == design::DesignStatus::invalid_sample_time);
    }
}
