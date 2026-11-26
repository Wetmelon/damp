// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <limits>
#include <numbers>
#include <type_traits>
#include <utility>

#include "damp/backend.hpp"
#include "damp/controllers/lead_lag.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for Lead-Lag compensator (lead_lag.hpp)
 *
 * Validates lead/lag design, discretization, and step response.
 * Reference values generated with scipy.signal.
 */

namespace {
constexpr double kTs = 1e-3; // every design is deployed at a Ts; the continuous maps do not depend on it
} // namespace

TEST_SUITE("Lead-Lag Compensator") {

    TEST_CASE("Lead design - 45 degree phase boost at wc=1000") {
        constexpr double phi = std::numbers::pi / 4.0; // 45 degrees
        constexpr double wc = 1000.0;
        constexpr auto   result = design::lead({.phase_max = phi, .wc = wc}, kTs);

        // Check zero, pole, gain against analytical values
        // alpha = (1-sin(45))/(1+sin(45)) = 0.17157...
        // sqrt(alpha) = 0.41421...
        static_assert(result.gains[0].z < result.gains[0].p); // Lead: zero < pole
        CHECK(result.gains[0].z == doctest::Approx(414.2135623730950).epsilon(1e-12));
        CHECK(result.gains[0].p == doctest::Approx(2414.213562373095).epsilon(1e-12));
        CHECK(result.gains[0].K == doctest::Approx(2.414213562373095).epsilon(1e-12));
    }

    TEST_CASE("Lead design - phase at wc is exactly phi_max") {
        // C(jw) = K*(jw + z)/(jw + p)
        // Phase at wc should be phi_max, magnitude should be 1.0
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr auto   r = design::lead({.phase_max = phi, .wc = wc}, kTs).gains[0];

        // Phase = atan(wc/z) - atan(wc/p)
        double phase_z = std::atan(wc / r.z);
        double phase_p = std::atan(wc / r.p);
        double phase_deg = (phase_z - phase_p) * 180.0 / std::numbers::pi;
        CHECK(phase_deg == doctest::Approx(45.0).epsilon(1e-10));

        // Magnitude = K * sqrt(wc²+z²) / sqrt(wc²+p²)
        double mag = r.K * std::sqrt((wc * wc) + (r.z * r.z)) / std::sqrt((wc * wc) + (r.p * r.p));
        CHECK(mag == doctest::Approx(1.0).epsilon(1e-12));
    }

    TEST_CASE("Lead design - 30 degree boost at wc=500") {
        constexpr double phi = 30.0 * std::numbers::pi / 180.0;
        constexpr double wc = 500.0;
        constexpr auto   r = design::lead({.phase_max = phi, .wc = wc}, kTs).gains[0];

        CHECK(r.z == doctest::Approx(288.6751345948128).epsilon(1e-12));
        CHECK(r.p == doctest::Approx(866.0254037844387).epsilon(1e-12));
        CHECK(r.K == doctest::Approx(1.732050807568877).epsilon(1e-12));

        // Verify phase and magnitude at wc
        double phase_z = std::atan(wc / r.z);
        double phase_p = std::atan(wc / r.p);
        double phase_deg = (phase_z - phase_p) * 180.0 / std::numbers::pi;
        CHECK(phase_deg == doctest::Approx(30.0).epsilon(1e-10));

        double mag = r.K * std::sqrt((wc * wc) + (r.z * r.z)) / std::sqrt((wc * wc) + (r.p * r.p));
        CHECK(mag == doctest::Approx(1.0).epsilon(1e-12));
    }

    TEST_CASE("Lag design - 10x DC gain boost") {
        constexpr double dc_boost = 10.0;
        constexpr double wc = 1000.0;
        constexpr auto   r = design::lag({.dc_gain_boost = dc_boost, .wc = wc}, kTs).gains[0];

        // z = wc/margin_factor = 1000/10 = 100
        // p = z/dc_gain_boost = 100/10 = 10
        static_assert(r.z > r.p); // Lag: zero > pole
        CHECK(r.z == doctest::Approx(100.0).epsilon(1e-12));
        CHECK(r.p == doctest::Approx(10.0).epsilon(1e-12));
        CHECK(r.K == doctest::Approx(1.0).epsilon(1e-12));

        // DC gain = K * z / p = 10
        double dc_gain = r.K * r.z / r.p;
        CHECK(dc_gain == doctest::Approx(10.0).epsilon(1e-12));

        // Phase at wc should be small (< 6 degrees lag)
        double phase_z = std::atan(wc / r.z);
        double phase_p = std::atan(wc / r.p);
        double phase_deg = (phase_z - phase_p) * 180.0 / std::numbers::pi;
        CHECK(phase_deg == doctest::Approx(-5.137654439816157).epsilon(1e-6));
    }

    TEST_CASE("Lead-lag direct specification") {
        constexpr auto r = design::lead_lag_direct({.K = 2.0, .z = 100.0, .p = 500.0}, kTs);
        static_assert(r.gains[0].K == 2.0);
        static_assert(r.gains[0].z == 100.0);
        static_assert(r.gains[0].p == 500.0);
        static_assert(r.Ts == kTs);
    }

    TEST_CASE("LeadLagResult to_tf") {
        constexpr auto r = design::lead_lag_direct({.K = 2.0, .z = 100.0, .p = 500.0}, kTs);
        constexpr auto tf = r.to_tf();

        // C(s) = 2*(s+100)/(s+500) = (2s + 200)/(s + 500)
        CHECK(tf.num[0] == doctest::Approx(2.0).epsilon(1e-12));   // K
        CHECK(tf.num[1] == doctest::Approx(200.0).epsilon(1e-12)); // K*z
        CHECK(tf.den[0] == doctest::Approx(1.0).epsilon(1e-12));   // 1
        CHECK(tf.den[1] == doctest::Approx(500.0).epsilon(1e-12)); // p
    }

    TEST_CASE("LeadLagResult to_ss") {
        constexpr auto r = design::lead_lag_direct({.K = 2.0, .z = 100.0, .p = 500.0}, kTs);
        constexpr auto ss = r.to_ss();

        // A = [-p] = [-500], B = [1], C = [K*(z-p)] = [2*(100-500)] = [-800], D = [K] = [2]
        CHECK(ss.A(0, 0) == doctest::Approx(-500.0).epsilon(1e-12));
        CHECK(ss.B(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(ss.C(0, 0) == doctest::Approx(-800.0).epsilon(1e-12));
        CHECK(ss.D(0, 0) == doctest::Approx(2.0).epsilon(1e-12));
    }

    TEST_CASE("Tustin discretization matches scipy") {
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr double Ts = 0.001;

        constexpr auto            r = design::lead({.phase_max = phi, .wc = wc}, Ts);
        LeadLagController<double> ctrl(r);

        // Coefficients verified against scipy.signal.cont2discrete bilinear, read off the impulse
        // response u = b0·e + b1·e_prev − a1·u_prev: u0 = b0, u1 = b1 − a1·b0, u2 = −a1·u1.
        const double u0 = ctrl.step(1.0);
        const double u1 = ctrl.step(0.0);
        const double u2 = ctrl.step(0.0);
        const double a1 = -u2 / u1;
        CHECK(u0 == doctest::Approx(1.320377241017041).epsilon(1e-12));
        CHECK(u1 + (a1 * u0) == doctest::Approx(-0.8672954016950677).epsilon(1e-12));
        CHECK(a1 == doctest::Approx(0.09383632135605431).epsilon(1e-12));
    }

    TEST_CASE("Discrete step response matches scipy") {
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr double Ts = 0.001;

        constexpr auto            r = design::lead({.phase_max = phi, .wc = wc}, Ts);
        LeadLagController<double> ctrl(r);

        // Step response reference from scipy simulation
        constexpr double expected[] = {
            1.320377241017041e+00,
            3.291824962226776e-01,
            4.221925648216337e-01,
            4.134648421352334e-01,
            4.142838195259410e-01,
        };

        for (int i = 0; i < 5; ++i) {
            double y = ctrl.step(1.0);
            CHECK(y == doctest::Approx(expected[i]).epsilon(1e-12));
        }
    }

    TEST_CASE("DC gain consistency") {
        // Lead: DC gain = K*z/p < 1 (attenuates DC)
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr auto   lead_r = design::lead({.phase_max = phi, .wc = wc}, kTs).gains[0];
        double           dc_lead = lead_r.K * lead_r.z / lead_r.p;
        CHECK(dc_lead == doctest::Approx(0.4142135623730950).epsilon(1e-12));
        CHECK(dc_lead < 1.0); // Lead attenuates DC

        // Lag: DC gain = K*z/p > 1 (boosts DC)
        constexpr auto lag_r = design::lag({.dc_gain_boost = 10.0, .wc = wc}, kTs).gains[0];
        double         dc_lag = lag_r.K * lag_r.z / lag_r.p;
        CHECK(dc_lag == doctest::Approx(10.0).epsilon(1e-12));
        CHECK(dc_lag > 1.0); // Lag boosts DC
    }

    TEST_CASE("Controller reset clears state") {
        constexpr double          phi = std::numbers::pi / 4.0;
        constexpr auto            r = design::lead({.phase_max = phi, .wc = 1000.0}, 0.001);
        LeadLagController<double> ctrl(r);

        // Run some steps
        (void)ctrl.step(1.0);
        (void)ctrl.step(1.0);
        (void)ctrl.step(1.0);

        ctrl.reset();

        // After reset, first step should be same as fresh controller
        LeadLagController<double> fresh(r);
        CHECK(ctrl.step(1.0) == doctest::Approx(fresh.step(1.0)).epsilon(1e-15));
    }

    TEST_CASE("Online lead design matches design") {
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr auto   design_r = design::lead({.phase_max = phi, .wc = wc}, kTs);
        auto             online_r = design::lead({.phase_max = phi, .wc = wc}, kTs);

        CHECK(online_r.gains[0].K == doctest::Approx(design_r.gains[0].K).epsilon(1e-12));
        CHECK(online_r.gains[0].z == doctest::Approx(design_r.gains[0].z).epsilon(1e-12));
        CHECK(online_r.gains[0].p == doctest::Approx(design_r.gains[0].p).epsilon(1e-12));
    }

    TEST_CASE("Online lag design matches design") {
        constexpr double wc = 1000.0;
        constexpr auto   design_r = design::lag({.dc_gain_boost = 10.0, .wc = wc}, kTs);
        auto             online_r = design::lag({.dc_gain_boost = 10.0, .wc = wc}, kTs);

        CHECK(online_r.gains[0].K == doctest::Approx(design_r.gains[0].K).epsilon(1e-12));
        CHECK(online_r.gains[0].z == doctest::Approx(design_r.gains[0].z).epsilon(1e-12));
        CHECK(online_r.gains[0].p == doctest::Approx(design_r.gains[0].p).epsilon(1e-12));
    }

    TEST_CASE("Lead-lag combined (2nd order)") {
        // 45-degree lead + 10x lag at wc=1000
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr auto   result = design::lead_lag({.phase_max = phi, .wc = wc, .dc_gain_boost = 10.0}, kTs);
        static_assert(result.success);
        constexpr auto ss = result.to_ss(); // the continuous series

        // Should be a 2x2 state-space system (series of two 1st-order)
        // Check dimensions compile
        static_assert(ss.A.rows() == 2);
        static_assert(ss.A.cols() == 2);
        static_assert(ss.B.rows() == 2);
        static_assert(ss.B.cols() == 1);
        static_assert(ss.C.rows() == 1);
        static_assert(ss.C.cols() == 2);

        // DC gain of combined = lead_dc * lag_dc = 0.41421... * 10 = 4.1421...
        // Compute via state-space: -C * A^-1 * B + D
        auto A_inv = ss.A.inverse();
        REQUIRE(A_inv.has_value());
        auto dc = (ss.C * (-A_inv.value()) * ss.B + ss.D);
        CHECK(dc(0, 0) == doctest::Approx(0.4142135623730950 * 10.0).epsilon(1e-10));
        // And from the transfer function: K·z_lead·z_lag/(p_lead·p_lag).
        const auto tf = result.to_tf();
        CHECK(tf.num[2] / tf.den[2] == doctest::Approx(0.4142135623730950 * 10.0).epsilon(1e-10));
    }

    TEST_CASE("Lead-lag series propagates section failure") {
        // Invalid lead phi -> success false, zero system
        static_assert(!design::lead_lag({.phase_max = 0.0, .wc = 1000.0, .dc_gain_boost = 10.0}, kTs).success);
        static_assert(!design::lead_lag({.phase_max = std::numbers::pi / 4.0, .wc = 1000.0, .dc_gain_boost = 1.0}, kTs).success); // lag needs boost > 1
        constexpr auto bad = design::lead_lag({.phase_max = -1.0, .wc = 1000.0, .dc_gain_boost = 10.0}, kTs);
        static_assert(!bad.success);
        CHECK(bad.to_ss().A(0, 0) == doctest::Approx(0.0));
        CHECK(bad.to_ss().D(0, 0) == doctest::Approx(0.0));
    }

    TEST_CASE("Type conversion via .as<float>()") {
        constexpr auto r = design::lead({.phase_max = std::numbers::pi / 4.0, .wc = 1000.0}, 0.001);
        constexpr auto rf = r.as<float>();

        LeadLagController<float> ctrl(rf);
        float                    y = ctrl.step(1.0f);
        CHECK(y == doctest::Approx(1.320377f).epsilon(1e-4));
    }

    TEST_CASE("control(r, y) equals control(r - y)") {
        constexpr double          phi = std::numbers::pi / 4.0;
        constexpr auto            r = design::lead({.phase_max = phi, .wc = 1000.0}, 0.001);
        LeadLagController<double> a(r);
        LeadLagController<double> b(r);

        // Feed the same error through both forms; outputs must match step for step.
        const double refs[] = {1.0, 1.0, 0.5, -0.3, 0.0};
        const double meas[] = {0.0, 0.4, 0.6, 0.1, -0.2};
        for (int i = 0; i < 5; ++i) {
            double ya = a.step(refs[i], meas[i]);
            double yb = b.step(refs[i] - meas[i]);
            CHECK(ya == doctest::Approx(yb).epsilon(1e-15));
        }
    }

    TEST_CASE("Invalid specs report success = false") {
        // Lead requires 0 < phi < pi/2 and wc > 0.
        static_assert(!design::lead({.phase_max = 0.0, .wc = 1000.0}, kTs).success);
        static_assert(!design::lead({.phase_max = std::numbers::pi / 2.0, .wc = 1000.0}, kTs).success);
        static_assert(!design::lead({.phase_max = std::numbers::pi / 4.0, .wc = -1.0}, kTs).success);
        static_assert(design::lead({.phase_max = std::numbers::pi / 4.0, .wc = 1000.0}, kTs).success);

        // Lag requires dc_gain_boost > 1, wc > 0, margin > 0.
        static_assert(!design::lag({.dc_gain_boost = 1.0, .wc = 1000.0}, kTs).success);
        static_assert(!design::lag({.dc_gain_boost = 10.0, .wc = 0.0}, kTs).success);
        static_assert(design::lag({.dc_gain_boost = 10.0, .wc = 1000.0}, kTs).success);

        static_assert(!design::lead_lag_direct({.K = 1.0, .z = -1.0, .p = 5.0}, kTs).success);
        static_assert(design::lead_lag_direct({.K = 1.0, .z = 100.0, .p = 500.0}, kTs).success);
    }

    TEST_CASE("Failed design constructs fail-closed controller (emits zero)") {
        constexpr auto bad = design::lead({.phase_max = 0.0, .wc = 1000.0}, 0.001);
        static_assert(!bad.success);
        LeadLagController<double> ctrl(bad);
        CHECK_FALSE(ctrl.valid());
        CHECK(ctrl.step(1.0) == doctest::Approx(0.0));
        CHECK(ctrl.step(2.0, 1.0) == doctest::Approx(0.0));
        ctrl.reset();
        CHECK(ctrl.step(-3.0) == doctest::Approx(0.0));
    }

    TEST_CASE("to_discrete_ss roundtrip (dimensions and sample time)") {
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr double Ts = 0.001;
        constexpr auto   r = design::lead({.phase_max = phi, .wc = wc}, Ts);

        // Get discrete SS
        constexpr auto dss = r.to_discrete_ss();

        // Verify dimensions
        static_assert(dss.A.rows() == 1);
        static_assert(dss.A.cols() == 1);
        static_assert(dss.B.rows() == 1);
        static_assert(dss.B.cols() == 1);
        static_assert(dss.C.rows() == 1);
        static_assert(dss.C.cols() == 1);

        // Verify sampling time set
        CHECK(dss.Ts == doctest::Approx(Ts).epsilon(1e-15));
    }

    TEST_CASE("to_discrete_ss matches controller step response") {
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr double Ts = 0.001;
        constexpr auto   r = design::lead({.phase_max = phi, .wc = wc}, Ts);

        // Discrete SS via Tustin
        constexpr auto dss = r.to_discrete_ss();

        // Verify SS matrices against scipy (bilinear)
        CHECK(dss.A(0, 0) == doctest::Approx(-9.383632135605437e-02).epsilon(1e-12));
        CHECK(dss.B(0, 0) == doctest::Approx(4.530818393219728e-04).epsilon(1e-12));
        CHECK(dss.C(0, 0) == doctest::Approx(-2.187672642712108e+03).epsilon(1e-6));
        CHECK(dss.D(0, 0) == doctest::Approx(1.320377241017041e+00).epsilon(1e-12));

        // Step response from state-space should match IIR controller
        LeadLagController<double> ctrl(r);

        Matrix<1, 1, double> x = Matrix<1, 1, double>::zeros();
        Matrix               u{{1.0}};

        for (int i = 0; i < 5; ++i) {
            auto   y_ss = dss.C * x + dss.D * u;
            double y_ctrl = ctrl.step(1.0);
            CHECK(y_ss(0, 0) == doctest::Approx(y_ctrl).epsilon(1e-10));
            x = dss.A * x + dss.B * u;
        }
    }

    TEST_CASE("Ts is refused unless finite and positive; NaN specs are refused") {
        // Ts = 0 was a continuous design, which the runtime refused (it used to emit the bare gain K).
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (const double Ts : {0.0, -1e-3, std::numeric_limits<double>::infinity(), nan}) {
            CAPTURE(Ts);
            const auto lead_r = design::lead({.phase_max = 0.5, .wc = 100.0}, Ts);
            CHECK(lead_r.status == design::DesignStatus::invalid_sample_time);
            LeadLagController<double> lead_c{lead_r};
            CHECK_FALSE(lead_c.valid());
            CHECK(lead_c.step(1.0) == 0.0);
            CHECK(design::lag({.dc_gain_boost = 10.0, .wc = 100.0}, Ts).status == design::DesignStatus::invalid_sample_time);
            CHECK(design::lead_lag({.phase_max = 0.5, .wc = 100.0, .dc_gain_boost = 10.0}, Ts).status == design::DesignStatus::invalid_sample_time);
            CHECK(design::lead_lag_direct({.K = 1.0, .z = 5.0, .p = 50.0}, Ts).status == design::DesignStatus::invalid_sample_time);
        }
        auto edited = design::lag({.dc_gain_boost = 10.0, .wc = 100.0}, 1e-3);
        edited.Ts = 0.0;
        CHECK_FALSE(LeadLagController<double>{edited}.valid());

        // Ts is the second argument of every member of the family.
        constexpr auto lag_d = design::lag({.dc_gain_boost = 10.0, .wc = 100.0}, 1e-3);
        static_assert(lag_d.success && lag_d.Ts == 1e-3 && lag_d.gains[0].z == 10.0);

        CHECK(design::lead({.phase_max = nan, .wc = 100.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lead({.phase_max = 0.5, .wc = nan}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lag({.dc_gain_boost = nan, .wc = 100.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lead_lag_direct({.K = 1.0, .z = nan, .p = 5.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE_TEMPLATE("InputLimits: a lag section holds against the clamp, its own or one reported", T, float, double) {
        //! The recursion ran on its own unclamped output: behind a clamp the slow pole kept charging
        //! towards the DC gain, and the output stayed on the rail long after the error reversed.
        constexpr double Ts = 1e-3;
        const auto       res = design::lag({.dc_gain_boost = 10.0, .wc = 20.0}, Ts, design::InputLimits{.lower = {-1.0}, .upper = {1.0}}); // zero at 2 rad/s, pole at 0.2: ten times the gain at DC
        const auto       open = design::lag({.dc_gain_boost = 10.0, .wc = 20.0}, Ts);
        REQUIRE(res.success);
        LeadLagController<T> bounded{res.template as<T>()};
        LeadLagController<T> outside{open.template as<T>()};
        LeadLagController<T> untold{open.template as<T>()};
        int                  bounded_on_rail = 0;
        int                  untold_on_rail = 0;
        for (int k = 0; k < 6000; ++k) {
            const T e = (k < 3000) ? T{1} : static_cast<T>(-0.02);
            const T ua = bounded.step(e);
            const T ub = damp::clamp(outside.step(e), T{-1}, T{1});
            outside.back_calculate(ub);
            const T uc = damp::clamp(untold.step(e), T{-1}, T{1});
            CHECK(damp::abs(ua - ub) <= T{64} * eps<T>()); // one tick, whichever stage clamps
            if (k < 3000) {
                CHECK(ua == T{1});
            } else {
                bounded_on_rail += (ua >= T{1}) ? 1 : 0;
                untold_on_rail += (uc >= T{1}) ? 1 : 0;
            }
        }
        CHECK(bounded_on_rail == 0);   // off the rail on the tick the error reverses
        CHECK(untold_on_rail == 3000); // the unreported clamp: still charged three seconds later
        // The command, not the applied output: back_calculate moves the recursion, not last_command.
        const T u = outside.step(T{1});
        outside.back_calculate(T{0});
        CHECK(outside.last_command() == u);
    }

    TEST_CASE("a limit that is not a box is refused: by the design, and by a runtime built from an edited result") {
        CHECK(design::lag({.dc_gain_boost = 10.0, .wc = 20.0}, 1e-3, design::InputLimits{.lower = {1.0}, .upper = {-1.0}}).status == design::DesignStatus::invalid_parameter);
        auto res = design::lag({.dc_gain_boost = 10.0, .wc = 20.0}, 1e-3);
        res.limits.input = {.lower = 1.0, .upper = -1.0};
        LeadLagController<double> inert{res};
        CHECK_FALSE(inert.valid());
        CHECK(inert.step(1.0) == 0.0);
    }
}

namespace {
constexpr auto kLeadLag = design::lead_lag({.phase_max = 0.6, .wc = 20.0, .dc_gain_boost = 10.0}, 1e-3, design::InputLimits{.lower = {-2.0}, .upper = {2.0}});
} // namespace

TEST_SUITE("Lead-Lag Compensator: the runtime contract") {
    TEST_CASE("the lead-lag series runs: LeadLagController<T, 2> is its Tustin model tick for tick") {
        // design::lead_lag returned a StateSpace that no runtime ran.
        static_assert(kLeadLag.success);
        const auto                   dss = kLeadLag.to_discrete_ss();
        const auto                   open = design::lead_lag({.phase_max = 0.6, .wc = 20.0, .dc_gain_boost = 10.0}, 1e-3);
        LeadLagController<double, 2> ctrl{open};
        Matrix<2, 1, double>         x{};
        for (int k = 0; k < 300; ++k) {
            const double e = std::sin(0.05 * k) + ((k % 50 < 25) ? 0.5 : -0.5);
            const double u_model = (dss.C * x)(0, 0) + (dss.D(0, 0) * e);
            CHECK(ctrl.step(e) == doctest::Approx(u_model).epsilon(1e-10).scale(1.0));
            x = (dss.A * x) + (dss.B * Matrix<1, 1, double>{{e}});
        }
        // The lead first, then the lag, whose output is the command.
        CHECK(kLeadLag.gains[0].z < kLeadLag.gains[0].p);
        CHECK(kLeadLag.gains[1].z > kLeadLag.gains[1].p);
    }

    TEST_CASE_TEMPLATE("track then step returns u, clamped into the limits, for a lead, a lag and the series", T, float, double) {
        const auto lead = design::lead({.phase_max = 0.6, .wc = 20.0}, 1e-3, design::InputLimits{.lower = {-2.0}, .upper = {2.0}});
        const auto lag = design::lag({.dc_gain_boost = 10.0, .wc = 20.0}, 1e-3, design::InputLimits{.lower = {-2.0}, .upper = {2.0}});
        const T    tol = T{64} * eps<T>();

        const auto lands = [&](auto ctrl) {
            for (const T e : {static_cast<T>(0.3), static_cast<T>(-0.7)}) {
                for (const T u : {static_cast<T>(1.25), static_cast<T>(-0.4), T{9}}) {
                    CAPTURE(e);
                    CAPTURE(u);
                    for (int k = 0; k < 20; ++k) {
                        ctrl.track(u, e); // the input applied, the error read
                    }
                    CHECK(ctrl.last_command() == u);
                    ctrl.back_calculate(u); // the input track recorded: no clamp to report
                    const T want = damp::clamp(u, T{-2}, T{2});
                    CHECK(damp::abs(ctrl.step(e) - want) <= tol * T{4});
                }
            }
        };
        lands(LeadLagController<T>{lead.template as<T>()});
        lands(LeadLagController<T>{lag.template as<T>()});
        lands(LeadLagController<T, 2>{kLeadLag.template as<T>()});

        // track(u, r, y) is track(u, r − y).
        LeadLagController<T> a{lag.template as<T>()};
        LeadLagController<T> b{lag.template as<T>()};
        a.track(T{1}, T{2}, static_cast<T>(1.5));
        b.track(T{1}, static_cast<T>(0.5));
        CHECK(damp::abs(a.step(T{2}, static_cast<T>(1.5)) - b.step(static_cast<T>(0.5))) <= T{8} * eps<T>());
    }

    TEST_CASE("a lag that tracked a constant input hands over from it, then follows its own dynamics") {
        // Without track, the first step after an operator held u = 1.5 is the lag's response from rest.
        const auto                lag = design::lag({.dc_gain_boost = 10.0, .wc = 20.0}, 1e-3);
        LeadLagController<double> handed{lag};
        LeadLagController<double> cold{lag};
        for (int k = 0; k < 100; ++k) {
            handed.track(1.5, 0.1);
        }
        const double first = handed.step(0.1);
        CHECK(first == doctest::Approx(1.5).epsilon(1e-12));
        CHECK(cold.step(0.1) == doctest::Approx(0.1).epsilon(1e-2)); // K = 1 at high frequency: a kick from 1.5 to 0.1
        double u = first;
        for (int k = 0; k < 80000; ++k) { // 16 time constants of the lag's pole at 0.2 rad/s
            u = handed.step(0.1);
        }
        CHECK(u == doctest::Approx(1.0).epsilon(1e-6)); // settled on the DC gain, 10·0.1
    }

    TEST_CASE_TEMPLATE("set_limits and set_gains take effect at the next tick; a value out of its domain is refused and changes nothing", T, float, double) {
        LeadLagController<T, 2> ll{kLeadLag.template as<T>()};
        REQUIRE(ll.valid());

        auto lim = ll.limits();
        CHECK(lim.input.upper(0) == T{2});
        lim.input = {.lower = static_cast<T>(-0.5), .upper = static_cast<T>(0.5)};
        REQUIRE(ll.set_limits(lim));
        CHECK(ll.step(T{10}) == static_cast<T>(0.5));
        CHECK_FALSE(ll.set_limits({.input = {.lower = T{1}, .upper = T{-1}}}));
        CHECK(ll.limits().input.lower(0) == static_cast<T>(-0.5));

        auto g = ll.gains();
        g[1].z = T{4};
        REQUIRE(ll.set_gains(g));
        CHECK(ll.gains()[1].z == T{4});
        for (const auto& [i, bad] : {std::pair{0, design::LeadLagGains<T>{.K = T{1}, .z = T{-1}, .p = T{5}}}, std::pair{1, design::LeadLagGains<T>{.K = T{1}, .z = T{1}, .p = T{0}}}}) {
            auto refused = ll.gains();
            refused[static_cast<size_t>(i)] = bad;
            CHECK_FALSE(ll.set_gains(refused));
        }
        CHECK(ll.gains()[1].z == T{4});
        CHECK(ll.gains()[0].p == static_cast<T>(kLeadLag.gains[0].p));
    }

    TEST_CASE("Law<kRes, T> of a lead-lag design runs as its LeadLagController") {
        static constinit Law<kLeadLag, float> law;
        static_assert(std::is_base_of_v<LeadLagController<float, 2>, Law<kLeadLag, float>>);
        LeadLagController<float, 2> same{kLeadLag};
        CHECK(law.valid());
        CHECK(law.limits().input.upper(0) == 2.0f);
        for (int k = 0; k < 50; ++k) {
            const float e = static_cast<float>(std::sin(0.1 * k));
            CHECK(law.step(e) == doctest::Approx(same.step(e)).epsilon(1e-6)); // the same law: a few ulps of float
        }
    }
}
