// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// A NaN parameter, or a +Inf one where the parameter must be finite, is refused as documented
// under -ffast-math too. -ffinite-math-only lets the compiler assume no NaN or Inf, so it may
// answer any comparison on one either way: g++ -O3 folds !(x > 0) into x <= 0 and a NaN passes.
// Each guard here tests finiteness by the bit pattern (damp::isfinite) before it compares. The
// values arrive at run time, from bits through a volatile, as a sensor's or a config file's would:
// a NaN or Inf constant could be folded at compile time.

#include <bit>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "damp/analysis/disk_margin.hpp"
#include "damp/analysis/frequency.hpp"
#include "damp/analysis/linspace.hpp"
#include "damp/analysis/norms.hpp"
#include "damp/analysis/time_response.hpp"
#include "damp/backend.hpp"
#include "damp/controllers/adrc.hpp"
#include "damp/controllers/esc.hpp"
#include "damp/controllers/lead_lag.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/controllers/smc.hpp"
#include "damp/controllers/smith_predictor.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/controllers/stsmc.hpp"
#include "damp/design/model_reduction.hpp"
#include "damp/design/options.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/stability.hpp"
#include "damp/design/status.hpp"
#include "damp/encoder/gate.hpp"
#include "damp/estimation/dob.hpp"
#include "damp/estimation/eskf.hpp"
#include "damp/estimation/eso.hpp"
#include "damp/estimation/excitation/prbs.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/frequency_response.hpp"
#include "damp/estimation/ins_eskf.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/parameter_estimation.hpp"
#include "damp/estimation/relay_autotune.hpp"
#include "damp/estimation/rls.hpp"
#include "damp/estimation/sensor_fusion.hpp"
#include "damp/estimation/state_variable_filter.hpp"
#include "damp/estimation/ukf.hpp"
#include "damp/filters/differentiator.hpp"
#include "damp/filters/fir.hpp"
#include "damp/filters/iir_design.hpp"
#include "damp/filters/lowpass.hpp"
#include "damp/filters/moving_average.hpp"
#include "damp/filters/sogi.hpp"
#include "damp/filters/spectral.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "damp/math/split_real.hpp"
#include "damp/math/transforms.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/motor/modulation.hpp"
#include "damp/motor/spm.hpp"
#include "damp/simulation/cached_zoh.hpp"
#include "damp/simulation/hybrid.hpp"
#include "damp/simulation/integrator.hpp"
#include "damp/simulation/quantize.hpp"
#include "damp/simulation/simulate.hpp"
#include "damp/simulation/solver.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/toolbox/bounds.hpp"
#include "damp/toolbox/conditioning.hpp"
#include "damp/toolbox/thermistor.hpp"
#include "damp/toolbox/timing.hpp"
#include "damp/trajectory/trajectory_types.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using design::DesignStatus;

namespace {

/// A quiet NaN from its bits, read through a volatile so the compiler cannot fold it.
template<typename T>
T runtime_nan() {
    if constexpr (sizeof(T) == 4) {
        volatile std::uint32_t bits = 0x7fc00000u;
        return std::bit_cast<float>(static_cast<std::uint32_t>(bits));
    } else {
        volatile std::uint64_t bits = 0x7ff8000000000000ull;
        return std::bit_cast<double>(static_cast<std::uint64_t>(bits));
    }
}

/// +Inf from its bits, read through a volatile.
template<typename T>
T runtime_inf() {
    if constexpr (sizeof(T) == 4) {
        volatile std::uint32_t bits = 0x7f800000u;
        return std::bit_cast<float>(static_cast<std::uint32_t>(bits));
    } else {
        volatile std::uint64_t bits = 0x7ff0000000000000ull;
        return std::bit_cast<double>(static_cast<std::uint64_t>(bits));
    }
}

/// The values a finite parameter refuses: NaN and +Inf (index 0 and 1).
template<typename T>
std::vector<T> non_finite() {
    return {runtime_nan<T>(), runtime_inf<T>()};
}

/// The bits of a float or double, to check a value was left exactly as it was.
template<typename T>
auto bits_of(T x) {
    if constexpr (sizeof(T) == 4) {
        return std::bit_cast<std::uint32_t>(x);
    } else {
        return std::bit_cast<std::uint64_t>(x);
    }
}

} // namespace

TEST_SUITE("NaN guards") {
    TEST_CASE_TEMPLATE("IIR designs refuse a non-finite cutoff, Q or sample time", T, float, double) {
        const T ts{static_cast<T>(1e-3)};
        for (const T bad : non_finite<T>()) {
            CHECK(design::lowpass_1st(bad).status == DesignStatus::invalid_parameter);
            CHECK(design::highpass_1st(bad).status == DesignStatus::invalid_parameter);
            CHECK(design::lowpass_2nd(bad).status == DesignStatus::invalid_parameter);
            CHECK(design::lowpass_2nd(T{10}, bad).status == DesignStatus::invalid_parameter);
            CHECK(design::highpass_2nd(bad).status == DesignStatus::invalid_parameter);
            CHECK(design::highpass_2nd(T{10}, bad).status == DesignStatus::invalid_parameter);
            CHECK(design::bandpass(bad, T{1}).status == DesignStatus::invalid_parameter);
            CHECK(design::bandpass(T{10}, bad).status == DesignStatus::invalid_parameter);
            CHECK(design::notch(bad, T{1}).status == DesignStatus::invalid_parameter);
            CHECK(design::notch(T{10}, bad).status == DesignStatus::invalid_parameter);
            CHECK(design::allpass(bad, T{1}).status == DesignStatus::invalid_parameter);
            CHECK(design::allpass(T{10}, bad).status == DesignStatus::invalid_parameter);

            CHECK(design::lowpass_1st(bad, ts).status == DesignStatus::invalid_parameter);
            CHECK(design::lowpass_1st(T{10}, bad).status == DesignStatus::invalid_sample_time);
            CHECK(design::lowpass_2nd(T{10}, bad, ts).status == DesignStatus::invalid_parameter);
            CHECK(design::notch(T{10}, T{1}, bad).status == DesignStatus::invalid_sample_time);
            CHECK(design::bandpass(bad, T{1}, ts).status == DesignStatus::invalid_parameter);

            const TransferFunction<2, 2, T> first{.num = {T{0}, T{1}}, .den = {T{1}, T{1}}};
            const TransferFunction<1, 2, T> first_pole{.num = {T{1}}, .den = {T{1}, T{1}}};
            const TransferFunction<3, 3, T> second{.num = {T{0}, T{0}, T{1}}, .den = {T{1}, T{1}, T{1}}};
            CHECK(design::to_coeffs(first, bad).status == DesignStatus::invalid_sample_time);
            CHECK(design::to_coeffs(first_pole, bad).status == DesignStatus::invalid_sample_time);
            CHECK(design::to_coeffs(second, bad).status == DesignStatus::invalid_sample_time);
            CHECK(design::to_coeffs(StateSpace<1, 1, 1, T>{}, bad).status == DesignStatus::invalid_sample_time);
            CHECK(design::to_coeffs(StateSpace<2, 1, 1, T>{}, bad).status == DesignStatus::invalid_sample_time);
        }
        CHECK(design::lowpass_2nd(T{10}, T{1}, ts).success); // the finite design is unchanged
    }

    TEST_CASE_TEMPLATE("FIR, moving-average and continuous low-pass helpers refuse a non-finite parameter", T, float, double) {
        for (const T bad : non_finite<T>()) {
            CHECK(design::fir1<16, T>(9, bad).status == DesignStatus::invalid_parameter);
            CHECK(comb_notch_window(bad, static_cast<T>(1e-3)) == 0);
            CHECK(comb_notch_window(T{50}, bad) == 0);
            CHECK(continuous_lpf_exact_step(T{1}, T{3}, bad, static_cast<T>(0.01)) == T{1});
            CHECK(continuous_lpf_exact_step(T{1}, T{3}, T{10}, bad) == T{1});
        }
        CHECK(design::fir1<16, T>(9, static_cast<T>(0.2)).success);
        CHECK(comb_notch_window(T{50}, static_cast<T>(1e-3)) == 20);
    }

    TEST_CASE_TEMPLATE("SOGI-FLL, Goertzel and the robust differentiator refuse a non-finite parameter", T, float, double) {
        const T ts{static_cast<T>(1e-4)};
        for (const T bad : non_finite<T>()) {
            CHECK_FALSE(SogiFll<T>(bad, ts).valid());
            CHECK_FALSE(SogiFll<T>(T{50}, ts, bad).valid());
            CHECK(Goertzel<T>(bad, ts, 64).block_size() == 0);
            CHECK_FALSE(RobustExactDifferentiator<T>(bad, static_cast<T>(0.01)).valid());
            CHECK_FALSE(RobustExactDifferentiator<T>(T{10}, bad).valid());
        }
        // A NaN end of the band falls back to its default, as an end ≤ 0 does.
        const T    nan = runtime_nan<T>();
        SogiFll<T> fll(T{50}, ts, numbers::sqrt2_v<T>, T{2}, Bounds<1, T>{.lower = ColVec<1, T>{nan}, .upper = ColVec<1, T>{nan}});
        CHECK(fll.valid());
        for (int k = 0; k < 100; ++k) {
            fll.step(damp::sin(T{2} * numbers::pi_v<T> * T{50} * static_cast<T>(k) * ts));
        }
        CHECK(damp::isfinite(fll.frequency_hz()));
        CHECK(SogiFll<T>(T{50}, ts).valid());
        CHECK(Goertzel<T>(T{50}, ts, 64).block_size() == 64);
    }

    TEST_CASE_TEMPLATE("controllers built from a result with a non-finite sample time, pole, zero or gain are invalid", T, float, double) {
        const T ts{static_cast<T>(1e-3)};
        for (const T bad : non_finite<T>()) {
            auto lead = design::lead(static_cast<T>(0.5), T{10}, ts);
            REQUIRE(lead.success);
            CHECK(LeadLagController<T>{lead}.valid());
            for (T design::LeadLagResult<T>::* field : {&design::LeadLagResult<T>::Ts, &design::LeadLagResult<T>::z, &design::LeadLagResult<T>::p}) {
                auto broken = lead;
                broken.*field = bad;
                LeadLagController<T> c{broken};
                CHECK_FALSE(c.valid());
                CHECK(c.control(T{1}, T{0}) == T{0});
            }

            auto pr = design::pr(design::PRSpec<T>{.Kp = T{1}, .Ki = T{10}, .w0 = T{314}}, static_cast<T>(1e-4));
            REQUIRE(pr.success);
            pr.Ts = bad;
            CHECK_FALSE(PRController<T>{pr}.valid());

            auto smc = design::smc(T{5}, T{2}, T{1}, ts, static_cast<T>(0.1));
            REQUIRE(smc.success);
            auto smc_ts = smc;
            smc_ts.Ts = bad;
            CHECK_FALSE(SMCController<T>{smc_ts}.valid());
            auto smc_phi = smc;
            smc_phi.phi = bad;
            CHECK(smc_phi.to_tf().num[0] == T{0}); // the zero TF of a design without a boundary layer
            auto smc_b0 = smc;
            smc_b0.b0 = bad;
            CHECK(smc_b0.to_tf().num[0] == T{0});

            auto stsmc = design::stsmc(T{1}, ts);
            REQUIRE(stsmc.success);
            stsmc.Ts = bad;
            CHECK_FALSE(STSMCController<T>{stsmc}.valid());

            auto smith = design::smith_predictor_from_fopdt(T{2}, T{1}, static_cast<T>(0.1), static_cast<T>(0.01));
            REQUIRE(smith.success);
            auto smith_ts = smith;
            smith_ts.Ts = bad;
            CHECK(smith_ts.primary_pid().Kp == T{0});
            SmithPredictor<8, T> hand{PIController<T>{smith.pid}, static_cast<T>(0.9), bad, 2};
            CHECK_FALSE(hand.valid());
        }
    }

    TEST_CASE_TEMPLATE("ADRC refuses a non-finite b0 in its result and its runtime", T, float, double) {
        const auto good = design::adrc<1>(T{10}, T{40}, T{2}, static_cast<T>(1e-3));
        REQUIRE(good.success);
        for (const T bad : non_finite<T>()) {
            auto broken = good;
            broken.b0 = bad;
            const auto sys = broken.to_ss_2dof();
            CHECK(sys.C(0, 0) == T{0}); // the zero system
            CHECK(sys.D(0, 0) == T{0});
            REQUIRE(ADRCController<1, T>{good}.valid());
            ADRCController<1, T> c{broken};
            CHECK_FALSE(c.valid());
            CHECK(c.control(T{1}, T{0}) == T{0});
        }
    }

    TEST_CASE_TEMPLATE("PID helpers refuse a non-finite filter or damping", T, float, double) {
        for (const T bad : non_finite<T>()) {
            design::PIDGains<T> g{};
            g.Kp = T{1};
            g.Ki = T{1};
            g.Kd = static_cast<T>(0.1);
            g.Tf = bad;
            CHECK_FALSE(g.to_ss_2dof().has_value());
            CHECK(design::phase_margin_from_damping_ratio(bad) == T{30});
        }
        CHECK(design::phase_margin_from_damping_ratio(T{1}) > T{70});
    }

    TEST_CASE("AMIGO refuses a relay point whose finite inputs overflow its model") {
        // Ks·Ku overflows to +Inf, so τ is +Inf and the gains NaN: they were returned as a success.
        volatile double big = 1e200;
        const auto      gains = design::amigo_kappa_tau(static_cast<double>(big), static_cast<double>(big), 1.0);
        CHECK(gains.status == DesignStatus::invalid_parameter);
        CHECK(design::amigo_kappa_tau(0.5, 4.0, 1.0).success);
    }

    TEST_CASE_TEMPLATE("MPC horizon suggestions refuse a non-finite settling time or sample time", T, float, double) {
        for (const T bad : non_finite<T>()) {
            CHECK_FALSE(design::suggest_mpc_horizon(bad, static_cast<T>(0.01)).success);
            CHECK_FALSE(design::suggest_mpc_horizon(T{1}, bad).success);
            const StateSpace<1, 1, 1, T> plant{.A = {{static_cast<T>(0.9)}}, .B = {{T{1}}}, .C = {{T{1}}}, .Ts = bad};
            CHECK_FALSE(design::suggest_mpc_horizon(plant).success);
        }
        CHECK(design::suggest_mpc_horizon(T{1}, static_cast<T>(0.01)).success);
    }

    TEST_CASE("place_jordan refuses a NaN robustness weight with bad_weights") {
        const StateSpace<2, 1, 1, double> sys{
            .A = {{1.0, 0.1}, {0.0, 1.0}},
            .B = {{0.005}, {0.1}},
            .C = {{1.0, 0.0}},
            .Ts = 0.1,
        };
        const damp::array<design::JordanBlock<double>, 1> blocks{{{{0.5, 0.0}, 2}}};
        CHECK(design::place_jordan(sys, blocks, design::JordanOptions<double>{.alpha = runtime_nan<double>()}).status == DesignStatus::bad_weights);
        CHECK(design::place_jordan(sys, blocks, design::JordanOptions<double>{.alpha = runtime_inf<double>()}).status == DesignStatus::bad_weights);
        CHECK(design::place_jordan(sys, blocks).success);
    }

    TEST_CASE("branch-and-bound floor passes a NaN through instead of casting it") {
        CHECK_FALSE(damp::isfinite(design::detail::floor_scalar(runtime_nan<double>())));
        CHECK(design::detail::floor_scalar(2.5) == 2.0);
        CHECK(design::detail::floor_scalar(-2.5) == -3.0);
    }

    TEST_CASE_TEMPLATE("observers built from a result with a non-finite sample time are invalid", T, float, double) {
        for (const T bad : non_finite<T>()) {
            auto eso = design::eso<1>(T{40}, T{2}, static_cast<T>(1e-3));
            REQUIRE(eso.success);
            eso.Ts = bad;
            CHECK_FALSE(ESO<1, T>{eso}.valid());

            const StateSpace<2, 1, 1, T> sys{
                .A = {{T{1}, static_cast<T>(0.1)}, {T{0}, T{1}}},
                .B = {{static_cast<T>(0.005)}, {static_cast<T>(0.1)}},
                .C = {{T{1}, T{0}}},
                .Ts = static_cast<T>(0.1),
            };
            auto reduced = design::reduced_order_observer(sys, {static_cast<T>(0.5)});
            REQUIRE(reduced.success);
            reduced.Ts = bad;
            CHECK_FALSE(ReducedOrderObserver<2, 1, T>(reduced).valid());

            auto placed = design::place_observer(sys, {static_cast<T>(0.3), static_cast<T>(0.4)});
            REQUIRE(placed.success);
            CHECK(FixedGainObserver<2, 1, 1, T>{placed}.valid());
            placed.sys.Ts = bad;
            CHECK_FALSE(FixedGainObserver<2, 1, 1, T>{placed}.valid());
        }
    }

    TEST_CASE_TEMPLATE("DOB configs and the classical DOB refuse a non-finite entry", T, float, double) {
        for (const T bad : non_finite<T>()) {
            for (T estimation::DOBConfig<T>::* field : {&estimation::DOBConfig<T>::gain, &estimation::DOBConfig<T>::leak, &estimation::DOBConfig<T>::innovation_deadband, &estimation::DOBConfig<T>::max_disturbance_magnitude}) {
                estimation::DOBConfig<T> cfg{};
                cfg.gain = static_cast<T>(0.5);
                REQUIRE(cfg.valid());
                cfg.*field = bad;
                CHECK_FALSE(cfg.valid());
            }
            const damp::array<T, 2> Bn{T{1}, static_cast<T>(0.5)};
            const damp::array<T, 2> An{T{1}, static_cast<T>(-0.9)};
            const damp::array<T, 1> Qn{static_cast<T>(0.1)};
            const damp::array<T, 2> Qd{T{1}, static_cast<T>(-0.9)};
            REQUIRE(estimation::classical_dob(Bn, An, Qn, Qd).success);
            CHECK_FALSE(estimation::classical_dob(damp::array<T, 2>{bad, static_cast<T>(0.5)}, An, Qn, Qd).success);
            CHECK_FALSE(estimation::classical_dob(Bn, damp::array<T, 2>{bad, static_cast<T>(-0.9)}, Qn, Qd).success);
            CHECK_FALSE(estimation::classical_dob(Bn, An, Qn, damp::array<T, 2>{T{1}, bad}).success);
            CHECK_FALSE(tf_detail::schur_stable(damp::array<T, 3>{T{1}, bad, static_cast<T>(0.25)}));
        }
    }

    TEST_CASE_TEMPLATE("frequency-response bins refuse a non-finite frequency or sample rate", T, float, double) {
        for (const T bad : non_finite<T>()) {
            CHECK(FrequencyResponseEstimator<4, T>::coherent_block_size(bad, T{1000}, 10) == 1);
            CHECK(FrequencyResponseEstimator<4, T>::coherent_block_size(T{50}, bad, 10) == 1);
            FrequencyResponseEstimator<4, T> fre{T{1000}};
            CHECK_FALSE(fre.configure_bin(0, bad, 100));
            CHECK_FALSE(fre.begin_bin(0, bad, 100));
            FrequencyResponseEstimator<4, T> no_rate{bad};
            CHECK_FALSE(no_rate.configure_bin(0, T{10}, 100));
        }
        FrequencyResponseEstimator<4, T> fre{T{1000}};
        CHECK(fre.configure_bin(0, T{10}, 100));
    }

    TEST_CASE_TEMPLATE("INS and attitude steps hold on a non-finite dt", T, float, double) {
        for (const T bad : non_finite<T>()) {
            InsState<T> x{};
            x.v = Vec3<T>{T{1}, T{2}, T{3}};
            const Vec3<T>     omega{static_cast<T>(0.1), T{0}, T{0}};
            const Vec3<T>     a_b{T{0}, T{0}, static_cast<T>(-9.8)};
            const InsState<T> y = mechanize_step_from_corrected(x, omega, a_b, bad, NavFrame::NED);
            CHECK(bits_of(y.v[0]) == bits_of(x.v[0]));
            CHECK(bits_of(y.p[2]) == bits_of(x.p[2]));
            const auto F = ins_error_jacobian(x.q, omega, a_b, bad);
            CHECK(F(0, 0) == T{1}); // the identity
            CHECK(F(3, 6) == T{0});

            ErrorStateKalmanFilter<6, T> eskf{{.Q = Matrix<6, 6, T>::identity(), .P0 = Matrix<6, 6, T>::identity()}};
            Quaternion<T>                q = Quaternion<T>::identity();
            Vec3<T>                      b_g{};
            detail::eskf_attitude_predict(eskf, q, b_g, omega, bad);
            CHECK(q.w() == T{1});
            CHECK(eskf.covariance()(0, 0) == T{1});
            CHECK_FALSE(detail::eskf_accel_gate_ok(Vec3<T>{T{0}, T{0}, T{1}}, Vec3<T>{T{0}, T{0}, bad}, static_cast<T>(0.1)));
        }
    }

    TEST_CASE_TEMPLATE("RLS, parameter-estimator and relay-autotune configs refuse a non-finite entry", T, float, double) {
        for (const T bad : non_finite<T>()) {
            CHECK_FALSE(estimation::RlsConfig<T>{.lambda = bad}.valid());
            CHECK_FALSE(estimation::RlsConfig<T>{.p0 = bad}.valid());

            for (T design::FirstOrderPlantEstimatorConfig<T>::* field :
                 {&design::FirstOrderPlantEstimatorConfig<T>::Ts, &design::FirstOrderPlantEstimatorConfig<T>::forgetting,
                  &design::FirstOrderPlantEstimatorConfig<T>::initial_covariance, &design::FirstOrderPlantEstimatorConfig<T>::confidence_scale}) {
                design::FirstOrderPlantEstimatorConfig<T> cfg{.Ts = static_cast<T>(0.01)};
                REQUIRE(cfg.valid());
                cfg.*field = bad;
                CHECK_FALSE(cfg.valid());
            }

            for (T design::RelayAutotuneConfig<T>::* field :
                 {&design::RelayAutotuneConfig<T>::amplitude, &design::RelayAutotuneConfig<T>::amplitude_neg, &design::RelayAutotuneConfig<T>::hysteresis,
                  &design::RelayAutotuneConfig<T>::u_bias, &design::RelayAutotuneConfig<T>::period_tolerance, &design::RelayAutotuneConfig<T>::max_duration}) {
                design::RelayAutotuneConfig<T> cfg{};
                REQUIRE(cfg.valid());
                cfg.*field = bad;
                CHECK_FALSE(cfg.valid());
            }
            design::RelayAutotuneConfig<T> nan_lim{};
            nan_lim.u_lim.upper[0] = runtime_nan<T>();
            CHECK_FALSE(nan_lim.valid());

            auto relay = design::relay_autotune(design::RelayAutotuneConfig<T>{}, static_cast<T>(0.01));
            REQUIRE(relay.success);
            relay.Ts = bad;
            CHECK(RelayAutotuner<T>{relay}.status() == RelayAutotuneStatus::Failed);
        }
    }

    TEST_CASE_TEMPLATE("RLS skips a reading whose regressor is NaN and keeps its estimate", T, float, double) {
        estimation::RlsVector<1, T> rls(estimation::rls_vector(estimation::RlsConfig<T>{}));
        REQUIRE(rls.update(ColVec<1, T>{T{1}}, T{2}));
        const T theta = rls.state().theta[0];
        const T p = rls.state().covariance(0, 0);
        CHECK_FALSE(rls.update(ColVec<1, T>{runtime_nan<T>()}, T{2}));
        CHECK(bits_of(rls.state().theta[0]) == bits_of(theta));
        CHECK(bits_of(rls.state().covariance(0, 0)) == bits_of(p));
    }

    TEST_CASE_TEMPLATE("state-variable filter, PRBS and UKF weights refuse a non-finite parameter", T, float, double) {
        for (const T bad : non_finite<T>()) {
            CHECK_FALSE(estimation::state_variable_filter(bad, static_cast<T>(0.01)).success);
            CHECK_FALSE(estimation::state_variable_filter(T{10}, bad).success);

            const design::PRBSConfig<T> cfg{.amplitude = T{1}, .lfsr_order = 5, .clock_period_s = static_cast<T>(0.1), .seed = 0x1u};
            auto                        prbs = design::prbs(cfg, static_cast<T>(0.1));
            REQUIRE(prbs.success);
            prbs.Ts = bad;
            PRBS<T> seq{prbs};
            for (int k = 0; k < 4; ++k) {
                CHECK(seq.step() == T{0});
            }
        }
        // κ = NaN makes λ = α²(n + κ) − n NaN: the weights fall back to the plain unscented transform.
        UnscentedParams<T> p{};
        p.kappa = runtime_nan<T>();
        struct Still {
            [[nodiscard]] ColVec<2, T> operator()(const ColVec<2, T>& x, const ColVec<1, T>& /*u*/) const { return x; }
        };
        struct First {
            [[nodiscard]] ColVec<1, T> operator()(const ColVec<2, T>& x, const ColVec<1, T>& /*u*/) const { return ColVec<1, T>{x[0]}; }
        };
        UnscentedKalmanFilter<Still> ukf{Still{}, {.Ts = static_cast<T>(0.01), .Q = Matrix<2, 2, T>::identity() * static_cast<T>(0.01), .x0 = ColVec<2, T>{T{1}, T{2}}, .unscented = p}};
        CHECK(ukf.predict({}));
        CHECK(ukf.update(First{}, ColVec<1, T>{static_cast<T>(1.5)}, static_cast<T>(0.1)));
        CHECK(damp::isfinite(ukf.state()[0]));
        CHECK(damp::isfinite(ukf.state()[1]));
    }

    TEST_CASE_TEMPLATE("sample-time matching and polynomial stability refuse NaN", T, float, double) {
        const T nan = runtime_nan<T>();
        CHECK_FALSE(detail::same_sample_time(nan, nan));
        CHECK_FALSE(detail::same_sample_time(nan, static_cast<T>(0.1)));
        CHECK_FALSE(detail::same_sample_time(runtime_inf<T>(), runtime_inf<T>()));
        CHECK(detail::same_sample_time(static_cast<T>(0.1), static_cast<T>(0.1)));
        CHECK_FALSE(tf_detail::schur_stable(damp::array<T, 2>{nan, static_cast<T>(0.5)}));
        CHECK(tf_detail::schur_stable(damp::array<T, 2>{T{1}, static_cast<T>(0.5)}));
    }

    TEST_CASE_TEMPLATE("designs for a continuous plant refuse a plant whose Ts is NaN", T, float, double) {
        // x86 under -ffast-math reads NaN == 0 as true: a test sys.Ts != 0 took a NaN for continuous.
        StateSpace<2, 1, 1, T> plant{.A = {{T{0}, T{1}}, {T{-2}, T{-3}}}, .B = {{T{0}}, {T{1}}}, .C = {{T{1}, T{0}}}};
        const T                ts{static_cast<T>(0.01)};
        REQUIRE(design::pidtune(plant, design::PIDTuneSpec<T>{.wc = T{2}}).success);
        REQUIRE(design::kalman(plant, ts, design::KalmanNoise{.process = Matrix<2, 2, T>::identity(), .measurement = static_cast<T>(0.01)}).success);
        REQUIRE(design::state_feedback(plant, ts, design::LQRWeights{.Q = Matrix<2, 2, T>::identity(), .R = {T{1}}}).success);
        plant.Ts = runtime_nan<T>();
        CHECK(design::pidtune(plant, design::PIDTuneSpec<T>{.wc = T{2}}).status == DesignStatus::invalid_sample_time);
        CHECK(design::kalman(plant, ts, design::KalmanNoise{.process = Matrix<2, 2, T>::identity(), .measurement = static_cast<T>(0.01)}).status == DesignStatus::invalid_sample_time);
        CHECK(design::state_feedback(plant, ts, design::LQRWeights{.Q = Matrix<2, 2, T>::identity(), .R = {T{1}}}).status == DesignStatus::invalid_sample_time);
        CHECK_FALSE(discretize(plant, ts, InputDelay{static_cast<T>(0.005)}).success);
    }

    TEST_CASE_TEMPLATE("ESC for MPPT refuses a NaN low-pass corner instead of taking the default", T, float, double) {
        design::ESCSpec<T> spec{.dither_amplitude = static_cast<T>(0.1), .dither_freq_hz = T{50}, .gain = T{20}};
        REQUIRE(design::esc_mppt(spec, static_cast<T>(1e-3)).success);
        spec.lp_cutoff_hz = runtime_nan<T>();
        CHECK(design::esc_mppt(spec, static_cast<T>(1e-3)).status == DesignStatus::invalid_parameter);
    }

    TEST_CASE_TEMPLATE("toolbox blocks hold or fall back on a non-finite parameter", T, float, double) {
        const Thermistor<T> ntc{design::thermistor_beta(T{10000}, static_cast<T>(298.15), T{3950})};
        for (const T bad : non_finite<T>()) {
            CHECK(scaled_deadband(static_cast<T>(0.5), bad) == T{0});

            SlewLimiter<T> slew{T{1}};
            (void)slew(T{0}, static_cast<T>(0.1));
            CHECK(slew(T{5}, bad) == T{0}); // held

            CHECK(ntc.kelvin(bad) == numeric_limits<T>::max()); // reads hottest, so protection trips

            Periodic<T> tick{bad};
            for (int k = 0; k < 4; ++k) {
                CHECK_FALSE(tick(T{1}));
            }

            CHECK(smooth_sign(static_cast<T>(0.3), bad) == T{1}); // the plain sign
            CHECK(wrap(T{1}, T{0}, bad) == T{0});
            CHECK_FALSE(TrajectoryLimits<T>{.max_velocity = bad, .max_acceleration = T{1}, .max_deceleration = T{1}}.valid());
            CHECK_FALSE(TrajectoryLimits<T>{.max_velocity = T{1}, .max_acceleration = T{1}, .max_deceleration = T{1}, .max_jerk = bad}.valid_jerk_limited());
        }
        CHECK(scaled_deadband(T{1}, static_cast<T>(0.5)) == doctest::Approx(1.0));
        CHECK(ntc.kelvin(T{10000}) == doctest::Approx(298.15));
        CHECK(TrajectoryLimits<T>{.max_velocity = T{1}, .max_acceleration = T{1}, .max_deceleration = T{1}, .max_jerk = T{1}}.valid_jerk_limited());
    }

    TEST_CASE_TEMPLATE("encoder gate with a non-finite tolerance is off", T, float, double) {
        for (const T bad : non_finite<T>()) {
            EncoderGate<T> gate{{.tolerance = bad, .resync_count = 1}};
            REQUIRE(gate.check(SplitReal<T>{0, T{0}}, T{0}, static_cast<T>(1e-3)) == EncoderVerdict::Resync);
            for (int k = 1; k < 4; ++k) {
                // A jump of 100 counts at rest: a finite tolerance drops it, an unusable one is no window.
                CHECK(gate.check(SplitReal<T>{100 * k, T{0}}, T{0}, static_cast<T>(1e-3)) == EncoderVerdict::Accepted);
            }
        }
    }

    TEST_CASE_TEMPLATE("modulation parks mid-rail on a non-finite bus voltage or command", T, float, double) {
        for (const T bad : non_finite<T>()) {
            const auto duties = motor::modulation_duty_cycles<T>(AlphaBeta<T>{T{1}, T{0}}, bad);
            CHECK(duties.is_clipped);
            CHECK(duties.duties[0] == static_cast<T>(0.5));
            CHECK(motor::modulation_index<T>(AlphaBeta<T>{T{1}, T{0}}, bad) == T{0});
        }
        const T    nan = runtime_nan<T>();
        const auto six = motor::six_step_duty_cycles<T>(AlphaBeta<T>{nan, T{0}});
        CHECK(six.is_clipped);
        CHECK(six.duties[1] == static_cast<T>(0.5));
        CHECK(motor::thipwm_zero_sequence<T>(AlphaBeta<T>{nan, T{0}}) == T{0});
    }

    TEST_CASE_TEMPLATE("SPM envelope helpers return no current for a non-finite parameter or speed", T, float, double) {
        const T we{T{100}};
        const T v{T{24}};
        const T i_max{T{10}};
        const T lq{static_cast<T>(1e-3)};
        const T r{static_cast<T>(0.1)};
        const T psi{static_cast<T>(0.01)};
        REQUIRE(motor::spm::max_iq(we, v, i_max, lq, r, psi) > T{0});
        for (const T bad : non_finite<T>()) {
            CHECK(motor::spm::max_iq(bad, v, i_max, lq, r, psi) == T{0});
            CHECK(motor::spm::max_iq(we, bad, i_max, lq, r, psi) == T{0});
            CHECK(motor::spm::max_iq(we, v, bad, lq, r, psi) == T{0});
            CHECK(motor::spm::max_iq(we, v, i_max, bad, r, psi) == T{0});
            CHECK(motor::spm::max_iq(we, v, i_max, lq, r, bad) == T{0});
            CHECK(motor::spm::max_torque_at_speed(T{10}, v, i_max, lq, r, psi, bad).iq == T{0});
        }
    }

    TEST_CASE_TEMPLATE("host helpers refuse a non-finite step, tolerance or span", T, float, double) {
        for (const T bad : non_finite<T>()) {
            CHECK(analysis::arange(T{0}, T{1}, bad).empty()); // terminates: the loop's end test never sees NaN
            CHECK(analysis::arange(T{0}, bad, static_cast<T>(0.25)).empty());
            CHECK(analysis::arange(bad, T{1}, static_cast<T>(0.25)).empty());
            CHECK(analysis::geomspace(bad, T{10}, 5).empty());
            CHECK(analysis::geomspace(T{1}, bad, 5).empty());

            const StateSpace<1, 1, 1, T> sys{.A = {{T{-1}}}, .B = {{T{1}}}, .C = {{T{1}}}};
            CHECK_FALSE(analysis::norm_hinf(sys, bad).success);
            CHECK_FALSE(analysis::detail::discretize_on_grid(sys, std::vector<T>{T{0}, bad}).has_value());
            CHECK_FALSE(analysis::detail::discretize_on_grid(sys, std::vector<T>{T{0}, static_cast<T>(0.1), bad}).has_value());

            CHECK(sim::quantize(static_cast<T>(0.3), bad, T{-1}, T{1}) == static_cast<T>(0.3));
            CHECK(sim::detail::step_count(T{0}, T{1}, bad) == 0);
            const auto zoh = sim::CachedZoh<1, 1, T>::from(sys.A, sys.B, bad);
            CHECK(zoh.Ad(0, 0) == T{1}); // identity: the state holds
            CHECK(zoh.Bd(0, 0) == T{0});
        }
        CHECK(analysis::arange(T{0}, T{1}, static_cast<T>(0.25)).size() == 4);
        CHECK(sim::detail::step_count(T{0}, T{1}, static_cast<T>(0.25)) == 4);
    }

    TEST_CASE_TEMPLATE("parameters that had no guard are refused when not finite", T, float, double) {
        const T ts{static_cast<T>(1e-4)};
        auto    smith = design::smith_predictor_from_fopdt(T{2}, T{1}, static_cast<T>(0.1), static_cast<T>(0.01));
        REQUIRE(smith.success);
        for (const T bad : non_finite<T>()) {
            CHECK_FALSE(SogiFll<T>(T{50}, ts, numbers::sqrt2_v<T>, bad).valid());
            CHECK_FALSE(RobustExactDifferentiator<T>(T{10}, static_cast<T>(0.01), bad).valid());
            CHECK_FALSE(RobustExactDifferentiator<T>(T{10}, static_cast<T>(0.01), static_cast<T>(1.5), bad).valid());
            CHECK_FALSE(SmithPredictor<8, T>(PIController<T>{smith.pid}, bad, static_cast<T>(0.1), 2).valid());

            design::RelayAutotuneConfig<T> relay{};
            relay.setpoint = bad;
            CHECK_FALSE(relay.valid());

            estimation::RlsConfig<T> rls{.projection_enabled = true};
            REQUIRE(rls.valid());
            rls.theta_lim.upper[0] = bad;
            CHECK_FALSE(rls.valid());

            // Every gain a discrete PID runtime reads; a NaN Kd must not read as no derivative.
            design::DiscretePIDGains<T> g{};
            g.Kp = bad;
            CHECK_FALSE(PIDController<T, PIDMode::PID>{g}.valid());
            CHECK_FALSE(PIDController<T, PIDMode::PI>{g}.valid());
            CHECK_FALSE(PIDController<T, PIDMode::PD>{g}.valid());
            CHECK_FALSE(PIDController<T, PIDMode::P>{g}.valid());
            design::DiscretePIDGains<T> kd{};
            kd.Kp = T{1};
            kd.Kd = bad;
            CHECK_FALSE(PIDController<T, PIDMode::PI>{kd}.valid());
            CHECK_FALSE(PIDController<T, PIDMode::P>{kd}.valid());
        }
        CHECK(SogiFll<T>(T{50}, ts, numbers::sqrt2_v<T>, T{0}).valid()); // Γ = 0: a fixed frequency
        CHECK(SmithPredictor<8, T>(PIController<T>{smith.pid}, smith.a, smith.b, 2).valid());
        design::DiscretePIDGains<T> p{};
        p.Kp = T{1};
        CHECK(PIDController<T, PIDMode::P>{p}.valid());
    }

    TEST_CASE_TEMPLATE("disk margins refuse a non-finite skew", T, float, double) {
        const StateSpace<1, 1, 1, T> loop{.A = {{T{-1}}}, .B = {{T{1}}}, .C = {{T{2}}}};
        REQUIRE(analysis::disk_margin(loop).success);
        for (const T bad : non_finite<T>()) {
            const auto dm = analysis::disk_margin(loop, bad);
            CHECK_FALSE(dm.success);
            CHECK(dm.stable); // the loop is stable; only the disk has no size
            const auto mimo = analysis::loop_at_a_time_disk_margin(loop, bad);
            CHECK_FALSE(mimo.success);
            CHECK(mimo.stable);
        }
    }

    TEST_CASE("the adaptive ODE solver refuses options that are not finite") {
        const auto   decay = [](double /*t*/, const ColVec<1>& x) { return ColVec<1>{-x(0, 0)}; };
        sim::DP45<1> dp;
        REQUIRE(sim::AdaptiveStepSolver(dp, sim::AdaptiveOptions<1>{}).solve(decay, ColVec<1>{1.0}, {0.0, 1.0}).success);
        const double nan = runtime_nan<double>();
        const double inf = runtime_inf<double>();
        for (int field = 0; field < 8; ++field) {
            sim::AdaptiveOptions<1> o{};
            switch (field) {
                case 0: o.atol = nan; break;
                case 1: o.rtol = nan; break;
                case 2: o.first_step = nan; break;
                case 3: o.min_step = nan; break;
                case 4: o.max_step = inf; break; // unbounded is max(), the default
                case 5: o.safety = nan; break;
                case 6: o.min_factor = nan; break;
                default: o.max_factor = nan; break;
            }
            INFO("field ", field);
            CHECK_FALSE(sim::AdaptiveStepSolver(dp, o).solve(decay, ColVec<1>{1.0}, {0.0, 1.0}).success);
        }
    }

    TEST_CASE_TEMPLATE("UKF weights fall back to the plain transform on a non-finite alpha or beta", T, float, double) {
        struct Still {
            [[nodiscard]] ColVec<2, T> operator()(const ColVec<2, T>& x, const ColVec<1, T>& /*u*/) const { return x; }
        };
        struct First {
            [[nodiscard]] ColVec<1, T> operator()(const ColVec<2, T>& x, const ColVec<1, T>& /*u*/) const { return ColVec<1, T>{x[0]}; }
        };
        for (int which = 0; which < 2; ++which) {
            UnscentedParams<T> p{};
            (which == 0 ? p.alpha : p.beta) = runtime_nan<T>();
            UnscentedKalmanFilter<Still> ukf{Still{}, {.Ts = static_cast<T>(0.01), .Q = Matrix<2, 2, T>::identity() * static_cast<T>(0.01), .x0 = ColVec<2, T>{T{1}, T{2}}, .unscented = p}};
            CHECK(ukf.predict({}));
            CHECK(ukf.update(First{}, ColVec<1, T>{static_cast<T>(1.5)}, static_cast<T>(0.1)));
            CHECK(damp::isfinite(ukf.state()[0]));
            CHECK(damp::isfinite(ukf.covariance()(0, 0)));
        }
    }

    TEST_CASE_TEMPLATE("encoder gate turns a non-finite acceleration limit or stale time off", T, float, double) {
        for (const T bad : non_finite<T>()) {
            EncoderGate<T> gate{{.tolerance = T{4}, .max_accel = bad, .stale_time = bad, .resync_count = 1}};
            REQUIRE(gate.check(SplitReal<T>{0, T{0}}, T{0}, static_cast<T>(1e-3)) == EncoderVerdict::Resync);
            CHECK(gate.check(SplitReal<T>{3, T{0}}, T{0}, T{1}) == EncoderVerdict::Accepted); // a second on: never stale
            CHECK_FALSE(gate.stale());
            CHECK(gate.check(SplitReal<T>{10, T{0}}, T{0}, T{1}) != EncoderVerdict::Accepted); // the window does not widen
        }
    }

    TEST_CASE_TEMPLATE("discrete and event-driven simulations refuse a non-finite step", T, float, double) {
        for (const T bad : non_finite<T>()) {
            const auto run = sim::simulate_discrete_nonlinear<1, 1, 1>(
                [](size_t /*k*/, const ColVec<1, T>& x, const ColVec<1, T>& u) { return ColVec<1, T>{x[0] + u[0]}; },
                [](const ColVec<1, T>& x) { return x; },
                [](const ColVec<1, T>& /*y*/) { return ColVec<1, T>{T{1}}; },
                ColVec<1, T>{},
                bad,
                4
            );
            CHECK_FALSE(run.success);
        }
        // A NaN fine-log step logs the events only; with no finite step the run still ends.
        const StateSpace              integ{.A = Matrix<1, 1, double>{{0.0}}, .B = Matrix<1, 1, double>{{1.0}}, .C = Matrix<1, 1, double>{{1.0}}};
        sim::PiecewiseLTI<1, 1, 1, 1> plant{.modes = {integ}};
        const auto                    schedule = [](double, std::size_t, const auto&, const auto&) { return 1.0; };
        const auto                    on_event = [](double, std::size_t m, const auto& x, const auto&, const auto& u) {
            return sim::HybridEventAction<1, 1>{.mode = m, .u = u, .x = x};
        };
        const sim::EventSimConfig<> events{.log = sim::LogPolicy::EventsOnly};
        const sim::EventSimConfig<> nan_stride{.log = sim::LogPolicy::FixedStride, .h_fine = runtime_nan<double>()};
        const auto                  a = sim::simulate_piecewise_lti(plant, 0, ColVec<1>{0.0}, ColVec<1>{1.0}, {0.0, 1.0}, schedule, on_event, events);
        const auto                  b = sim::simulate_piecewise_lti(plant, 0, ColVec<1>{0.0}, ColVec<1>{1.0}, {0.0, 1.0}, schedule, on_event, nan_stride);
        CHECK(b.t.size() == a.t.size());
    }

    TEST_CASE_TEMPLATE("a StateSpace whose Ts is NaN, infinite or negative has no time domain", T, float, double) {
        StateSpace<1, 1, 1, T> sys{.A = {{T{-1}}}, .B = {{T{1}}}, .C = {{T{1}}}};
        CHECK(sys.is_continuous());
        CHECK_FALSE(sys.is_discrete());
        sys.Ts = static_cast<T>(0.1);
        CHECK(sys.is_discrete());
        CHECK_FALSE(sys.is_continuous());
        for (const T bad : {runtime_nan<T>(), runtime_inf<T>(), -runtime_inf<T>(), static_cast<T>(-0.1)}) {
            sys.Ts = bad;
            CHECK_FALSE(sys.is_discrete());
            CHECK_FALSE(sys.is_continuous());
        }
    }

    TEST_CASE_TEMPLATE("designs and conversions refuse a plant with no time domain", T, float, double) {
        using C = damp::complex<T>;
        StateSpace<2, 1, 1, T>                       plant{.A = {{T{0}, T{1}}, {T{-2}, T{-3}}}, .B = {{T{0}}, {T{1}}}, .C = {{T{1}, T{0}}}};
        const damp::array<C, 2>                      poles{C{T{-4}, T{0}}, C{T{-5}, T{0}}};
        const damp::array<design::JordanBlock<T>, 1> blocks{{{{T{-4}, T{0}}, 2}}};
        const design::LQRWeights                     weights{.Q = Matrix<2, 2, T>::identity(), .R = {T{1}}};
        const T                                      ts{static_cast<T>(0.01)};
        REQUIRE(design::to_coeffs(plant, ts).success);
        REQUIRE(design::state_feedback(plant, weights).success);
        REQUIRE(design::place(plant, poles).success);
        REQUIRE(design::place_jordan(plant, blocks).success);
        REQUIRE(design::place_observer(plant, poles).success);
        REQUIRE(design::reduced_order_observer(plant, {C{T{-6}, T{0}}}).success);
        REQUIRE(design::balanced_realization(plant).success);
        REQUIRE(design::detail::residualize_states<1>(plant).has_value());
        auto discrete = design::state_feedback(plant, ts, weights);
        REQUIRE(discrete.success);
        for (const T bad : {runtime_nan<T>(), runtime_inf<T>()}) {
            plant.Ts = bad;
            CHECK(design::to_coeffs(plant, ts).status == DesignStatus::invalid_sample_time);
            CHECK(design::state_feedback(plant, weights).status == DesignStatus::invalid_sample_time);
            CHECK(design::place(plant, poles).status == DesignStatus::invalid_sample_time);
            CHECK(design::place_jordan(plant, blocks).status == DesignStatus::invalid_sample_time);
            CHECK(design::place_observer(plant, poles).status == DesignStatus::invalid_sample_time);
            CHECK(design::reduced_order_observer(plant, {C{T{-6}, T{0}}}).status == DesignStatus::invalid_sample_time);
            CHECK(design::balanced_realization(plant).status == DesignStatus::invalid_sample_time);
            CHECK_FALSE(design::detail::residualize_states<1>(plant).has_value());
            auto broken = discrete;
            broken.Ts = bad;
            CHECK_FALSE(LQR{broken}.valid());
        }
    }

    TEST_CASE_TEMPLATE("analysis of a model with no time domain gives no result", T, float, double) {
        StateSpace<1, 1, 1, T> sys{.A = {{T{-1}}}, .B = {{T{1}}}, .C = {{T{1}}}};
        const std::vector<T>   w{T{1}, T{10}};
        REQUIRE(analysis::bode(sys, w).points.size() == 2);
        REQUIRE(stability::is_stable(sys));
        for (const T bad : {runtime_nan<T>(), runtime_inf<T>()}) {
            sys.Ts = bad;
            CHECK(analysis::bode(sys, w).points.empty());
            CHECK(analysis::nyquist(sys, w).points.empty());
            CHECK(analysis::sigma(sys, w).points.empty());
            CHECK(analysis::loop_response(sys, w).open_loop.points.empty());
            CHECK(analysis::impedance(sys, w).points.empty());
            CHECK(analysis::impedance_from_admittance(sys, w).points.empty());
            CHECK_FALSE(analysis::dcgain(sys).has_value());
            CHECK_FALSE(analysis::norm_h2(sys).has_value());
            CHECK_FALSE(analysis::norm_hinf(sys).success);
            CHECK_FALSE(analysis::detail::discretize_on_grid(sys, std::vector<T>{T{0}, static_cast<T>(0.1)}).has_value());
            CHECK_FALSE(stability::is_stable(sys));
            CHECK_FALSE(stability::is_stabilizable(sys));
            CHECK(stability::stability_margin(sys) == T{-1});
            CHECK_FALSE(matlab::gram(sys).has_value());
            CHECK_FALSE(sim::simulate_discrete(sys, [](const ColVec<1, T>& /*y*/) { return ColVec<1, T>{}; }, ColVec<1, T>{}, 4).success);
        }
    }
}
