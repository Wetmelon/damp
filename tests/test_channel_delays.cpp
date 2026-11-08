// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Delay per channel, and on the measurement. discretize(plant, Ts, OutputDelay<N>{θ}) holds a late
// reading in N output registers, and InputDelay / OutputDelay take a ColVec with one delay per
// channel. The oracle is the continuous plant integrated on a fine grid with each command arriving
// θ_j late and each reading taken θ_i before it arrives; the discrete model has to read the same.
// The designs then run on that model like on any other plant.

#include <cmath>
#include <cstddef>
#include <vector>

#include "damp/controllers/lqg.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/ofmpc.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

constexpr double kTs = 0.1;
constexpr int    kGrid = 200; // fine steps per period: every delay below is a whole number of them
constexpr double kDt = kTs / kGrid;

/// Three states, two inputs, two outputs, a feedthrough on the second output and noise on every state.
constexpr StateSpace<3, 2, 2, double, 3, 2> two_by_two() {
    StateSpace<3, 2, 2, double, 3, 2> s{};
    s.A = Matrix<3, 3>{{-1.0, 1.0, 0.0}, {0.0, -2.0, 1.0}, {0.0, 0.0, -0.5}};
    s.B = Matrix<3, 2>{{1.0, 0.0}, {0.0, 0.5}, {0.3, 1.0}};
    s.C = Matrix<2, 3>{{1.0, 0.0, 0.0}, {0.0, 1.0, 1.0}};
    s.D = Matrix<2, 2>{{0.0, 0.0}, {0.2, -0.1}};
    return s; // G = I, H = I
}

/// Lightly damped mass-spring (ω = 2 rad/s, ζ = 0.1): position, speed; force in; position out.
constexpr StateSpace<2, 1, 1, double, 2, 1> mass_spring() {
    StateSpace<2, 1, 1, double, 2, 1> s{};
    s.A(0, 1) = 1.0;
    s.A(1, 0) = -4.0;
    s.A(1, 1) = -0.4;
    s.B(1, 0) = 1.0;
    s.C(0, 0) = 1.0;
    return s;
}

int floor_div(int a, int b) {
    return (a >= 0) ? a / b : -((-a + b - 1) / b);
}

double command(int k, size_t j) {
    return (k < 0) ? 0.0 : std::sin((0.7 * k) + static_cast<double>(j)) + ((k / 9) % 2 == 0 ? 0.4 : -0.4);
}

double held_noise(int k, size_t w) {
    return (k < 0) ? 0.0 : 0.3 * std::cos((1.3 * k) + (2.0 * static_cast<double>(w)));
}

/**
 * The readings of the continuous plant over K periods: input j arrives in_steps[j] fine steps after
 * it is sent, and the reading of output i that arrives at sample k was taken out_steps[i] fine
 * steps earlier. RK4 on the fine grid; the inputs and the noise only switch on it.
 */
template<size_t NX, size_t NU, size_t NY, size_t NW, size_t NV>
std::vector<ColVec<NY>> continuous_readings(
    const StateSpace<NX, NU, NY, double, NW, NV>& plant,
    const damp::array<int, NU>&                   in_steps,
    const damp::array<int, NY>&                   out_steps,
    int                                           K
) {
    const auto u_at = [&](int n) { // the input the plant holds over fine step n (and from its start)
        ColVec<NU> u{};
        for (size_t j = 0; j < NU; ++j) {
            u[j] = command(floor_div(n - in_steps[j], kGrid), j);
        }
        return u;
    };
    const auto w_at = [&](int n) {
        ColVec<NW> w{};
        for (size_t i = 0; i < NW; ++i) {
            w[i] = held_noise(floor_div(n, kGrid), i);
        }
        return w;
    };
    std::vector<ColVec<NX>> x(static_cast<size_t>(K * kGrid) + 1);
    for (int n = 0; n < K * kGrid; ++n) {
        const ColVec<NX> drive = ColVec<NX>((plant.B * u_at(n)) + (plant.G * w_at(n)));
        const auto       f = [&](const ColVec<NX>& v) { return ColVec<NX>((plant.A * v) + drive); };
        const ColVec<NX> x0 = x[static_cast<size_t>(n)];
        const ColVec<NX> k1 = f(x0);
        const ColVec<NX> k2 = f(ColVec<NX>(x0 + (k1 * (0.5 * kDt))));
        const ColVec<NX> k3 = f(ColVec<NX>(x0 + (k2 * (0.5 * kDt))));
        const ColVec<NX> k4 = f(ColVec<NX>(x0 + (k3 * kDt)));
        x[static_cast<size_t>(n) + 1] = ColVec<NX>(x0 + ((k1 + (k2 * 2.0) + (k3 * 2.0) + k4) * (kDt / 6.0)));
    }
    std::vector<ColVec<NY>> y(static_cast<size_t>(K));
    for (int k = 0; k < K; ++k) {
        for (size_t i = 0; i < NY; ++i) {
            const int m = (k * kGrid) - out_steps[i];
            if (m < 0) {
                continue; // before the start, at rest
            }
            const ColVec<NX> xm = x[static_cast<size_t>(m)];
            const ColVec<NU> um = u_at(m);
            double           v = 0.0;
            for (size_t l = 0; l < NX; ++l) {
                v += plant.C(i, l) * xm[l];
            }
            for (size_t j = 0; j < NU; ++j) {
                v += plant.D(i, j) * um[j];
            }
            y[static_cast<size_t>(k)][i] = v;
        }
    }
    return y;
}

/// The readings of a discrete model driven by the same commands and held noise, from rest.
template<size_t NS, size_t NU, size_t NY, size_t NW, size_t NV>
std::vector<ColVec<NY>> model_readings(const StateSpace<NS, NU, NY, double, NW, NV>& model, int K) {
    std::vector<ColVec<NY>> y(static_cast<size_t>(K));
    ColVec<NS>              x{};
    for (int k = 0; k < K; ++k) {
        ColVec<NU> u{};
        for (size_t j = 0; j < NU; ++j) {
            u[j] = command(k, j);
        }
        ColVec<NW> w{};
        for (size_t i = 0; i < NW; ++i) {
            w[i] = held_noise(k, i);
        }
        y[static_cast<size_t>(k)] = ColVec<NY>((model.C * x) + (model.D * u));
        x = ColVec<NS>((model.A * x) + (model.B * u) + (model.G * w));
    }
    return y;
}

/// Largest gap between two runs of readings, relative to the largest reading.
template<size_t NY>
double relative_gap(const std::vector<ColVec<NY>>& a, const std::vector<ColVec<NY>>& b) {
    double gap = 0.0;
    double scale = 1e-300;
    for (size_t k = 0; k < a.size(); ++k) {
        for (size_t i = 0; i < NY; ++i) {
            gap = std::fmax(gap, std::abs(a[k][i] - b[k][i]));
            scale = std::fmax(scale, std::abs(b[k][i]));
        }
    }
    return gap / scale;
}

/// The model against the continuous oracle: delays in fine steps per channel.
template<size_t NS, size_t NW, size_t NV>
void check_against_continuous(const StateSpace<NS, 2, 2, double, NW, NV>& model, damp::array<int, 2> in_steps, damp::array<int, 2> out_steps) {
    constexpr int K = 60;
    const auto    truth = continuous_readings(two_by_two(), in_steps, out_steps, K);
    CAPTURE(in_steps[0]);
    CAPTURE(in_steps[1]);
    CAPTURE(out_steps[0]);
    CAPTURE(out_steps[1]);
    CHECK(relative_gap(model_readings(model, K), truth) < 1e-9);
}

template<size_t R, size_t C>
bool same(const Matrix<R, C>& a, const Matrix<R, C>& b) {
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            if (a(i, j) != b(i, j)) {
                return false;
            }
        }
    }
    return true;
}

template<size_t NS, size_t NU, size_t NY, size_t NW, size_t NV>
bool same_model(const StateSpace<NS, NU, NY, double, NW, NV>& a, const StateSpace<NS, NU, NY, double, NW, NV>& b) {
    return same(a.A, b.A) && same(a.B, b.B) && same(a.C, b.C) && same(a.D, b.D) && same(a.G, b.G) && same(a.H, b.H) && a.Ts == b.Ts;
}

/// The true plant, the delayed model itself, driven by u plus a load at the input.
template<size_t NS>
struct Truth {
    StateSpace<NS, 1, 1, double, 2, 1> m{};
    ColVec<NS>                         x{};
    [[nodiscard]] ColVec<1>            y() const { return ColVec<1>(m.C * x); }
    void                               step(double u, double load = 0.0) { x = ColVec<NS>((m.A * x) + (m.B * ColVec<1>{{u + load}})); }
};

/// Largest |value| over the last quarter of a run.
struct Tail {
    double worst = 0.0;
    void   add(int k, int n, double v) {
        if (4 * k >= 3 * n) {
            worst = std::fmax(worst, std::abs(v));
        }
    }
};

double u_at(int k) {
    return (k < 0) ? 0.0 : std::sin(0.037 * k) + ((k / 150) % 2 == 0 ? 0.5 : -0.5);
}

constexpr double              kTsSpring = 0.01;
inline constexpr auto         kLate = *discretize(mass_spring(), kTsSpring, OutputDelay<2>{1.5 * kTsSpring}); // [x; z(k−1); z(k−2)]
inline constexpr Matrix<2, 2> kQn{{1e-6, 0.0}, {0.0, 1e-4}};
inline constexpr auto         kLateKalman = design::kalman(kLate, {.process = kQn, .measurement = 1e-6});

} // namespace

TEST_SUITE("Delay per channel and on the measurement") {
    TEST_CASE("the delay's type: a scalar for every channel, a ColVec for one per channel") {
        static_assert(damp::is_same_v<decltype(InputDelay{0.4}), InputDelay<1, double>>);
        static_assert(damp::is_same_v<decltype(OutputDelay{0.4}), OutputDelay<1, double>>);
        static_assert(damp::is_same_v<decltype(OutputDelay{ColVec<2>{0.1, 0.2}}), OutputDelay<1, ColVec<2>>>);
        using Both = decltype(discretize(two_by_two(), kTs, InputDelay<2>{0.1}, OutputDelay<3>{0.1}));
        static_assert(damp::is_same_v<Both, DiscretizationResult<3 + (2 * 2) + (3 * 2), 2, 2, double, 3, 2>>);
        using Late = decltype(discretize(two_by_two(), kTs, OutputDelay<2>{0.1}));
        static_assert(damp::is_same_v<Late, DiscretizationResult<3 + (2 * 2), 2, 2, double, 3, 2>>);
    }

    TEST_CASE("an output delay reads what the continuous plant read that long ago: whole samples and fractions") {
        check_against_continuous(*discretize(two_by_two(), kTs, OutputDelay<2>{2.0 * kTs}), {0, 0}, {400, 400});
        check_against_continuous(*discretize(two_by_two(), kTs, OutputDelay<2>{1.5 * kTs}), {0, 0}, {300, 300});
        check_against_continuous(*discretize(two_by_two(), kTs, OutputDelay{0.25 * kTs}), {0, 0}, {50, 50});
        check_against_continuous(*discretize(two_by_two(), kTs, OutputDelay<2, ColVec<2>>{{0.0, 1.2 * kTs}}), {0, 0}, {0, 240});
        check_against_continuous(*discretize(two_by_two(), kTs, OutputDelay<2, ColVec<2>>{{0.45 * kTs, 2.0 * kTs}}), {0, 0}, {90, 400});
    }

    TEST_CASE("one delay per input: each column moves by its own delay") {
        check_against_continuous(*discretize(two_by_two(), kTs, InputDelay<2, ColVec<2>>{{1.0 * kTs, 2.0 * kTs}}), {200, 400}, {0, 0});
        check_against_continuous(*discretize(two_by_two(), kTs, InputDelay<2, ColVec<2>>{{0.3 * kTs, 1.6 * kTs}}), {60, 320}, {0, 0});
        check_against_continuous(*discretize(two_by_two(), kTs, InputDelay<1, ColVec<2>>{{0.0, 0.7 * kTs}}), {0, 140}, {0, 0});
    }

    TEST_CASE("input and output delays together: each reading sees the inputs the plant held when it was taken") {
        check_against_continuous(*discretize(two_by_two(), kTs, InputDelay<2, ColVec<2>>{{0.3 * kTs, 1.6 * kTs}}, OutputDelay<2, ColVec<2>>{{0.45 * kTs, 1.2 * kTs}}), {60, 320}, {90, 240});
        check_against_continuous(*discretize(two_by_two(), kTs, InputDelay<2>{0.5 * kTs}, OutputDelay<2>{0.5 * kTs}), {100, 100}, {100, 100}); // the reading at the switch: the new input
        check_against_continuous(*discretize(two_by_two(), kTs, InputDelay<2>{1.0 * kTs}, OutputDelay<2, ColVec<2>>{{0.0, 0.8 * kTs}}), {200, 200}, {0, 160});
    }

    TEST_CASE("a delay per input is that input's single-input model; one value for every channel is the same value per channel") {
        const auto both = *discretize(two_by_two(), kTs, InputDelay<2, ColVec<2>>{{0.3 * kTs, 1.6 * kTs}});
        for (size_t j = 0; j < 2; ++j) {
            CAPTURE(j);
            StateSpace<3, 1, 2, double, 3, 2> alone{};
            alone.A = two_by_two().A;
            alone.C = two_by_two().C;
            for (size_t r = 0; r < 3; ++r) {
                alone.B(r, 0) = two_by_two().B(r, j);
            }
            for (size_t i = 0; i < 2; ++i) {
                alone.D(i, 0) = two_by_two().D(i, j);
            }
            const auto single = *discretize(alone, kTs, InputDelay<2>{(j == 0) ? 0.3 * kTs : 1.6 * kTs});
            ColVec<7>  xb{};
            ColVec<5>  xs{};
            double     gap = 0.0;
            for (int k = 0; k < 40; ++k) {
                ColVec<2> u{};
                u[j] = command(k, j);
                const ColVec<1> uj{{u[j]}};
                const ColVec<2> yb = (both.C * xb) + (both.D * u);
                const ColVec<2> ys = (single.C * xs) + (single.D * uj);
                gap = std::fmax(gap, (yb - ys).norm());
                xb = (both.A * xb) + (both.B * u);
                xs = (single.A * xs) + (single.B * uj);
            }
            CHECK(gap < 1e-13);
        }
        CHECK(same_model(*discretize(two_by_two(), kTs, InputDelay<2>{1.3 * kTs}), *discretize(two_by_two(), kTs, InputDelay<2, ColVec<2>>{{1.3 * kTs, 1.3 * kTs}})));
        CHECK(same_model(*discretize(two_by_two(), kTs, OutputDelay<2>{1.3 * kTs}), *discretize(two_by_two(), kTs, OutputDelay<2, ColVec<2>>{{1.3 * kTs, 1.3 * kTs}})));
    }

    TEST_CASE("whole samples: any method, and a model already discrete") {
        const auto discrete = *discretize(two_by_two(), kTs);
        CHECK(same_model(*discretize(discrete, kTs, OutputDelay<2>{2.0 * kTs}), *discretize(two_by_two(), kTs, OutputDelay<2>{2.0 * kTs})));
        CHECK(same_model(*discretize(discrete, kTs, InputDelay<2, ColVec<2>>{{kTs, 2.0 * kTs}}, OutputDelay{kTs}), *discretize(two_by_two(), kTs, InputDelay<2, ColVec<2>>{{kTs, 2.0 * kTs}}, OutputDelay{kTs})));
        CHECK(discretize(two_by_two(), kTs, OutputDelay<2>{2.0 * kTs}, DiscretizationMethod::Tustin).success);
        const auto none = *discretize(two_by_two(), kTs, OutputDelay<2>{0.0});
        CHECK(same(Matrix<3, 3>(none.A.block<3, 3>(0, 0)), discrete.A));
        CHECK(same(Matrix<2, 3>(none.C.block<2, 3>(0, 0)), discrete.C)); // no delay: the plant's own reading, the registers unread
        CHECK(same(none.D, discrete.D));
    }

    TEST_CASE("refusals: a delay that is negative, not finite or past N·Ts, and a fraction where it cannot be exact") {
        using design::DesignStatus;
        const auto discrete = *discretize(two_by_two(), kTs);
        CHECK(discretize(two_by_two(), kTs, OutputDelay<2>{-0.01}).status == DesignStatus::invalid_parameter);
        CHECK(discretize(two_by_two(), kTs, OutputDelay<2>{std::nan("")}).status == DesignStatus::invalid_parameter);
        CHECK(discretize(two_by_two(), kTs, OutputDelay<2>{HUGE_VAL}).status == DesignStatus::invalid_parameter);
        CHECK(discretize(two_by_two(), kTs, OutputDelay<2>{2.0 * kTs + 1e-6}).status == DesignStatus::invalid_parameter);
        CHECK(discretize(two_by_two(), kTs, OutputDelay<2, ColVec<2>>{{0.0, 0.5}}).status == DesignStatus::invalid_parameter);
        CHECK(discretize(two_by_two(), kTs, InputDelay<2, ColVec<2>>{{0.05, -1e-3}}).status == DesignStatus::invalid_parameter);
        CHECK(discretize(two_by_two(), kTs, InputDelay<2>{0.1}, OutputDelay<1>{0.3}).status == DesignStatus::invalid_parameter);
        CHECK(discretize(discrete, kTs, OutputDelay{0.5 * kTs}).status == DesignStatus::invalid_sample_time);
        CHECK(discretize(discrete, kTs, InputDelay<2, ColVec<2>>{{kTs, 0.5 * kTs}}).status == DesignStatus::invalid_sample_time);
        CHECK(discretize(two_by_two(), kTs, OutputDelay{0.5 * kTs}, DiscretizationMethod::Tustin).status == DesignStatus::invalid_parameter);
        CHECK(discretize(two_by_two(), 0.0, OutputDelay{0.0}).status == DesignStatus::invalid_sample_time);
        CHECK(discretize(two_by_two(), kTs, OutputDelay<2>{2.0 * kTs}).success); // N·Ts exactly fits
    }

    TEST_CASE("a Kalman filter on the late reading estimates the current state; the delay-blind filter lags") {
        constexpr int n = 3000;
        const auto    undelayed = design::kalman(*discretize(mass_spring(), kTsSpring), {.process = kQn, .measurement = 1e-6});
        REQUIRE(kLateKalman.success);
        REQUIRE(undelayed.success);
        const auto run = [&](auto& est) {
            Truth<4> t{kLate, ColVec<4>{1.0, 0.0, 0.0, 0.0}};
            Tail     e;
            for (int k = 0; k < n; ++k) {
                static_cast<void>(est.update(t.y()));
                const auto x_hat = [&] {
                    if constexpr (requires { est.state_estimate(); }) {
                        return est.state_estimate();
                    } else {
                        return est.state();
                    }
                }();
                e.add(k, n, std::abs(x_hat[0] - t.x[0]) + std::abs(x_hat[1] - t.x[1]));
                est.predict(ColVec<1>{{u_at(k)}});
                t.step(u_at(k));
            }
            return e.worst;
        };
        FixedGainObserver fixed{kLateKalman};
        CHECK(run(fixed) < 1e-6);
        KalmanFilter recursive{kLate, {.Q = kQn, .R = 1e-6}};
        CHECK(run(recursive) < 1e-6);
        FixedGainObserver blind{undelayed};
        CHECK(run(blind) > 1e-3);
    }

    TEST_CASE("LQG, offset-free LQG and MPC on the late reading hold the position, under a load for the offset-free ones") {
        constexpr size_t     NS = 4;
        const Matrix<NS, NS> Q = blkdiag(Matrix<2, 2>{{1e3, 0.0}, {0.0, 0.0}}, Matrix<2, 2>{});
        const auto           lqg = design::lqg(kLate, {.Q = Q, .R = Matrix<1, 1>{{1e-4}}}, design::KalmanNoise<2, 1>{.process = kQn, .measurement = 1e-6});
        REQUIRE(lqg.success);
        LQG      regulating{lqg};
        Truth<4> r{kLate, ColVec<4>{1.0, 0.0, 0.0, 0.0}};
        Tail     tr;
        for (int k = 0; k < 2000; ++k) {
            r.step(regulating.step(r.y())[0]);
            tr.add(k, 2000, r.x[0]);
        }
        CHECK(tr.worst < 1e-6);

        const design::KalmanNoise noise{Matrix<2, 2>{{1e-4, 0.0}, {0.0, 1e-2}}, Matrix<1, 1>{{1e-6}}, Matrix<1, 1>{{1e-2}}};
        const auto                oflqg = design::oflqg(kLate, {.Q = Q, .R = Matrix<1, 1>{{1e-4}}}, noise, design::InputDisturbance{});
        REQUIRE(oflqg.success);
        OFLQG    holding{oflqg};
        Truth<4> a{kLate};
        Tail     ta;
        for (int k = 0; k < 2000; ++k) {
            a.step(holding.step(ColVec<1>{{0.3}}, a.y())[0], (k > 500) ? 0.5 : 0.0);
            ta.add(k, 2000, a.x[0] - 0.3);
        }
        CHECK(ta.worst < 1e-4);

        design::MPCWeights<1, 1> w{};
        w.Qy = 1e3;
        w.Rdu = 1e-4;
        const auto ofmpc = design::ofmpc<25, 5>(kLate, w, {}, noise, design::InputDisturbance{});
        REQUIRE(ofmpc.success);
        OFMPC    planner{ofmpc};
        Truth<4> b{kLate};
        Tail     tb;
        for (int k = 0; k < 2000; ++k) {
            b.step(planner.step(ColVec<1>{{0.2}}, b.y())[0], (k > 500) ? 0.5 : 0.0);
            tb.add(k, 2000, b.x[0] - 0.2);
        }
        CHECK(tb.worst < 1e-4);
    }

    TEST_CASE("MPC counts an output delay's dead time: a horizon the move cannot reach is refused, the suggestion covers it") {
        constexpr StateSpace<1, 1, 1> lag{.A = {{0.9}}, .B = {{0.1}}, .C = {{1.0}}, .Ts = kTs};
        constexpr auto                model = discretize(lag, kTs, OutputDelay<5>{5.0 * kTs}).value();
        design::MPCWeights<1, 1>      w{};
        w.Rdu = 0.05;
        CHECK(design::mpc<5, 2>(model, w, {}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::mpc<6, 2>(model, w, {}).success);
        constexpr auto h = design::suggest_mpc_horizon(model);
        static_assert(h.success);
        CHECK(h.prediction_horizon == 43); // 4τ = 3.80 s and the 0.5 s delay, as for an input delay
        CHECK(design::mpc<h.prediction_horizon, h.control_horizon>(model, w, {}).success);
    }

    TEST_CASE("a placed observer on the output registers: one pole per output pinned at 0, the error gone") {
        const auto                      model = discretize(mass_spring(), kTsSpring, OutputDelay<2>{2.0 * kTsSpring}).value();
        damp::array<complex<double>, 4> poles{complex<double>{0.9, 0.0}, complex<double>{0.92, 0.0}, complex<double>{0.0, 0.0}, complex<double>{1e-3, 0.0}};
        const auto                      placed = design::place_observer(model, poles);
        REQUIRE(placed.success);
        CHECK(placed.is_stable());
        FixedGainObserver observer{placed};
        Truth<4>          t{model, ColVec<4>{1.0, 0.0, 0.0, 0.0}};
        Tail              e;
        for (int k = 0; k < 600; ++k) {
            observer.update(t.y());
            e.add(k, 600, std::abs(observer.state_estimate()[0] - t.x[0]) + std::abs(observer.state_estimate()[1] - t.x[1]));
            observer.predict(ColVec<1>{{u_at(k)}});
            t.step(u_at(k));
        }
        CHECK(e.worst < 1e-9);
    }

    TEST_CASE("a constexpr design on the late reading deploys as Law on its constants, in float") {
        static constinit Law<kLateKalman, float> deployed;
        FixedGainObserver                        held{kLateKalman.as<float>()};
        double                                   gap = 0.0;
        double                                   scale = 0.0;
        Truth<4>                                 t{kLate, ColVec<4>{0.5, 0.0, 0.0, 0.0}};
        for (int k = 0; k < 400; ++k) {
            const ColVec<1, float> y{{static_cast<float>(t.y()[0])}};
            deployed.update(y);
            held.update(y);
            for (size_t i = 0; i < 4; ++i) {
                gap = std::fmax(gap, static_cast<double>(std::abs(deployed.state_estimate()[i] - held.state_estimate()[i])));
                scale = std::fmax(scale, static_cast<double>(std::abs(held.state_estimate()[i])));
            }
            const ColVec<1, float> u{{static_cast<float>(u_at(k))}};
            deployed.predict(u);
            held.predict(u);
            t.step(u_at(k));
        }
        CHECK(gap <= 1e-5 * scale);
        static_assert(sizeof(Law<kLateKalman, float>) < sizeof(FixedGainObserver<4, 1, 1, float>));
    }

    TEST_CASE("a float plant takes the delay in double and gives the double model to float precision") {
        const auto single = *discretize(mass_spring().as<float>(), 0.01F, OutputDelay<2>{0.015});
        const auto twice = discretize(mass_spring(), 0.01, OutputDelay<2>{0.015})->as<float>();
        float      gap = 0.0F;
        const auto widen = [&](const auto& a, const auto& b) {
            for (size_t i = 0; i < a.rows(); ++i) {
                for (size_t j = 0; j < a.cols(); ++j) {
                    gap = std::fmax(gap, std::abs(a(i, j) - b(i, j)));
                }
            }
        };
        widen(single.A, twice.A);
        widen(single.B, twice.B);
        widen(single.C, twice.C);
        widen(single.G, twice.G);
        CHECK(gap < 1e-5F);
    }
}
