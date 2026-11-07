// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Law<kRes, T> of a constexpr design: the matrices constants of the type. It has to be the same law
// as the runtime with members, tick for tick, on every family, with no term for an exact zero.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/const_matrix.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

constexpr double kTs = 0.01;

constexpr StateSpace<2, 1, 1, double, 2, 1> mass_spring() {
    StateSpace<2, 1, 1, double, 2, 1> s{};
    s.A(0, 1) = 1.0;
    s.A(1, 0) = -4.0;
    s.A(1, 1) = -0.4;
    s.B(1, 0) = 1.0;
    s.C(0, 0) = 1.0;
    return s; // G = I, H = I by default
}

inline constexpr auto kPlant = *discretize(mass_spring(), kTs);
inline constexpr auto kDelayed = *discretize(mass_spring(), kTs, InputDelay{kTs}); // [x; u(k−1)]

inline constexpr Matrix<2, 2> kQ2{{1e3, 0.0}, {0.0, 1.0}};
inline constexpr Matrix<3, 3> kQ3{{1e3, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}};
inline constexpr Matrix<2, 2> kQn{{1e-4, 0.0}, {0.0, 1e-2}};

// Every family, with and without the delay register.
inline constexpr auto kLqr = design::lqr(kPlant, {.Q = kQ2, .R = 1e-2});
inline constexpr auto kLqrDelayed = design::lqr(kDelayed, {.Q = kQ3, .R = 1e-2});
inline constexpr auto kLqrTrack = design::lqr(kPlant, {.Q = kQ2, .R = 1e-2}, design::TrackedOutputs<1>{{0}});
inline constexpr auto kLqi = design::lqi(kPlant, {.Q = kQ2, .R = 1e-2, .Q_integral = 1e2});
inline constexpr auto kLqiDelayed = design::lqi(kDelayed, {.Q = kQ3, .R = 1e-2, .Q_integral = 1e2});
inline constexpr auto kLqg = design::lqg(kPlant, {.Q = kQ2, .R = 1e-2}, {.process = kQn, .measurement = 1e-4});
inline constexpr auto kLqgDelayed = design::lqg(kDelayed, {.Q = kQ3, .R = 1e-2}, {.process = kQn, .measurement = 1e-4});
inline constexpr auto kLqgTrack = design::lqg(kPlant, {.Q = kQ2, .R = 1e-2}, {.process = kQn, .measurement = 1e-4}, design::TrackedOutputs<1>{{0}});
inline constexpr auto kLqgi = design::lqgi(kPlant, {.Q = kQ2, .R = 1e-2, .Q_integral = 1e2}, {.process = kQn, .measurement = 1e-4});
inline constexpr auto kLqgiDelayed = design::lqgi(kDelayed, {.Q = kQ3, .R = 1e-2, .Q_integral = 1e2}, {.process = kQn, .measurement = 1e-4});
inline constexpr auto kOflqg = design::oflqg(kPlant, {.Q = kQ2, .R = 1e-2}, {.process = kQn, .measurement = 1e-4, .disturbance = 1e-2}, design::InputDisturbance{});
inline constexpr auto kOflqgDelayed = design::oflqg(kDelayed, {.Q = kQ3, .R = 1e-2}, {.process = kQn, .measurement = 1e-4, .disturbance = 1e-2}, design::InputDisturbance{});

/// The same design on the runtime with members.
template<const auto& Res, typename T>
[[nodiscard]] auto held() {
    return typename Law<Res, T>::Runtime{Res.template as<T>()};
}

/// A reading with some structure.
double reading(int k) {
    return (0.5 * std::sin(0.037 * k)) + (((k / 150) % 2 == 0) ? 0.3 : -0.2);
}

/// A second input: a reference, or a state to feed back.
double second(int k) {
    return (0.2 * std::cos(0.021 * k)) + (((k / 90) % 2 == 0) ? 0.1 : -0.1);
}

/// Close enough for the same arithmetic in another order: -ffast-math may associate either form.
template<typename T>
double tolerance() {
    return std::is_same_v<T, float> ? 1e-4 : 1e-11;
}

/// The whole-tick form, and the ISR's order with a clamp, on both runtimes: @p tick maps a law and k to its command.
template<typename T, typename Held, typename Fixed, typename Tick>
void check_same_commands(Held held, Fixed fixed, Tick tick) {
    double worst = 0.0;
    double scale = 0.0;
    for (int k = 0; k < 3000; ++k) {
        const double a = static_cast<double>(tick(held, k));
        const double b = static_cast<double>(tick(fixed, k));
        worst = std::max(worst, std::abs(a - b));
        scale = std::max(scale, std::abs(a));
    }
    CAPTURE(worst);
    CAPTURE(scale);
    CHECK(scale > 0.0);
    CHECK(worst <= tolerance<T>() * std::max(scale, 1.0));
}

/// Update, control, clamp, back_calculate, predict, against the design's own runtime.
template<typename T, typename Law, typename Control>
[[nodiscard]] T isr_tick(Law& law, int k, const ColVec<Law::Core::NM, T>& m, Control control) {
    law.update(m);
    const T requested = control(law)(0);
    const T applied = std::clamp(requested, T{-5}, T{5});
    law.back_calculate(ColVec<1, T>{applied});
    law.predict();
    static_cast<void>(k);
    return applied;
}

} // namespace

TEST_SUITE("ConstMatrix") {
    static constexpr Matrix<3, 4, double> kSparse{{1.0, 0.0, 0.5, 0.0}, {0.0, 0.0, 0.0, 0.0}, {-1.0, 2.0, 0.0, 3.0}};
    static constexpr Matrix<2, 2, float>  kDense{{1.0F, 2.0F}, {3.0F, 4.0F}};

    TEST_CASE("it takes no storage and has no term for an exact zero") {
        static_assert(std::is_empty_v<ConstMatrix<kSparse>>);
        static_assert(detail::nonzero_columns<kSparse, 0>.n == 2);
        static_assert(detail::nonzero_columns<kSparse, 0>.col[1] == 2);
        static_assert(detail::nonzero_columns<kSparse, 1>.n == 0);
        static_assert(detail::nonzero_columns<kSparse, 2>.n == 3);
        static_assert(detail::nonzero_columns<kDense, 1>.n == 2);
    }

    TEST_CASE("its product is the matrix's, an all-zero row included") {
        const ColVec<4> x{1.0, -2.0, 3.0, 0.5};
        const ColVec<3> y = ConstMatrix<kSparse>{} * x;
        const ColVec<3> dense = kSparse * x;
        CHECK((y - dense).norm() == 0.0);
        CHECK(y[1] == 0.0);
        constexpr ColVec<2, float> folded = ConstMatrix<kDense>{} * ColVec<2, float>{1.0F, 1.0F}; // constant-evaluated too
        static_assert(folded[0] == 3.0F && folded[1] == 7.0F);
    }

    TEST_CASE("it reads as the matrix it is") {
        constexpr ConstMatrix<kSparse> m{};
        static_assert(decltype(m)::rows() == 3 && decltype(m)::cols() == 4);
        CHECK(m(2, 3) == 3.0);
        const Matrix<3, 4, double>& held = m;
        CHECK(&held == &kSparse);
    }
}

TEST_SUITE("Constant indices past the unroll limit") {
    TEST_CASE("the matrix operators give the same answer folded (up to 64 entries) and looped") {
        const auto fill = []<size_t R, size_t C>(Matrix<R, C> m, double scale) {
            for (size_t i = 0; i < R; ++i) {
                for (size_t j = 0; j < C; ++j) {
                    m(i, j) = scale * (1.0 + static_cast<double>(i) + (0.1 * static_cast<double>(j)));
                }
            }
            return m;
        };
        const auto check = [&]<size_t R, size_t C>(Matrix<R, C> zero) {
            const Matrix<R, C> a = fill(zero, 1.0);
            const Matrix<R, C> b = fill(zero, 0.5);
            const Matrix<R, C> sum = a + b;
            const Matrix<R, C> diff = a - b;
            const Matrix<R, C> neg = -a;
            Matrix<R, C>       inc = a;
            inc += b;
            Matrix<R, C> dec = a;
            dec -= b;
            for (size_t i = 0; i < R; ++i) {
                for (size_t j = 0; j < C; ++j) {
                    CHECK(sum(i, j) == a(i, j) + b(i, j));
                    CHECK(diff(i, j) == a(i, j) - b(i, j));
                    CHECK(neg(i, j) == -a(i, j));
                    CHECK(inc(i, j) == sum(i, j));
                    CHECK(dec(i, j) == diff(i, j));
                }
            }
        };
        check(Matrix<18, 1>{}); // past 16 iterations, folded
        check(Matrix<1, 18>{});
        check(Matrix<8, 8>{}); // the largest folded
        check(Matrix<9, 8>{}); // the first looped
    }

    TEST_CASE("an aliased += reads before it writes, and constant evaluation still works") {
        Matrix<2, 2> a{{1.0, 2.0}, {3.0, 4.0}};
        a += a.transpose();
        CHECK(a(0, 1) == 5.0);
        CHECK(a(1, 0) == 5.0);
        static constexpr Matrix<2, 2> folded = Matrix<2, 2>{{1.0, 2.0}, {3.0, 4.0}} + Matrix<2, 2>{{1.0, 1.0}, {1.0, 1.0}};
        static_assert(folded(1, 1) == 5.0);
    }
}

TEST_SUITE("Law of a constexpr design") {
    TEST_CASE("no matrix takes room in the object; the same laws run as the design's runtimes") {
        static_assert(sizeof(Law<kLqgDelayed, float>) < sizeof(Law<kLqgDelayed, float>::Runtime) / 2);
        static_assert(sizeof(Law<kLqr, float>) < sizeof(Law<kLqr, float>::Runtime));
        static_assert(sizeof(Law<kOflqgDelayed, float>) < sizeof(Law<kOflqgDelayed, float>::Runtime) / 2);
        static_assert(std::is_empty_v<detail::ConstantObserver<kLqg, float, 2>::Model>);

        // The delay register's row of A, and of M, is no code at all.
        using Observer = detail::ConstantObserver<kLqgDelayed, float, 3>;
        static_assert(detail::nonzero_columns<Observer::kA, 2>.n == 0);
        static_assert(detail::nonzero_columns<Observer::kM, 2>.n == 0);
        static_assert(detail::nonzero_columns<Observer::kC, 0>.n == 1);
        static_assert(detail::nonzero_columns<Observer::kD, 0>.n == 0);
    }

    TEST_CASE("a constexpr design runs on constants; one initialized at start-up keeps members and can be replaced") {
        static_assert(detail::constant_design<kLqg>);
        static const auto kAtInit = design::lqg(kPlant, {.Q = kQ2, .R = 1e-2}, {.process = kQn, .measurement = 1e-4});
        static_assert(!detail::constant_design<kAtInit>);
        static_assert(std::is_same_v<Law<kLqg, float>::Runtime, Law<kAtInit, float>::Runtime>);
        static_assert(sizeof(Law<kAtInit, float>) == sizeof(Law<kAtInit, float>::Runtime));
        static_assert(sizeof(Law<kLqg, float>) < sizeof(Law<kAtInit, float>) / 2);

        Law<kAtInit, float> law;
        CHECK(law.K(0, 0) == static_cast<float>(kAtInit.K(0, 0)));
        const float before = law.K(0, 0);
        law = LQG{design::lqg(kPlant, {.Q = kQ2 * 10.0, .R = 1e-2}, {.process = kQn, .measurement = 1e-4}).as<float>()};
        CHECK(law.K(0, 0) > before);

        // The two forms are the same law.
        check_same_commands<double>(Law<kAtInit, double>{}, Law<kLqg, double>{}, [](auto& l, int k) { return l.step(ColVec<1, double>{reading(k)})[0]; });

        static const auto         kKalmanAtInit = design::kalman(kPlant, {.process = kQn, .measurement = 1e-4});
        Law<kKalmanAtInit, float> observer;
        observer.update(ColVec<1, float>{0.5F});
        CHECK(observer.valid());
    }

    TEST_CASE("constinit, and in float by default construction") {
        static constinit Law<kLqgDelayed, float>   lqg;
        static constinit Law<kLqr, float>          lqr;
        static constinit Law<kOflqgDelayed, float> oflqg;
        CHECK(std::isfinite(lqg.step(ColVec<1, float>{0.5F})[0]));
        CHECK(std::isfinite(lqr.step(ColVec<2, float>{0.5F, 0.0F})[0]));
        CHECK(std::isfinite(oflqg.step(ColVec<1, float>{0.5F}, ColVec<1, float>{0.1F})[0]));
    }

    TEST_CASE("what a law exposes reads the design's numbers") {
        Law<kLqgDelayed, double> law;
        for (size_t j = 0; j < 3; ++j) {
            CHECK(law.K(0, j) == kLqgDelayed.K(0, j));
            CHECK(law.getK()(0, j) == kLqgDelayed.K(0, j));
            CHECK(law.estimator.gain()(j, 0) == kLqgDelayed.kalman.M(j, 0));
        }
        CHECK(law.estimator.model().A(0, 1) == kLqgDelayed.kalman.sys.A(0, 1));
        CHECK(law.estimator.valid());
    }

    TEST_CASE("LQR: the same commands") {
        const auto tick = [](auto& law, int k) { return law.step(ColVec<2, double>{reading(k), second(k)})[0]; };
        check_same_commands<double>(held<kLqr, double>(), Law<kLqr, double>{}, tick);
        const auto delayed = [](auto& law, int k) { return law.step(ColVec<3, double>{reading(k), second(k), reading(k - 1)})[0]; };
        check_same_commands<double>(held<kLqrDelayed, double>(), Law<kLqrDelayed, double>{}, delayed);
        const auto tracked = [](auto& law, int k) { return law.step(ColVec<1, double>{second(k)}, ColVec<2, double>{reading(k), 0.1})[0]; };
        check_same_commands<double>(held<kLqrTrack, double>(), Law<kLqrTrack, double>{}, tracked);
    }

    TEST_CASE("LQI: the same commands, the integrators held against a clamp") {
        const auto tick = [](auto& law, int k) {
            using L = std::remove_cvref_t<decltype(law)>;
            return isr_tick<double>(law, k, ColVec<2, double>{reading(k), second(k)}, [&](L& l) { return l.control(ColVec<1, double>{second(k)}); });
        };
        check_same_commands<double>(held<kLqi, double>(), Law<kLqi, double>{}, tick);
        const auto delayed = [](auto& law, int k) {
            using L = std::remove_cvref_t<decltype(law)>;
            return isr_tick<double>(law, k, ColVec<3, double>{reading(k), second(k), reading(k - 1)}, [&](L& l) { return l.control(ColVec<1, double>{second(k)}); });
        };
        check_same_commands<double>(held<kLqiDelayed, double>(), Law<kLqiDelayed, double>{}, delayed);
    }

    TEST_CASE("LQG: the same commands, the filter on the input applied") {
        const auto whole = [](auto& law, int k) { return law.step(ColVec<1, double>{reading(k)})[0]; };
        check_same_commands<double>(held<kLqg, double>(), Law<kLqg, double>{}, whole);
        check_same_commands<double>(held<kLqgDelayed, double>(), Law<kLqgDelayed, double>{}, whole);
        const auto isr = [](auto& law, int k) {
            using L = std::remove_cvref_t<decltype(law)>;
            return isr_tick<double>(law, k, ColVec<1, double>{reading(k)}, [](L& l) { return l.control(); });
        };
        check_same_commands<double>(held<kLqgDelayed, double>(), Law<kLqgDelayed, double>{}, isr);
        const auto feedforward = [](auto& law, int k) { return law.step(ColVec<3, double>{second(k), 0.0, second(k - 1)}, ColVec<1, double>{0.2}, ColVec<1, double>{reading(k)})[0]; };
        check_same_commands<double>(held<kLqgDelayed, double>(), Law<kLqgDelayed, double>{}, feedforward);
        const auto tracked = [](auto& law, int k) { return law.step(ColVec<1, double>{second(k)}, ColVec<1, double>{reading(k)})[0]; };
        check_same_commands<double>(held<kLqgTrack, double>(), Law<kLqgTrack, double>{}, tracked);
    }

    TEST_CASE("LQGI: the same commands, the integrators held against a clamp") {
        const auto isr = [](auto& law, int k) {
            using L = std::remove_cvref_t<decltype(law)>;
            return isr_tick<double>(law, k, ColVec<1, double>{reading(k)}, [&](L& l) { return l.control(ColVec<1, double>{second(k)}); });
        };
        check_same_commands<double>(held<kLqgi, double>(), Law<kLqgi, double>{}, isr);
        check_same_commands<double>(held<kLqgiDelayed, double>(), Law<kLqgiDelayed, double>{}, isr);
    }

    TEST_CASE("OFLQG: the same commands and disturbance estimate, the anchor moving and a reference point") {
        const auto isr = [](auto& law, int k) {
            using L = std::remove_cvref_t<decltype(law)>;
            return isr_tick<double>(law, k, ColVec<1, double>{reading(k)}, [&](L& l) { return l.control(ColVec<1, double>{second(k)}); });
        };
        check_same_commands<double>(held<kOflqg, double>(), Law<kOflqg, double>{}, isr);
        check_same_commands<double>(held<kOflqgDelayed, double>(), Law<kOflqgDelayed, double>{}, isr);
        const auto point = [](auto& law, int k) {
            return law.step(ColVec<1, double>{second(k)}, ColVec<3, double>{0.1 * second(k), 0.0, 0.05}, ColVec<1, double>{0.05}, ColVec<1, double>{reading(k)})[0];
        };
        check_same_commands<double>(held<kOflqgDelayed, double>(), Law<kOflqgDelayed, double>{}, point);

        auto                       held = ::held<kOflqgDelayed, double>();
        Law<kOflqgDelayed, double> fixed;
        for (int k = 0; k < 500; ++k) {
            static_cast<void>(held.step(ColVec<1, double>{0.2}, ColVec<1, double>{reading(k)}));
            static_cast<void>(fixed.step(ColVec<1, double>{0.2}, ColVec<1, double>{reading(k)}));
        }
        CHECK(std::abs(held.disturbance_estimate()[0] - fixed.disturbance_estimate()[0]) < 1e-10);
        CHECK(std::abs(held.state_estimate()[0] - fixed.state_estimate()[0]) < 1e-10);
    }

    TEST_CASE("an observer design: the same estimate as Luenberger, a Kalman gain and a placed one") {
        static constexpr auto kKalman = design::kalman(kDelayed, {.process = kQn, .measurement = 1e-4});
        static constexpr auto kPlaced = design::luenberger(kPlant, {0.9, 0.92});
        static constexpr auto kAnchored = design::kalman(kDelayed, {.process = kQn, .measurement = 1e-4, .disturbance = 1e-2}, design::InputDisturbance{});
        static_assert(sizeof(Law<kKalman, float>) < sizeof(Law<kKalman, float>::Runtime) / 2);

        const auto run = [](auto held, auto fixed) {
            double worst = 0.0;
            for (int k = 0; k < 2000; ++k) {
                const ColVec<1> u{{second(k)}};
                const ColVec<1> y{{reading(k)}};
                held.predict(u);
                fixed.predict(u);
                held.update(y);
                fixed.update(y);
                worst = std::max(worst, (held.state() - fixed.state()).norm());
                worst = std::max(worst, (held.innovation() - fixed.innovation()).norm());
            }
            CHECK(worst < 1e-10);
        };
        run(Luenberger{kKalman}, Law<kKalman, double>{});
        run(Luenberger{kPlaced}, Law<kPlaced, double>{});

        // The anchor: the deviation shifts with it, the reference point's motion carried.
        Luenberger             held{kAnchored};
        Law<kAnchored, double> fixed;
        double                 worst = 0.0;
        for (int k = 0; k < 1000; ++k) {
            const ColVec<3> x_t{0.1 * second(k), 0.0, 0.2};
            const ColVec<1> u_t{{0.2}};
            const ColVec<3> x_ref{0.0, 0.05 * second(k), 0.0};
            const ColVec<1> u_ref{{0.01}};
            held.anchor(x_t, u_t, x_ref, u_ref);
            fixed.anchor(x_t, u_t, x_ref, u_ref);
            held.predict(ColVec<1>{{second(k)}});
            fixed.predict(ColVec<1>{{second(k)}});
            held.update(ColVec<1>{{reading(k)}});
            fixed.update(ColVec<1>{{reading(k)}});
            worst = std::max(worst, (held.state() - fixed.state()).norm());
            worst = std::max(worst, std::abs(held.disturbance_estimate()[0] - fixed.disturbance_estimate()[0]));
        }
        CHECK(worst < 1e-10);
        static constinit Law<kKalman, float> deployed;
        deployed.predict(ColVec<1, float>{0.1F});
        deployed.update(ColVec<1, float>{0.5F});
        CHECK(deployed.valid());
        CHECK(deployed.model().A(0, 1) == static_cast<float>(kKalman.sys.A(0, 1)));
    }

    TEST_CASE("deployed in float: the same commands to float's precision") {
        const auto lqg = [](auto& law, int k) { return law.step(ColVec<1, float>{static_cast<float>(reading(k))})[0]; };
        check_same_commands<float>(held<kLqgDelayed, float>(), Law<kLqgDelayed, float>{}, lqg);
        const auto oflqg = [](auto& law, int k) { return law.step(ColVec<1, float>{0.2F}, ColVec<1, float>{static_cast<float>(reading(k))})[0]; };
        check_same_commands<float>(held<kOflqgDelayed, float>(), Law<kOflqgDelayed, float>{}, oflqg);
    }
}
