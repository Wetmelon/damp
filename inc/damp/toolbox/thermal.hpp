// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file thermal.hpp
 * @brief General thermal-network state-space models and Kalman observer
 *
 * Absolute-temperature RC ladders with ambient (or coolant) as a real input,
 * MIMO multi-node outputs, and a thin ThermalKalmanObserver that wraps
 * KalmanFilter for predict/update with @f$ u = [P,\,T_{\mathrm{amb}}] @f$
 * and measured NTC node(s).
 *
 * FET-specific loss models and open-loop junction estimators remain in
 * @ref damp/motor/thermal.hpp (design::foster_thermal_ss,
 * design::cauer_thermal_ss, JunctionEstimator).
 *
 * Example: absolute Cauer ladder + steady-state Kalman design
 * @code
 * #include "damp/toolbox/thermal.hpp"
 * #include "damp/systems/discretization.hpp"
 * #include "damp/estimation/kalman.hpp"
 *
 * using namespace damp;
 *
 * constexpr auto cont = design::cauer_thermal_ss_ambient<2, 1>( // measure case node
 *     damp::array<double, 2>{0.4, 0.6}, damp::array<double, 2>{0.02, 0.1});
 * constexpr auto disc = *discretize(cont, 1e-3, DiscretizationMethod::ZOH);
 * // G←I so Q is per-sample temperature variance (see thermal_kalman_plant).
 * constexpr auto kf   = design::kalman(
 *     design::thermal_kalman_plant(disc),
 *     Matrix<2, 2>::identity() * 0.1, Matrix<1, 1>{{1.0}});
 * static_assert(kf.success);
 *
 * ThermalKalmanObserver obs{kf.as<float>()};
 * // tick: obs.predict(P_w, T_amb); obs.update(t_case_ntc);
 * // float tj = obs.junction();
 * @endcode
 *
 * @see design::cauer_thermal_ss for the rise-above-reference (NU = 1) form
 * @see "Heat Transfer" (Incropera & DeWitt) for RC thermal-network modelling
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp::design {

namespace detail {

/**
 * @brief Fill the Cauer ladder A matrix and power/ambient B columns.
 *
 * Shared by the single-output and full-output ambient designers. Returns false
 * (and leaves @p A / @p B unchanged) if any R or C is non-positive.
 */
template<std::size_t N, typename T>
[[nodiscard]] constexpr bool fill_cauer_ambient_ab(
    const damp::array<T, N>& R,
    const damp::array<T, N>& C,
    Matrix<N, N, T>&         A,
    Matrix<N, 2, T>&         B
) {
    for (std::size_t i = 0; i < N; ++i) {
        if (!(R[i] > T{0}) || !(C[i] > T{0})) {
            return false;
        }
    }
    A = Matrix<N, N, T>{};
    B = Matrix<N, 2, T>{};
    // Couplings between adjacent nodes via R[i] (i = 0..N-2).
    for (std::size_t i = 0; i + 1 < N; ++i) {
        const T g = T{1} / R[i];
        A(i, i) -= g / C[i];
        A(i, i + 1) += g / C[i];
        A(i + 1, i) += g / C[i + 1];
        A(i + 1, i + 1) -= g / C[i + 1];
    }
    // Last node returns to ambient through R[N-1]:
    //   C_{N-1} \dot T_{N-1} += -(T_{N-1} - T_amb)/R_{N-1}
    // so A picks up -1/(R C) and B(:,1) picks up +1/(R C).
    const T g_amb = (T{1} / R[N - 1]) / C[N - 1];
    A(N - 1, N - 1) -= g_amb;
    B(0, 0) = T{1} / C[0]; // power [W] at the junction node
    B(N - 1, 1) = g_amb;   // ambient temperature input
    return true;
}

} // namespace detail

/**
 * @brief Continuous Cauer RC ladder in absolute temperature with ambient input.
 *
 * Physical ladder: series thermal resistance @f$ R_i @f$ from node @f$ i @f$ to
 * node @f$ i+1 @f$ (or to ambient for @f$ i = N-1 @f$) and shunt capacitance
 * @f$ C_i @f$ at each node. States are absolute node temperatures; inputs are
 * @f$ u = [P,\, T_{\mathrm{amb}}]^\top @f$ (power at node 0, ambient at the
 * far end). The measured output is node @p YNode (0 = junction, @f$ N-1 @f$ =
 * case / heatsink end).
 *
 * @f[
 *   C_i\dot T_i = \frac{T_{i-1}-T_i}{R_{i-1}}
 *     - \frac{T_i - T_{i+1}}{R_i}\ (+\,P\ \text{at node 0}),
 * @f]
 * with @f$ T_N \equiv T_{\mathrm{amb}} @f$. Steady state:
 * @f$ T_0 = T_{\mathrm{amb}} + P\sum_i R_i @f$.
 *
 * Noise channels default to continuous @f$ G = I_N @f$, @f$ H = I @f$. After
 * discretize, pass the plant through @ref thermal_kalman_plant before
 * design::kalman so process noise is per-sample on the states (ZOH shrinks
 * continuous @f$ G = I @f$ by about @f$ T_s @f$).
 *
 * Continuous (@f$ T_s = 0 @f$); discretize with discretize (ZOH is exact
 * for piecewise-constant power and ambient).
 *
 * @note Compare with the rise-form cauer_thermal_ss (NU = 1, reference at 0).
 *
 * @tparam N     Number of thermal nodes / states.
 * @tparam YNode Measured node index in @f$ [0, N) @f$ (default 0 = junction).
 * @tparam T     Scalar type.
 * @param R [K/W] stage resistances (R[N-1] to ambient; each > 0).
 * @param C [J/K] node capacitances (each > 0).
 * @return Continuous StateSpace<N, 2, 1, double, N, 1>. Zeroed if any R or C is non-positive.
 */
template<std::size_t N, std::size_t YNode = 0, typename T = double>
    requires(YNode < N)
[[nodiscard]] constexpr StateSpace<N, 2, 1, T, N, 1> cauer_thermal_ss_ambient(
    const damp::array<T, N>& R,
    const damp::array<T, N>& C
) {
    StateSpace<N, 2, 1, T, N, 1> sys{};
    if (!detail::fill_cauer_ambient_ab(R, C, sys.A, sys.B)) {
        return StateSpace<N, 2, 1, T, N, 1>{};
    }
    sys.C(0, YNode) = T{1};
    // G = I_N, H = I_1 by StateSpace default when NW == NX and NV == NY.
    return sys;
}

/**
 * @brief Continuous Cauer RC ladder, absolute temps, all nodes as outputs.
 *
 * Same dynamics as cauer_thermal_ss_ambient with @f$ C = I_N @f$ so every
 * node temperature is an output (MIMO sensing — multiple NTCs). Inputs remain
 * @f$ u = [P,\, T_{\mathrm{amb}}]^\top @f$.
 *
 * @tparam N Number of thermal nodes / states / outputs.
 * @tparam T Scalar type.
 * @param R [K/W] stage resistances (each > 0).
 * @param C [J/K] node capacitances (each > 0).
 * @return Continuous StateSpace<N, 2, N, double, N, N>. Zeroed if any R or C is non-positive.
 */
template<std::size_t N, typename T = double>
[[nodiscard]] constexpr StateSpace<N, 2, N, T, N, N> cauer_thermal_ss_ambient_mimo(
    const damp::array<T, N>& R,
    const damp::array<T, N>& C
) {
    StateSpace<N, 2, N, T, N, N> sys{};
    if (!detail::fill_cauer_ambient_ab(R, C, sys.A, sys.B)) {
        return StateSpace<N, 2, N, T, N, N>{};
    }
    sys.C = Matrix<N, N, T>::identity();
    return sys;
}

/**
 * @brief Prepare a discretized thermal plant for design::kalman.
 *
 * Continuous ambient Cauer models ship with @f$ G = I @f$ (white noise on
 * @f$ \dot T @f$). Zero-order hold maps that to @f$ G_d \approx T_s\,I @f$, so a
 * naive discrete @f$ Q @f$ in @f$ \mathrm{K}^2 @f$ is attenuated by about
 * @f$ T_s^2 @f$ and the filter ignores measurements. This helper sets
 * @f$ G = I_{N_x} @f$ on the discrete plant so @p Q is a per-sample state
 * covariance in temperature squared.
 *
 * Requires @f$ N_W = N_X @f$ (the ambient Cauer default). @f$ H @f$ is left
 * unchanged (identity for the single/multi-NTC models).
 *
 * @param discrete_sys Discretized thermal StateSpace (@f$ T_s > 0 @f$).
 * @return Copy with @f$ G = I @f$ for discrete process noise.
 *
 * @see design::kalman
 * @see cauer_thermal_ss_ambient
 */
template<std::size_t NX, std::size_t NU, std::size_t NY, std::size_t NW, std::size_t NV, typename T>
    requires(NW == NX)
[[nodiscard]] constexpr StateSpace<NX, NU, NY, T, NW, NV> thermal_kalman_plant(
    const StateSpace<NX, NU, NY, T, NW, NV>& discrete_sys
) {
    return StateSpace<NX, NU, NY, T, NW, NV>{
        .A = discrete_sys.A,
        .B = discrete_sys.B,
        .C = discrete_sys.C,
        .D = discrete_sys.D,
        .G = Matrix<NX, NX, T>::identity(),
        .H = discrete_sys.H,
        .Ts = discrete_sys.Ts,
    };
}

} // namespace damp::design

namespace damp {

/**
 * @brief Runtime thermal-network Kalman observer.
 *
 * Thin wrapper around KalmanFilter for a discrete thermal StateSpace:
 * predict with known power and ambient, update from measured node temperature(s)
 * (case / junction NTC). States are absolute node temperatures.
 *
 * Typical deploy path:
 * @f[
 *   \text{cont} = \mathtt{cauer\_thermal\_ss\_ambient}(R,C)
 *   \xrightarrow{\mathrm{ZOH}} \text{disc}
 *   \xrightarrow{\mathtt{thermal\_kalman\_plant}} \text{plant}
 *   \xrightarrow{\mathtt{design::kalman}} \text{result}
 *   \xrightarrow{.as<float>()} \mathtt{ThermalKalmanObserver}
 * @f]
 *
 * @tparam NX Number of thermal states (nodes).
 * @tparam NU Number of inputs (2 for ambient Cauer: power + ambient).
 * @tparam NY Number of measured nodes.
 * @tparam NW Process-noise channels (default NX).
 * @tparam NV Measurement-noise channels (default NY).
 * @tparam T  Scalar type (default float for embedded).
 *
 * @see design::cauer_thermal_ss_ambient
 * @see design::kalman
 * @see KalmanFilter
 * @see "Optimal State Estimation" (Simon, 2006), §5
 */
template<
    std::size_t NX,
    std::size_t NU = 2,
    std::size_t NY = 1,
    typename T = float,
    std::size_t NW = NX,
    std::size_t NV = NY>
class ThermalKalmanObserver {
public:
    constexpr ThermalKalmanObserver() = default;

    /**
     * @brief Construct from a steady-state Kalman design result.
     * @param result Discrete thermal plant + L, P, Q, R from design::kalman.
     */
    constexpr explicit ThermalKalmanObserver(const design::KalmanResult<NX, NU, NY, T, NW, NV>& result)
        : filt_(result) {}

    /**
     * @brief Construct from a discrete thermal model and noise covariances.
     * @param discrete_sys Discrete thermal SS (Ts > 0), e.g. ZOH of ambient Cauer.
     * @param Q            Process noise covariance (NW × NW).
     * @param R            Measurement noise covariance (NV × NV).
     * @param x0           Initial absolute node temperatures.
     * @param P0           Initial error covariance.
     */
    constexpr ThermalKalmanObserver(
        const StateSpace<NX, NU, NY, T, NW, NV>& discrete_sys,
        const Matrix<NW, NW, T>&                 Q,
        const Matrix<NV, NV, T>&                 R,
        const ColVec<NX, T>&                     x0 = ColVec<NX, T>{},
        const Matrix<NX, NX, T>&                 P0 = Matrix<NX, NX, T>::identity()
    )
        : filt_(discrete_sys, Q, R, x0, P0) {}

    /**
     * @brief Time update with the full input vector.
     * @param u Control input (ambient Cauer: @f$ [P,\, T_{\mathrm{amb}}] @f$).
     */
    constexpr void predict(const ColVec<NU, T>& u) { filt_.predict(u); }

    /**
     * @brief Time update for the ambient-Cauer input shape @f$ u = [P,\, T_{\mathrm{amb}}] @f$.
     * @param power Power injected at the modelled heat-source node [W].
     * @param t_amb Ambient (or coolant / case-reference) temperature [°C or K].
     */
    constexpr void predict(T power, T t_amb)
        requires(NU == 2)
    {
        filt_.predict(ColVec<NU, T>{power, t_amb});
    }

    /**
     * @brief Measurement update with the full output vector.
     * @param y Measured node temperature(s).
     * @return false if the innovation covariance was singular (prediction kept).
     */
    [[nodiscard]] constexpr bool update(const ColVec<NY, T>& y) { return filt_.update(y); }

    /**
     * @brief Measurement update for a single NTC reading.
     * @param y_meas Measured node temperature [°C or K].
     * @return false if the innovation covariance was singular.
     */
    [[nodiscard]] constexpr bool update(T y_meas)
        requires(NY == 1)
    {
        return filt_.update(ColVec<NY, T>{y_meas});
    }

    /**
     * @brief Fused tick: predict with @p u then correct with @p y.
     * @return Updated absolute node-temperature estimate.
     */
    constexpr const ColVec<NX, T>& estimate(const ColVec<NY, T>& y, const ColVec<NU, T>& u) {
        return filt_.estimate(y, u);
    }

    /**
     * @brief Fused tick for ambient Cauer + single NTC.
     * @param y_meas Measured node temperature.
     * @param power  Power [W].
     * @param t_amb  Ambient temperature.
     * @return Updated absolute node-temperature estimate.
     */
    constexpr const ColVec<NX, T>& estimate(T y_meas, T power, T t_amb)
        requires(NU == 2 && NY == 1)
    {
        return filt_.estimate(ColVec<NY, T>{y_meas}, ColVec<NU, T>{power, t_amb});
    }

    /// Absolute node temperatures (state estimate).
    [[nodiscard]] constexpr const ColVec<NX, T>& temperatures() const { return filt_.state(); }

    /// Node-0 (junction) temperature.
    [[nodiscard]] constexpr T junction() const { return filt_.state()[0]; }

    /// Temperature of node @p i.
    [[nodiscard]] constexpr T node(std::size_t i) const { return filt_.state()[i]; }

    [[nodiscard]] constexpr const auto& filter() const { return filt_; }
    [[nodiscard]] constexpr auto&       filter() { return filt_; }

    constexpr void reset(
        const ColVec<NX, T>&     x0 = ColVec<NX, T>{},
        const Matrix<NX, NX, T>& P0 = Matrix<NX, NX, T>::identity()
    ) {
        filt_.reset(x0, P0);
    }

    /// Seed all nodes to the same absolute temperature (e.g. ambient at start).
    constexpr void reset_uniform(T t0, const Matrix<NX, NX, T>& P0 = Matrix<NX, NX, T>::identity()) {
        ColVec<NX, T> x0{};
        for (std::size_t i = 0; i < NX; ++i) {
            x0[i] = t0;
        }
        filt_.reset(x0, P0);
    }

private:
    KalmanFilter<NX, NU, NY, T, NW, NV> filt_{};
};

} // namespace damp
