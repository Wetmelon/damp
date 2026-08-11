// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file doxygen_groups.hpp
 * @brief Doxygen module group definitions (documentation only; no code).
 *
 * Declares module groups used by ingroup tags across the library so Doxygen
 * does not warn about missing groups. Groups already defined with titles in
 * domain headers (FOC, ACIM, modulation, and similar) are not repeated here.
 *
 * Expand briefs here when a group is shared; domain-specific groups keep their
 * @c @defgroup next to the owning header.
 */

/**
 * @defgroup linear_algebra Linear algebra
 * @brief Fixed-size stack matrices, vectors, views, decompositions, and solves.
 *
 * Core types are `Matrix`, `ColVec`, and non-owning views. Free algorithms live
 * in `damp::mat` (prefer solve over explicit inverse). Fully `constexpr`-capable.
 */

/**
 * @defgroup controllers Controllers (general)
 * @brief Runtime control blocks beyond simple discrete state-feedback wrappers.
 *
 * Includes composition helpers, governors, ESC, harmonic suppression, repetitive
 * control, and other tick-rate objects that are not pure LQR/PID shells.
 */

/**
 * @defgroup discrete_controllers Discrete controllers
 * @brief Sampled-data runtime controllers (PID, PR, LQR family, ADRC, SMC, Smith, …).
 *
 * Typically constructed from `design::` Results (or pack maps) via `.as<float>()`
 * and driven with `.control()` in the ISR or RTOS loop.
 */

/**
 * @defgroup estimators Estimators (general)
 * @brief Observers and estimators spanning continuous-time and hybrid use cases.
 *
 * DOB and related blocks that are not exclusively discrete Luenberger-style.
 */

/**
 * @defgroup discrete_estimators Discrete estimators
 * @brief Discrete-time observers such as full- and reduced-order Luenberger forms.
 */

/**
 * @defgroup trajectory Trajectory generators
 * @brief Motion profiles and path generators for multi-axis motion.
 *
 * Trapezoidal and S-curve profiles, polynomials, splines, electronic cam, TOPP,
 * input shapers, and online OTG helpers. Pair with toolbox actuator bridges for
 * drive-native units.
 */
