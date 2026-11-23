// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Damp after Arduino.h, which an .ino sketch always gets first. The Arduino cores define these
// as macros (Renesas, SAMD and AVR alike; B0 … B11111111 are SAMD's and AVR's binary.h), and a
// Damp name one of them rewrites does not compile: this TU is the check.

#include <math.h>
#include <stdlib.h>

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

// clang-format off
#define abs(x) ((x)>0?(x):-(x))
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#define constrain(amt,low,high) ((amt)<(low)?(low):((amt)>(high)?(high):(amt)))
#define radians(deg) ((deg)*DEG_TO_RAD)
#define degrees(rad) ((rad)*RAD_TO_DEG)
#define sq(x) ((x)*(x))
#define bit(b) (1UL << (b))
#define lowByte(w) ((uint8_t) ((w) & 0xff))
#define highByte(w) ((uint8_t) ((w) >> 8))
#define bitRead(value, bit) (((value) >> (bit)) & 0x01)
#define F(string_literal) (string_literal)
#define PI 3.1415926535897932384626433832795
#define HALF_PI 1.5707963267948966192313216916398
#define TWO_PI 6.283185307179586476925286766559
#define DEG_TO_RAD 0.017453292519943295769236907684886
#define RAD_TO_DEG 57.295779513082320876798154814105
#define EULER 2.718281828459045235360287471352
#define HIGH 0x1
#define LOW 0x0
#define INPUT 0x0
#define OUTPUT 0x1
#define CHANGE 2
#define FALLING 3
#define RISING 4
#define DEFAULT 1
#define EXTERNAL 0
#define LSBFIRST 0
#define MSBFIRST 1
#define B0 0
#define B1 1
#define B10 2
#define B11 3
#define B100 4
#define B101 5
#define B110 6
#define B111 7
// clang-format on

#include "damp/control.hpp"

using namespace damp;

namespace {

constexpr auto sys = StateSpace{
    .A = Matrix<2, 2>{{0.0, 1.0}, {9.81, -0.1}},
    .B = Matrix<2, 1>{{0.0}, {1.0}},
    .C = Matrix<1, 2>{{1.0, 0.0}},
};
constexpr auto res = design::state_feedback(sys, 0.01, {.Q = Matrix<2, 2>::identity() * 10.0, .R = 1.0});

} // namespace

TEST_CASE("Damp after Arduino's macros: an LQR and a PI design and tick") {
    Law<res, float> law{};
    const float     u = law.step(ColVec<2, float>{0.1f, 0.0f})(0);
    CHECK(u < 0.0f);
    CHECK(damp::abs(u) > 0.0f);

    PIController<float> pi{design::pid(2.0f, 5.0f, 0.0f).discretize(0.01f)};
    CHECK(pi.control(1.0f, 0.0f) > 0.0f);
}
