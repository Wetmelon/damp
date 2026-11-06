#pragma once

/**
 * PlatformIO® / target compile-smoke profile (`targets/`).
 *
 * Product-ish embedded defaults:
 *   - DAMP_BACKEND_ETL          — damp::array/optional/... → etl::
 *   - DAMP_MATH_BACKEND_DAMP     — runtime float math via damp/math/trig.hpp
 *
 * Found via `-I.` in platformio.ini and `-Itargets` in `make target-smoke`.
 * Host SILs under examples/ keep their own damp_profile (stdlib + std:: math).
 *
 * Requires ETL headers on the include path (`-isystem../libs/etl/include` in PIO; any
 * copy of ETL works). DAMP_MATH_BACKEND_DAMP needs <cmath>, so this profile does not suit a
 * toolchain with no C++ standard headers: see "Using ETL" in README.md.
 */

#define DAMP_BACKEND_ETL
#define DAMP_MATH_BACKEND_DAMP
