// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file iec61131.hpp
 * @brief IEC 61131-3 function-block style helpers
 */

#include <cstddef>
#include <cstdint>
#include <limits>

#include "damp/backend.hpp"

namespace damp::plc {

typedef bool     BOOL;  ///< Boolean type
typedef uint8_t  BYTE;  ///< 8-bit unsigned integer
typedef uint16_t WORD;  ///< 16-bit unsigned integer
typedef uint32_t DWORD; ///< 32-bit unsigned integer
typedef uint64_t LWORD; ///< 64-bit unsigned integer

typedef int8_t  SINT; ///< Signed 8-bit integer
typedef int16_t INT;  ///< Signed 16-bit integer
typedef int32_t DINT; ///< Signed 32-bit integer
typedef int64_t LINT; ///< Signed 64-bit integer

typedef uint8_t  USINT; ///< Unsigned 8-bit integer
typedef uint16_t UINT;  ///< Unsigned 16-bit integer
typedef uint32_t UDINT; ///< Unsigned 32-bit integer
typedef uint64_t ULINT; ///< Unsigned 64-bit integer
typedef float    REAL;  ///< Single-precision floating point
typedef double   LREAL; ///< Double-precision floating point

typedef char     CHAR;    ///< Character type
typedef wchar_t  WCHAR;   ///< Wide character type
typedef char*    STRING;  ///< String type (pointer to char)
typedef wchar_t* WSTRING; ///< Wide string type (pointer to wchar_t)

/**
 * @defgroup iec61131 IEC61131-3 Functions
 * @brief Allocation-free PLC-style function blocks (TON/TOF, CTU/CTD, …)
 *
 * Implementation of IEC 61131-3 standard function blocks commonly used in
 * industrial control and embedded systems. Lives under @c damp::plc.
 */

/**
 * @brief SR Latch (Set-dominant Set-Reset Latch)
 *
 * Bistable function block with set and reset inputs. Per IEC 61131-3, SR is
 * set-dominant: when S and R are both true, Set wins and Q1 stays true.
 *
 *     Q1 = S OR (Q1 AND NOT R)
 *
 * Contrast @ref RS, which is reset-dominant. (Earlier versions of this block
 * were incorrectly reset-dominant — identical to RS.)
 */
struct SR {
    bool Q1{false}; ///< Current state

    /**
     * @brief Execute SR latch
     * @param S Set input (dominant over R)
     * @param R Reset input
     * @return Current state
     */
    constexpr bool operator()(bool S, bool R) {
        if (S) {
            Q1 = true;
        } else if (R) {
            Q1 = false;
        }
        return Q1;
    }

    /// @brief Reset the latch to false.
    constexpr void reset() { Q1 = false; }
};

/**
 * @brief RS Latch (Reset-Set Latch)
 *
 * Bistable function block with reset and set inputs.
 * Q1 = NOT R AND (S OR Q1)
 */
struct RS {
    bool Q1{false}; ///< Current state

    /**
     * @brief Execute RS latch
     * @param R Reset input (dominant over S)
     * @param S Set input
     * @return Current state
     */
    constexpr bool operator()(bool R, bool S) {
        if (R) {
            Q1 = false;
        } else if (S) {
            Q1 = true;
        }
        return Q1;
    }

    /// @brief Reset the latch to false.
    constexpr void reset() { Q1 = false; }
};

/**
 * @brief R_TRIG (Rising Edge Trigger)
 *
 * Detects rising edges on the CLK input.
 */
struct R_TRIG {
    bool CLK{false}; ///< Clock input
    bool Q{false};   ///< Output (true on rising edge)

    /**
     * @brief Execute edge detector
     * @param CLK_input Clock input
     * @return Edge detected
     */
    constexpr bool operator()(bool CLK_input) {
        Q = CLK_input && !CLK;
        CLK = CLK_input;
        return Q;
    }

    /**
     * @brief Reset detector
     */
    constexpr void reset() {
        CLK = false;
        Q = false;
    }
};

/**
 * @brief F_TRIG (Falling Edge Trigger)
 *
 * Detects falling edges on the CLK input. Per IEC 61131-3 the body is
 * `Q := NOT CLK AND NOT M; M := NOT CLK`, with the internal memory `M`
 * initialized to 0 on a cold restart. Because `M = NOT CLK_prev`, an `M` of 0
 * corresponds to a *previous* clock of 1 — so the standard's NOTE specifies that
 * an F_TRIG whose CLK is FALSE produces `Q = 1` on its first execution after a
 * cold restart. We store the previous raw clock in `CLK`, so it is initialized
 * to `true` (not `false`) to reproduce that conforming cold-start behavior.
 */
struct F_TRIG {
    bool CLK{true}; ///< Previous clock state; init true so cold start matches IEC (M=0).
    bool Q{false};  ///< Output (true on falling edge)

    /**
     * @brief Execute edge detector
     * @param CLK_input Clock input
     * @return Edge detected
     */
    constexpr bool operator()(bool CLK_input) {
        Q = !CLK_input && CLK;
        CLK = CLK_input;
        return Q;
    }

    /**
     * @brief Reset detector (to the conforming cold-start state).
     */
    constexpr void reset() {
        CLK = true;
        Q = false;
    }
};

/**
 * @brief TON Timer (Timer On Delay)
 *
 * On-delay timer that activates output after specified time.
 */
template<typename T = float>
class TON {
public:
    T    PT{0};    ///< Preset time
    bool Q{false}; ///< Timer output
    T    ET{0};    ///< Elapsed time

    /**
     * @brief Execute timer
     * @param IN Timer input
     * @param dt Time step
     * @return Timer output
     */
    constexpr bool operator()(bool IN, T dt) {
        if (IN) {
            if (ET < PT) {
                ET += dt;
                if (ET >= PT) {
                    ET = PT; // clamp: ET saturates at PT, never grows unbounded
                }
            }
            Q = (ET >= PT);
        } else {
            ET = 0;
            Q = false;
        }
        return Q;
    }

    /**
     * @brief Reset timer
     */
    constexpr void reset() {
        Q = false;
        ET = 0;
    }
};

/**
 * @brief TOF Timer (Timer Off Delay)
 *
 * Off-delay timer that deactivates output after specified time.
 */
template<typename T = float>
class TOF {
public:
    T    PT{0};    ///< Preset time
    bool Q{false}; ///< Timer output
    T    ET{0};    ///< Elapsed time

    /**
     * @brief Execute timer
     * @param IN Timer input
     * @param dt Time step
     * @return Timer output
     */
    constexpr bool operator()(bool IN, T dt) {
        if (IN) {
            ET = 0;
            Q = true;
        } else {
            if (ET < PT) {
                ET += dt;
                if (ET >= PT) {
                    ET = PT; // clamp: ET saturates at PT, never grows unbounded
                }
            }
            if (ET >= PT) {
                Q = false;
            }
        }
        return Q;
    }

    /**
     * @brief Reset timer
     */
    constexpr void reset() {
        Q = false;
        ET = 0;
    }
};

/**
 * @brief TP Timer (Timer Pulse)
 *
 * Pulse timer that generates a pulse of specified duration.
 */
template<typename T = float>
class TP {
public:
    T    PT{0};    ///< Pulse time
    bool Q{false}; ///< Timer output
    T    ET{0};    ///< Elapsed time

    /**
     * @brief Execute timer
     * @param IN Timer input (rising edge triggers pulse)
     * @param dt Time step
     * @return Timer output
     */
    constexpr bool operator()(bool IN, T dt) {
        // Start a pulse only on a rising edge while idle. The pulse is
        // non-retriggerable: once running (Q true) IN is ignored until PT
        // elapses, and re-arming requires the input to return low (ET back to 0).
        if (IN && !prev_IN && !Q && ET == T{0}) {
            Q = true;
        }

        if (Q) {
            ET += dt;
            if (ET >= PT) {
                ET = PT; // clamp: ET saturates at PT, never grows unbounded
                Q = false;
            }
        } else if (!IN) {
            ET = T{0}; // input released after the pulse -> ET returns to 0 (idle)
        }

        prev_IN = IN;
        return Q;
    }

    /**
     * @brief Reset timer
     */
    constexpr void reset() {
        Q = false;
        ET = 0;
        prev_IN = false;
    }

private:
    bool prev_IN{false}; ///< Previous input state
};

/**
 * @brief CTU Counter (Count Up)
 *
 * Up-counter with reset and count limit.
 */
template<typename T = uint32_t>
class CTU {
public:
    T    PV{0};     ///< Preset value (count limit)
    T    CV{0};     ///< Current value
    bool Q{false};  ///< Counter output (CV >= PV)
    bool CU{false}; ///< Count up input (previous state)

    /**
     * @brief Execute counter
     * @param CU_input Count up input
     * @param R Reset input
     * @return Counter output
     */
    constexpr bool operator()(bool CU_input, bool R) {
        if (R) {
            CV = 0;
        } else if (CU_input && !CU && CV < std::numeric_limits<T>::max()) {
            ++CV; // rising edge on CU, saturating at the type max (no overflow)
        }
        CU = CU_input;
        Q = (CV >= PV); // per IEC, Q is evaluated every invocation
        return Q;
    }

    /**
     * @brief Reset counter
     */
    constexpr void reset() {
        CV = 0;
        Q = false;
        CU = false;
    }
};

/**
 * @brief CTD Counter (Count Down)
 *
 * Down-counter with load and count limit.
 */
template<typename T = uint32_t>
class CTD {
public:
    T    PV{0};     ///< Preset value (load value)
    T    CV{0};     ///< Current value
    bool Q{false};  ///< Counter output (CV <= 0)
    bool CD{false}; ///< Count down input (previous state)

    /**
     * @brief Execute counter
     * @param CD_input Count down input
     * @param LD Load input
     * @return Counter output
     */
    constexpr bool operator()(bool CD_input, bool LD) {
        if (LD) {
            CV = PV;
        } else if (CD_input && !CD && CV > 0) {
            --CV; // rising edge on CD, floored at 0 (PVmin); no unsigned wrap
        }
        CD = CD_input;
        Q = (CV == 0); // per IEC, Q := (CV <= 0) evaluated every invocation
        return Q;
    }

    /**
     * @brief Reset counter
     */
    constexpr void reset() {
        CV = 0;
        Q = false;
        CD = false;
    }
};

/**
 * @brief CTUD Counter (Count Up Down)
 *
 * Up-down counter with reset and load.
 */
template<typename T = uint32_t>
class CTUD {
public:
    T    PV{0};     ///< Preset value (load value)
    T    CV{0};     ///< Current value
    bool QU{false}; ///< Up counter output (CV >= PV)
    bool QD{true};  ///< Down counter output (CV <= 0) - initially true since CV=0
    bool CU{false}; ///< Count up input (previous state)
    bool CD{false}; ///< Count down input (previous state)

    /**
     * @brief Execute counter
     * @param CU_input Count up input
     * @param CD_input Count down input
     * @param R Reset input
     * @param LD Load input
     */
    constexpr void operator()(bool CU_input, bool CD_input, bool R, bool LD) {
        const bool cu_edge = CU_input && !CU; // rising edge on CU
        const bool cd_edge = CD_input && !CD; // rising edge on CD

        if (R) {
            CV = 0;
        } else if (LD) {
            CV = PV;
        } else if (!(cu_edge && cd_edge)) {
            // Per IEC, simultaneous up/down edges cancel (do nothing); otherwise
            // up takes priority (ELSIF), each saturating at its limit (no wrap).
            if (cu_edge && CV < std::numeric_limits<T>::max()) {
                ++CV;
            } else if (cd_edge && CV > 0) {
                --CV;
            }
        }

        CU = CU_input;
        CD = CD_input;
        QU = (CV >= PV); // per IEC, outputs are evaluated every invocation
        QD = (CV == 0);
    }

    /**
     * @brief Reset counter
     */
    constexpr void reset() {
        CV = 0;
        QU = false;
        QD = true; // CV == 0
        CU = CD = false;
    }
};

/**
 * @brief D Flip-Flop (edge-triggered data latch)
 *
 * Captures the data input D on each rising edge of CLK and holds it on Q until
 * the next rising edge.
 */
struct DFF {
    bool Q{false};   ///< Stored output
    bool CLK{false}; ///< Previous clock state

    /**
     * @brief Execute the flip-flop.
     * @param D   Data input (captured on rising CLK).
     * @param CLK_input Clock input.
     * @return Stored output Q.
     */
    constexpr bool operator()(bool D, bool CLK_input) {
        if (CLK_input && !CLK) { // rising edge
            Q = D;
        }
        CLK = CLK_input;
        return Q;
    }

    constexpr void reset() {
        Q = false;
        CLK = false;
    }
};

/**
 * @brief D Latch (level-sensitive / transparent latch)
 *
 * While the enable input E is true the latch is transparent (Q follows D); when
 * E is false it holds the last value. Unlike @ref DFF this is level- not
 * edge-triggered.
 */
struct DLATCH {
    bool Q{false}; ///< Stored / pass-through output

    /**
     * @brief Execute the latch.
     * @param D Data input.
     * @param E Enable (transparent while true, hold while false).
     * @return Output Q.
     */
    constexpr bool operator()(bool D, bool E) {
        if (E) {
            Q = D;
        }
        return Q;
    }

    constexpr void reset() { Q = false; }
};

/**
 * @brief T Flip-Flop (toggle on rising edge)
 *
 * Toggles Q on each rising edge of the T input while enabled. Acts as a
 * divide-by-two on a clock, or a press-to-toggle latch on a button edge.
 */
struct TFF {
    bool Q{false}; ///< Toggled output
    bool T{false}; ///< Previous trigger state

    /**
     * @brief Execute the toggle.
     * @param T_input Trigger input (Q toggles on its rising edge).
     * @return Output Q.
     */
    constexpr bool operator()(bool T_input) {
        if (T_input && !T) { // rising edge
            Q = !Q;
        }
        T = T_input;
        return Q;
    }

    constexpr void reset() {
        Q = false;
        T = false;
    }
};

/**
 * @brief BLINK (free-running square-wave / flasher)
 *
 * Generates a periodic boolean while enabled, with independent on/off times so
 * the duty cycle is arbitrary (status LEDs, audible-alarm cadence, test
 * stimulus). Holds output false and resets its phase when disabled.
 *
 * @tparam T Time scalar type (default: float).
 */
template<typename T = float>
class BLINK {
public:
    T    time_on{0};  ///< Duration of the high phase [s].
    T    time_off{0}; ///< Duration of the low phase [s].
    bool Q{false};    ///< Current output.

    constexpr BLINK() = default;

    /// @param t_on High duration [s]; @param t_off Low duration [s].
    constexpr BLINK(T t_on, T t_off) : time_on(t_on), time_off(t_off) {}

    /**
     * @brief Execute the flasher.
     * @param enable Run while true; hold low and reset phase while false.
     * @param dt     Time step [s].
     * @return Output Q.
     */
    constexpr bool operator()(bool enable, T dt) {
        if (!enable) {
            Q = false;
            elapsed_ = T{0};
            return Q;
        }
        elapsed_ += dt;
        const T phase_time = Q ? time_on : time_off;
        if (elapsed_ >= phase_time) {
            Q = !Q;
            elapsed_ = T{0};
        }
        return Q;
    }

    constexpr void reset() {
        Q = false;
        elapsed_ = T{0};
    }

private:
    T elapsed_{0};
};

/**
 * @brief LIMIT (IEC 61131-3 selection function): clamp @p in to [@p mn, @p mx].
 *
 * Equivalent to `MIN(MAX(in, mn), mx)`. Argument order matches the standard:
 * `LIMIT(MN, IN, MX)`.
 */
template<typename T>
[[nodiscard]] constexpr T LIMIT(T mn, T in, T mx) {
    return damp::min(damp::max(in, mn), mx);
}

/**
 * @brief SEL (IEC 61131-3 binary selection): @p g ? @p in1 : @p in0.
 */
template<typename T>
[[nodiscard]] constexpr T SEL(bool g, T in0, T in1) {
    return g ? in1 : in0;
}

/**
 * @brief MUX (IEC 61131-3 multiplexer): select input @p k of N (0-based).
 *
 * `MUX(K, IN0, IN1, …)` returns the Kth input. An out-of-range @p k selects
 * IN0 (clamped, never UB), matching the conservative PLC convention.
 */
template<typename T, typename... Ts>
[[nodiscard]] constexpr T MUX(size_t k, T in0, Ts... ins) {
    const T arr[] = {in0, static_cast<T>(ins)...};
    return arr[k < (sizeof...(ins) + 1) ? k : 0];
}

} // namespace damp::plc