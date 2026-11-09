#!/bin/sh
# Host contract checks for the embeddable umbrella (damp/control.hpp).
#
# Usage (any CWD):
#   sh tools/contract_check.sh freestanding <g++> <stamp>
#   sh tools/contract_check.sh embedded     <g++> <stamp>
#   sh tools/contract_check.sh workbench    <g++> <stamp>
#   sh tools/contract_check.sh nostd        <g++> <stamp>
#
# <g++> is one argument and may carry flags: "clang++ --target=x86_64-w64-mingw32".
#
# freestanding: control.hpp must compile under ETL + freestanding math, with a constexpr
#               design and its tick instantiated.
# embedded:     no Damp header reachable from the default (stdlib) control.hpp includes a heap
#               container header.
# nostd:        checks/nostd/nostd_check.cpp must build and run with no C++ standard
#               headers at all (-nostdinc++, ETL_NO_STL): a constexpr design and its tick. An
#               unguarded hosted #include anywhere in the core fails it.
set -e

mode=$1
GXX=$2
OUT=$3

if [ -z "$mode" ] || [ -z "$GXX" ] || [ -z "$OUT" ]; then
    echo "usage: contract_check.sh <freestanding|embedded|workbench|nostd> <g++> <stamp>" >&2
    exit 2
fi

# Resolve repo root from this script's location (works from checks/ or repo root).
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
cd "$ROOT"

# Stamp path may be relative to caller's CWD — normalize after cd if needed.
case "$OUT" in
    /* | [A-Za-z]:*) ;;
    *)
        # Relative stamp: re-root from original caller's intent via ROOT/checks/...
        # Tup runs with CWD=checks/, so OUT is build/freestanding.ok → checks/build/...
        if [ -d checks ] && echo "$OUT" | grep -q '^build/'; then
            OUT="checks/$OUT"
        fi
        ;;
esac

mkdir -p "$(dirname "$OUT")"
probe="$(dirname "$OUT")/probe_$mode.cpp" # per mode: tup runs these checks in parallel

case "$mode" in
freestanding)
    printf '%s\n' \
        '#define DAMP_BACKEND_ETL' \
        '#define DAMP_MATH_BACKEND_FREESTANDING' \
        '#include "damp/control.hpp"' \
        'using namespace damp;' \
        'constexpr auto sys = StateSpace{.A = Matrix<2, 2>{{0.0, 1.0}, {9.81, -0.1}}, .B = Matrix<2, 1>{{0.0}, {1.0}}, .C = Matrix<1, 2>{{1.0, 0.0}}};' \
        'constexpr auto res = design::lqr(sys, 0.01, {.Q = Matrix<2, 2>::identity() * 10.0, .R = 1.0});' \
        'static_assert(res.success);' \
        'int main() { LQR law{res.as<float>()}; law.update(ColVec<2, float>{0.1F, 0.0F}); return damp::isfinite(law.control()(0, 0)) ? 0 : 1; }' >"$probe"

    if ! $GXX -std=c++20 -Iinc -Ilibs/etl/include -fsyntax-only "$probe" 2>"${probe}.log"; then
        cat "${probe}.log" >&2
        echo "freestanding-check FAILED: the umbrella, a design or its tick do not compile under ETL + freestanding math" >&2
        rm -f "$probe" "${probe}.log"
        exit 1
    fi

    rm -f "$probe" "${probe}.log"
    echo "freestanding-check OK: a design and its tick compile on ETL + series math"
    echo ok >"$OUT"
    ;;

embedded)
    printf '%s\n' '#include "damp/control.hpp"' 'int main(){}' >"$probe"

    if ! deps=$($GXX -std=c++20 -Iinc -M "$probe" 2>"${probe}.log"); then
        cat "${probe}.log" >&2
        rm -f "$probe" "${probe}.log"
        echo "embedded-check FAILED: damp/control.hpp does not compile" >&2
        exit 1
    fi

    # Damp's own headers only: what a standard header reaches varies by library (GCC 11's
    # <algorithm> reaches stl_vector.h), and libc++ has no stl_vector.h to look for.
    heap='^[[:space:]]*#[[:space:]]*include[[:space:]]*<(vector|string|deque|list|map|set|unordered_map|unordered_set|memory|functional)>'
    if printf '%s\n' "$deps" | tr '\\ ' '/\n' | grep -E 'damp/.*\.hpp$' | xargs grep -lE "$heap" >&2; then
        rm -f "$probe" "${probe}.log"
        echo "embedded-check FAILED: damp/control.hpp reaches the heap container headers in the files above" >&2
        exit 1
    fi

    rm -f "$probe" "${probe}.log"
    echo "embedded-check OK: damp/control.hpp is allocation-free"
    echo ok >"$OUT"
    ;;

workbench)
    printf '%s\n' '#include "damp/workbench.hpp"' 'int main(){}' >"$probe"

    if ! $GXX -std=c++20 -Iinc -fsyntax-only "$probe" 2>"${probe}.log"; then
        cat "${probe}.log" >&2
        rm -f "$probe" "${probe}.log"
        echo "workbench-check FAILED: damp/workbench.hpp does not compile on its own" >&2
        exit 1
    fi

    rm -f "$probe" "${probe}.log"
    echo "workbench-check OK: damp/workbench.hpp is a self-contained host umbrella"
    echo ok >"$OUT"
    ;;

nostd)
    exe="${probe%.cpp}.exe"
    if ! $GXX -std=c++20 -nostdinc++ -fno-exceptions -fno-rtti -Ichecks/nostd -Iinc -Ilibs/etl/include \
        checks/nostd/nostd_check.cpp -o "$exe" 2>"${probe}.log"; then
        cat "${probe}.log" >&2
        rm -f "$exe" "${probe}.log"
        echo "nostd-check FAILED: the core, a design or its tick need a C++ standard header" >&2
        exit 1
    fi

    if ! "$exe"; then
        rm -f "$exe" "${probe}.log"
        echo "nostd-check FAILED: the tick gave a wrong answer with no C++ standard headers" >&2
        exit 1
    fi

    rm -f "$exe" "${probe}.log"
    echo "nostd-check OK: a constexpr design and its tick run with no C++ standard headers"
    echo ok >"$OUT"
    ;;

*)
    echo "unknown mode: $mode" >&2
    exit 2
    ;;
esac
