#!/bin/sh
# Host contract checks for the embeddable umbrella (damp/control.hpp).
#
# Usage (any CWD):
#   sh tools/contract_check.sh freestanding <g++> <stamp>
#   sh tools/contract_check.sh embedded     <g++> <stamp>
#
# freestanding: ETL + freestanding math must compile; no damp/ header reachable
#               from control.hpp may unconditionally #include <hosted-std>.
# embedded:     default (stdlib) control.hpp must not pull <vector>.
set -e

mode=$1
GXX=$2
OUT=$3

if [ -z "$mode" ] || [ -z "$GXX" ] || [ -z "$OUT" ]; then
    echo "usage: contract_check.sh <freestanding|embedded> <g++> <stamp>" >&2
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
probe="$(dirname "$OUT")/probe.cpp"

case "$mode" in
freestanding)
    printf '%s\n' \
        '#define DAMP_BACKEND_ETL' \
        '#define DAMP_MATH_BACKEND_FREESTANDING' \
        '#include "damp/control.hpp"' \
        'int main(){}' >"$probe"

    if ! "$GXX" -std=c++20 -Iinc -Ilibs/etl/include -fsyntax-only "$probe" 2>"${probe}.log"; then
        cat "${probe}.log" >&2
        echo "freestanding-check FAILED: umbrella does not compile under ETL + freestanding math" >&2
        rm -f "$probe" "${probe}.log"
        exit 1
    fi

    deps=$("$GXX" -std=c++20 -Iinc -Ilibs/etl/include -M "$probe" 2>/dev/null \
        | tr ' \\' '\n' \
        | grep -E 'inc[/\\]damp[/\\].*\.hpp' \
        | sort -u) || true

    rm -f "${probe}.bad"
    # shellcheck disable=SC2086
    for f in $deps; do
        case "$f" in
            *[/\\]backend.hpp) continue ;;
        esac
        [ -f "$f" ] || continue

        hits=$(grep -nE '#include[[:space:]]*<(cmath|algorithm|vector|string|numbers|optional|tuple|memory|functional|complex|valarray|array|iostream|sstream|map|set|deque|list)>' "$f" 2>/dev/null || true)
        [ -z "$hits" ] && continue

        echo "$hits" | while IFS= read -r line; do
            h=$(printf '%s\n' "$line" | sed -n 's/.*#include[[:space:]]*<\([^>]*\)>.*/\1/p')
            [ -n "$h" ] || continue
            if ! grep -q "__has_include(<${h}>)" "$f" 2>/dev/null; then
                echo "  LEAK: $f unconditionally includes <$h>" >&2
                echo 1 >"${probe}.bad"
            fi
        done
    done

    if [ -f "${probe}.bad" ]; then
        rm -f "$probe" "${probe}.log" "${probe}.bad"
        echo "freestanding-check FAILED: hosted headers reachable from damp/control.hpp" >&2
        exit 1
    fi

    rm -f "$probe" "${probe}.log" "${probe}.bad"
    echo "freestanding-check OK: damp/control.hpp is freestanding-clean (ETL + series math)"
    echo ok >"$OUT"
    ;;

embedded)
    printf '%s\n' '#include "damp/control.hpp"' 'int main(){}' >"$probe"

    if "$GXX" -std=c++20 -Iinc -M "$probe" 2>/dev/null | grep -q stl_vector.h; then
        rm -f "$probe"
        echo "embedded-check FAILED: <vector> is reachable from damp/control.hpp" >&2
        exit 1
    fi

    rm -f "$probe"
    echo "embedded-check OK: damp/control.hpp is allocation-free"
    echo ok >"$OUT"
    ;;

*)
    echo "unknown mode: $mode" >&2
    exit 2
    ;;
esac
