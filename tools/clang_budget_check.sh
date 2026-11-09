#!/bin/sh
# Constexpr budget tripwire: a product sketch must compile with clang at clang's default
# -fconstexpr-steps (1048576). Its designs are constant-evaluated, and clang counts steps
# where GCC counts operations, so a design can fit GCC and still fail clang.
#
# Usage (tup, CWD = examples/):
#   sh ../tools/clang_budget_check.sh <c++ driver> <stamp> <sketch.cpp> [compile flags...]
#
# The driver names the target (-dumpmachine) and, on MinGW, the toolchain root, so clang
# analyzes against the same libstdc++ whatever is on PATH. Without clang++ 19 or newer on
# PATH (the oldest AGENTS.md supports) the check writes a "skipped" stamp and passes.
set -e

GXX=$1
OUT=$2
SRC=$3
shift 3

mkdir -p "$(dirname "$OUT")"

if ! command -v clang++ >/dev/null 2>&1; then
    echo "skipped: no clang++ on PATH" >"$OUT"
    exit 0
fi
major=$(clang++ -dumpversion 2>/dev/null | cut -d. -f1)
case "$major" in
    '' | *[!0-9]*)
        echo "skipped: clang++ version unknown" >"$OUT"
        exit 0
        ;;
esac
if [ "$major" -lt 19 ]; then
    echo "skipped: clang++ $major is older than 19" >"$OUT"
    exit 0
fi

target=$("$GXX" -dumpmachine)
# clang's MinGW driver looks for a triple-prefixed gcc on PATH, not the build driver. Linux
# keeps clang's own detection: this derivation would give it --sysroot=/usr.
case "$target" in
    *mingw*) set -- --sysroot="$("$GXX" -print-libgcc-file-name | sed 's,/lib/gcc/.*,,')" "$@" ;;
esac
if ! clang++ --target="$target" -std=c++20 -fsyntax-only "$@" "$SRC"; then
    echo "clang-budget-check FAILED: $SRC does not compile with clang $major at the default -fconstexpr-steps" >&2
    exit 1
fi
echo "ok: clang $major" >"$OUT"
