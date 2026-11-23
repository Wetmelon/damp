#!/bin/sh
# clang++ for the clang build of the test suites (clang-tests/Tupfile.lua, make clang-tests), aimed
# at the build driver's target and, on MinGW, its toolchain root, as tools/clang_budget_check.sh is.
#
# Usage (tup, CWD = clang-tests/):
#   sh ../tools/clang_cxx.sh <c++ driver> <output> <clang++ arguments...>
#
# Without clang++ 19 or newer on PATH (the oldest AGENTS.md supports) it writes an empty <output>
# and passes, so a build without clang has nothing to run: make clang-tests reports the suites
# skipped. A link of an empty object (compiled before clang was installed) is skipped the same way.
set -e

GXX=$1
OUT=$2
shift 2

skip() {
    : >"$OUT"
    exit 0
}

command -v clang++ >/dev/null 2>&1 || skip
major=$(clang++ -dumpversion 2>/dev/null | cut -d. -f1)
case "$major" in
    '' | *[!0-9]*) skip ;;
esac
[ "$major" -ge 19 ] || skip
link=yes
for arg in "$@"; do
    if [ "$arg" = -c ]; then
        link=
    fi
done
# An empty object: skipped. Its size through a read, not test -s, whose stat tup sees as a write.
if [ -n "$link" ]; then
    for arg in "$@"; do
        case "$arg" in
            *.o)
                if [ "$(wc -c <"$arg")" -eq 0 ]; then
                    skip
                fi
                ;;
        esac
    done
fi

target=$("$GXX" -dumpmachine)
# clang's MinGW driver looks for a triple-prefixed gcc on PATH, not the build driver. Emulated TLS:
# Clang 19's native TLS crashes doctest's thread_local on that target (Tuprules.lua). Linux keeps
# clang's own detection: this derivation would give it --sysroot=/usr.
case "$target" in
    *mingw*) set -- --sysroot="$("$GXX" -print-libgcc-file-name | sed 's,/lib/gcc/.*,,')" -femulated-tls "$@" ;;
esac
exec clang++ --target="$target" "$@"
