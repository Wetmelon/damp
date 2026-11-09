#!/bin/sh
# Cortex-M7 disassembly for checks/cycle_guesser.py. It is for inspection, not a contract, so
# without arm-none-eabi-g++ on PATH it writes a "skipped" file and passes.
#
# Usage (tup, CWD = checks/):
#   sh ../tools/arm_asm_dump.sh <source.cpp> <out.asm> [compile flags...]
set -e

SRC=$1
OUT=$2
shift 2

if ! command -v arm-none-eabi-g++ >/dev/null 2>&1; then
    echo "skipped: no arm-none-eabi-g++ on PATH" >"$OUT"
    exit 0
fi

obj="${OUT%.asm}.o"
arm-none-eabi-g++ -c "$SRC" "$@" -o "$obj"
arm-none-eabi-objdump -drC "$obj" >"$OUT"
rm -f "$obj"
