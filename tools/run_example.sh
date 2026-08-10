#!/bin/sh
# Run one example binary, then write a tup stamp.
# Usage (CWD = examples/ as set by tup):
#   sh ../tools/run_example.sh <exe> <stamp>
#
# Portable: Linux / macOS / Windows+MSYS/Git-Bash. No Python, no .bat.
set -e

if [ "$#" -ne 2 ]; then
    echo "usage: run_example.sh <exe> <stamp>" >&2
    exit 2
fi

exe=$1
stamp=$2

if [ ! -f "$exe" ]; then
    echo "run_example: missing $exe" >&2
    exit 1
fi

base=$(basename "$exe")
echo "=== $base ==="

# Run via path as given (build/foo.exe). On Windows, .exe is fine under sh.
"./$exe" 2>/dev/null || "$exe"

mkdir -p "$(dirname "$stamp")"
echo ok >"$stamp"
