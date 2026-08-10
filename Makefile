.PHONY: all format build clean tests tidy iwyu fix_includes docs gui embedded-check freestanding-check doxygen-style-check reference-check ci targets target-smoke examples sils examples-force

# Compiler from tup.config (same keys as Tuprules.lua). Missing file is OK.
#   CONFIG_CXX                 full override (clang++, path/to/clang++)
#   CONFIG_COMPILER_PATH       directory of the toolchain
#   CONFIG_COMPILER_PREFIX     e.g. arm-none-eabi-
# Default: g++ on PATH.
#
# tup.config is gitignored. After `git clean -x` it is gone and every Tupfile
# no-ops unless CONFIG_DAMP_SELF_BUILD=y (submodule guard — docs/submodule_tup.md).
# Bootstrap from the tracked default so `make` works on a clean tree.
ifeq ($(wildcard tup.config),)
ifeq ($(wildcard tup.config.default),)
$(error tup.config missing and tup.config.default not found — cannot enable self-build)
endif
$(info [make] tup.config missing; copying tup.config.default (CONFIG_DAMP_SELF_BUILD=y))
$(shell cp tup.config.default tup.config)
endif

-include tup.config
CONFIG_COMPILER_PATH ?=
CONFIG_COMPILER_PREFIX ?=
CONFIG_CXX ?=
ifneq ($(CONFIG_CXX),)
GXX := $(CONFIG_CXX)
else ifneq ($(CONFIG_COMPILER_PATH),)
GXX := $(CONFIG_COMPILER_PATH)/$(CONFIG_COMPILER_PREFIX)g++
else
GXX := $(CONFIG_COMPILER_PREFIX)g++
endif

# Host tools (refs, doxygen-style-check, gui bootstrap). Override: make refs PYTHON=...
ifeq ($(OS),Windows_NT)
PYTHON ?= py -3
VENV_PY := .venv/Scripts/python.exe
else
PYTHON ?= python3
VENV_PY := .venv/bin/python
endif

# `build` runs full tup: tests, examples (run stamps), and contract checks
# (checks/build/{embedded,freestanding}.ok — only re-run when headers change).
all: format build refs
	@./tests/build/test_runner.exe
	@./tests/damp_backend/build/test_damp_backend.exe
	@./tests/etl_backend/build/test_etl_backend.exe

format:
	@clang-format -i $$(find inc $(wildcard sil) -name '*.hpp') \
		$$(find tests examples -name '*.cpp' -o -name '*.hpp')

compiledb:
	@tup --quiet --no-environ-check compiledb

build: format compiledb
	@tup --quiet --no-environ-check

# Build and run examples. Plot HTML under examples/plots/<domain>/ is a *real*
# tup output of each run rule (see examples/Tupfile.lua plot_outputs): up-to-date
# → skip; delete a plot (or its stamp) → that example re-runs on the next make.
examples: format compiledb
	@tup --quiet --no-environ-check examples

# SIL plots that exist in this tree (tup targets HTML — re-runs if missing/stale).
# Public v0.1 is sparse (no PE pack demos); monorepo may add power/* targets.
sils:
	@tup --quiet --no-environ-check \
		examples/plots/control/cart_pole_lqr.html \
		examples/plots/control/pendulum_sim.html \
		examples/plots/estimation/imu_pose_3d.html

# Wipe every example run product (stamps + plots), then rebuild through tup.
# Use when you want a full refresh without editing sources.
examples-force:
	@rm -f examples/build/run/*.ok
	@rm -rf examples/plots/power examples/plots/pmac examples/plots/motor \
		examples/plots/motion examples/plots/control examples/plots/estimation
	@tup --quiet --no-environ-check examples

tests: format compiledb
	@tup --quiet --no-environ-check tests
	@./tests/build/test_runner.exe
	@./tests/damp_backend/build/test_damp_backend.exe
	@./tests/etl_backend/build/test_etl_backend.exe

docs: refs
	@mkdir -p docs/html
	@doxygen Doxyfile

refs:
	@$(PYTHON) tools/gen_reference.py

# Contract stamps live in the tup DAG (checks/Tupfile.lua). These targets force
# that branch only — useful for a quick local check without a full rebuild.
embedded-check:
	@tup --quiet --no-environ-check checks/build/embedded.ok

freestanding-check:
	@tup --quiet --no-environ-check checks/build/freestanding.ok


# Doxygen house style under inc/damp (member tags, no **bold**, MATLAB notes, @file).
doxygen-style-check:
	@$(PYTHON) tools/doxygen_style_check.py

# Regenerate REFERENCE*.md and fail if the committed copies are stale.
reference-check: refs
	@git diff --exit-code -- REFERENCE.md REFERENCE_INDEX.md

# Lean host CI gate (matches .github/workflows/ci.yml). Does not format sources
# (local `make` still does) and does not build/run full examples/ or PIO boards.
# Covers: style, REFERENCE drift, tests + backends, embed/freestanding stamps,
# host target-smoke (product sketches). Requires tup init on a fresh clone.
ci: doxygen-style-check reference-check
	@tup --quiet --no-environ-check tests
	@./tests/build/test_runner.exe
	@./tests/damp_backend/build/test_damp_backend.exe
	@./tests/etl_backend/build/test_etl_backend.exe
	@tup --quiet --no-environ-check checks/build/embedded.ok checks/build/freestanding.ok
	@$(MAKE) --no-print-directory target-smoke




# Run clang-tidy with --fix over all .cpp files (and inc/ headers they pull in).
# Lives here, not in tup: clang-tidy --fix rewrites sources in place, which tup's
# input/output tracking forbids (same reason clang-format -i is a make step).
# clang on Windows targets MSVC by default, so we hand it the build compiler's
# target triple and system include paths or it won't find <cmath> et al.
tidy:
	@tup --quiet --no-environ-check compiledb
	@TGT=$$($(GXX) -dumpmachine); \
	ISYS=$$(echo | $(GXX) -std=c++20 -E -x c++ -v - 2>&1 \
		| sed -n '/search starts here:/,/End of search list/p' \
		| grep '^ ' | sed 's,^ ,-extra-arg=-isystem,'); \
	run-clang-tidy -p . -fix -header-filter='inc[/\\].*' \
		-extra-arg=--target=$$TGT -extra-arg=-fbracket-depth=16384 $$ISYS '\.cpp$$'

# Run include-what-you-use over every TU in compile_commands.json. The headers
# are already annotated with IWYU pragmas (keep/export); this drives the tool.
#
# IWYU is clang-based but the project builds with mingw-GCC, so we point clang's
# driver at the GCC toolchain (--gcc-toolchain) to analyze against the SAME
# libstdc++ the real build uses -- otherwise on Windows clang defaults to the
# MSVC STL and the include suggestions are wrong. The conda-forge IWYU package
# ships no clang resource headers (mm_malloc.h / intrinsics) and bakes a dead
# -resource-dir, so we borrow the system LLVM's resource dir instead.
#
# Needs: IWYU on PATH (conda-forge include-what-you-use) and a system clang
# (the one already used for tidy/format) for -print-resource-dir.
# Scope to tests/*.cpp: those TUs transitively pull all of inc/, so they cover
# the library, while keeping examples/'s third-party TUs (fmt/json/plotlypp) out
# of the sweep. IWYU reports each TU's filename relative to its compile dir
# (tests/), which is why fix_includes runs with `-p tests` below.
GCC_ROOT := $(dir $(CONFIG_COMPILER_PATH))
IWYU_TUS := $(wildcard tests/*.cpp)
iwyu:
	@tup --quiet --no-environ-check compiledb
	@export PYTHONUTF8=1; TGT=$$($(GXX) -dumpmachine); \
	RESDIR=$$(clang -print-resource-dir); \
	iwyu_tool.py -j 0 -p . $(IWYU_TUS) -- --target=$$TGT --gcc-toolchain="$(GCC_ROOT)" \
		-resource-dir="$$RESDIR" -fbracket-depth=16384 \
		-Xiwyu --no_fwd_decls -Xiwyu --mapping_file=$(CURDIR)/iwyu.imp

# Same IWYU sweep, but apply the suggestions: pipe IWYU's output straight into
# fix_includes.py (in place, like clang-tidy --fix -- a make step, not tup).
# --noreorder leaves include ordering to clang-format; `format` runs after to
# sort + align the edits.
fix_includes:
	@tup --quiet --no-environ-check compiledb
	@export PYTHONUTF8=1; TGT=$$($(GXX) -dumpmachine); \
	RESDIR=$$(clang -print-resource-dir); \
	iwyu_tool.py -j 0 -p . $(IWYU_TUS) -- --target=$$TGT --gcc-toolchain="$(GCC_ROOT)" \
		-resource-dir="$$RESDIR" -fbracket-depth=16384 \
		-Xiwyu --no_fwd_decls -Xiwyu --mapping_file=$(CURDIR)/iwyu.imp \
		| fix_includes.py --noreorder -p tests
	@$(MAKE) --no-print-directory format

# Host smoke: product *_sketch.cpp files with targets/damp_profile (ETL + Damp math).
# Arduino APIs stubbed when ARDUINO is undefined. Does not need PlatformIO.
# -Itargets before other profiles so config.hpp picks targets/damp_profile.hpp.
TARGET_SMOKE_INCLUDES = -Iinc -Itargets -Ilibs/etl/include
target-smoke:
	@mkdir -p targets/build
	@$(GXX) -std=c++20 -O3 -Wall -Wextra $(TARGET_SMOKE_INCLUDES) -Iexamples/control/pid \
		-o targets/build/pid_sketch_host.exe examples/control/pid/pid_sketch.cpp
	@$(GXX) -std=c++20 -O3 -Wall -Wextra $(TARGET_SMOKE_INCLUDES) -Iexamples/control/cart_pole \
		-o targets/build/cart_pole_sketch_host.exe examples/control/cart_pole/cart_pole_sketch.cpp
	@./targets/build/pid_sketch_host.exe
	@./targets/build/cart_pole_sketch_host.exe
	@echo "target-smoke: host compile+run of pid + cart_pole product sketches (ETL profile) OK"

# PlatformIO board compile-smoke (optional; needs `pio` on PATH).
# Default = full green matrix in targets/platformio.ini (Teensy 4.x, STM32, Pico, ESP32-S3).
# Subset: make targets ENVS="teensy41 nucleo_g474re"
# -s silences per-file "Compiling …" noise; failures still print.
ENVS ?=
targets:
ifeq ($(ENVS),)
	@cd targets && pio run -s
else
	@cd targets && pio run -s $(foreach e,$(ENVS),-e $(e))
endif

# Servo GUI companion lived under examples/servo_drive (not in public core tree).
gui:
	@echo "make gui: examples/servo_drive is not in this public tree."
	@exit 1