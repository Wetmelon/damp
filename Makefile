# Damp build driver. tup does the building; make sequences it with the steps tup
# can't own (formatting and other in-place source rewrites, docs, board builds).
#
#   make                 format, build everything, check docs, run all tests  ← final check
#   make tests           build and run the test suites only
#   make examples        build and run the examples
#   make ci              the lean gate: style, REFERENCE drift, tests, contracts, smoke
#   make tidy | iwyu | fix_includes    clang-based hygiene (see "Code hygiene")
#
# On Windows run make from PowerShell or cmd, not Git Bash (see AGENTS.md).

.PHONY: all format compiledb build tests examples examples-force sils \
        refs docs doxygen-style-check reference-check embedded-check freestanding-check ci \
        tidy iwyu fix_includes target-smoke targets

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# tup.config is gitignored. Without it every Tupfile no-ops (submodule guard,
# docs/submodule_tup.md), so bootstrap it from the tracked default.
ifeq ($(wildcard tup.config),)
ifeq ($(wildcard tup.config.default),)
$(error tup.config missing and tup.config.default not found — cannot enable self-build)
endif
$(info [make] tup.config missing; copying tup.config.default)
$(shell cp tup.config.default tup.config)
endif

# Compiler, from the same tup.config keys Tuprules.lua reads:
#   CONFIG_CXX               full override (clang++, path/to/clang++)
#   CONFIG_COMPILER_PATH     toolchain bin directory
#   CONFIG_COMPILER_PREFIX   cross prefix, e.g. arm-none-eabi-
-include tup.config
ifneq ($(CONFIG_CXX),)
GXX := $(CONFIG_CXX)
else ifneq ($(CONFIG_COMPILER_PATH),)
GXX := $(CONFIG_COMPILER_PATH)/$(CONFIG_COMPILER_PREFIX)g++
else
GXX := $(CONFIG_COMPILER_PREFIX)g++
endif

# Python for the tools/ scripts. Override: make refs PYTHON=python3
ifeq ($(OS),Windows_NT)
PYTHON ?= py -3
else
PYTHON ?= python3
endif

TUP := tup --quiet --no-environ-check

# The three doctest runners: main suite, damp math backend, ETL backend.
define run_tests
	@./tests/build/test_runner.exe
	@./tests/damp_backend/build/test_damp_backend.exe
	@./tests/etl_backend/build/test_etl_backend.exe
endef

# ---------------------------------------------------------------------------
# Everyday
# ---------------------------------------------------------------------------

# Full build = tests, examples (built and run) and the embedded/freestanding
# contract checks, all in the tup DAG.
all: format build refs doxygen-style-check
	$(run_tests)

build: format compiledb
	@$(TUP)

tests: format compiledb
	@$(TUP) tests
	$(run_tests)

format:
	@clang-format -i $$(find inc $(wildcard sil) -name '*.hpp') \
		$$(find tests examples -name '*.cpp' -o -name '*.hpp')

compiledb:
	@$(TUP) compiledb

# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------

# Each example run writes its plots under examples/plots/<domain>/ as real tup
# outputs: up to date → skipped; delete a plot or its stamp → it re-runs.
examples: format compiledb
	@$(TUP) examples

# Re-run every example: drop the run stamps and plots, then rebuild.
examples-force:
	@rm -f examples/build/run/*.ok
	@rm -rf examples/plots/control examples/plots/estimation examples/plots/motor
	@$(TUP) examples

# Just the headline SIL plots.
sils:
	@$(TUP) examples/plots/control/cart_pole_lqr.html \
		examples/plots/control/pendulum_sim.html \
		examples/plots/estimation/imu_pose_3d.html

# ---------------------------------------------------------------------------
# Docs and checks
# ---------------------------------------------------------------------------

refs:
	@$(PYTHON) tools/gen_reference.py

docs: refs
	@mkdir -p docs/html
	@doxygen Doxyfile

# Doxygen house style under inc/damp (member tags, no **bold**, MATLAB notes, @file).
doxygen-style-check:
	@$(PYTHON) tools/doxygen_style_check.py

# Regenerate REFERENCE*.md and fail if the committed copies are stale.
reference-check: refs
	@git diff --exit-code -- REFERENCE.md REFERENCE_INDEX.md

# The contract stamps are in the tup DAG (checks/Tupfile.lua); these build only them.
embedded-check:
	@$(TUP) checks/build/embedded.ok

freestanding-check:
	@$(TUP) checks/build/freestanding.ok

# Lean gate: no formatting, no examples, no boards. Needs `tup init` on a fresh clone.
ci: doxygen-style-check reference-check
	@$(TUP) tests checks/build/embedded.ok checks/build/freestanding.ok
	$(run_tests)
	@$(MAKE) --no-print-directory target-smoke

# ---------------------------------------------------------------------------
# Code hygiene (clang tools)
# ---------------------------------------------------------------------------
# These rewrite sources in place, which tup's input/output tracking forbids, so
# they live here like clang-format.
#
# clang on Windows targets MSVC by default. Point it at the build compiler's
# target and toolchain root (--sysroot) so it analyzes against the same
# libstdc++ regardless of PATH order. The constexpr budget matches .clangd:
# consteval designs in the tests exceed clang's default.
CLANG_TARGET   = $(shell $(GXX) -dumpmachine)
CLANG_SYSROOT  = $(shell $(GXX) -print-libgcc-file-name | sed 's,/lib/gcc/.*,,')
CLANG_ARGS     = --extra-arg=--target=$(CLANG_TARGET) \
                 --extra-arg=--sysroot=$(CLANG_SYSROOT) \
                 --extra-arg=-fbracket-depth=16384 \
                 --extra-arg=-fconstexpr-steps=100000000 \
                 --extra-arg=-Wno-unknown-warning-option

# clang-tidy --fix over every .cpp in compile_commands.json and the inc/ headers they use.
# (run-clang-tidy spells its option -extra-arg.)
tidy: compiledb
	@run-clang-tidy -p . -fix -header-filter='inc[/\\].*' \
		$(subst --extra-arg,-extra-arg,$(CLANG_ARGS)) '\.cpp$$'

# Include hygiene with clang-include-cleaner; IWYU keep/export pragmas apply.
# Scope: every inc/ header (analyzed on its own) plus tests/*.cpp. Excluded: the
# alternative math backends, which a profile swaps in for the default and so
# clash with it when parsed alone.
INCLUDE_CLEANER       = clang-include-cleaner -p . $(CLANG_ARGS)
INCLUDE_CLEANER_FILES = $(filter-out inc/damp/math/damp_backend.hpp inc/damp/math/series_backend.hpp, \
                            $(shell find inc -name '*.hpp')) \
                        $(wildcard tests/*.cpp)

# Print suggested include changes, one file per call (--print takes a single input).
iwyu: compiledb
	@for f in $(INCLUDE_CLEANER_FILES); do \
		$(INCLUDE_CLEANER) --print=changes "$$f" | sed "s|^|$$f: |"; \
	done

# Apply them, then re-sort with clang-format.
fix_includes: compiledb
	@$(INCLUDE_CLEANER) --edit $(INCLUDE_CLEANER_FILES)
	@$(MAKE) --no-print-directory format

# ---------------------------------------------------------------------------
# Targets (embedded)
# ---------------------------------------------------------------------------

# Host compile + run of product sketches against targets/damp_profile.hpp (ETL +
# Damp math); Arduino APIs are stubbed when ARDUINO is undefined. -Itargets comes
# first so config.hpp picks that profile. No PlatformIO needed.
SMOKE_SKETCHES = pid cart_pole
target-smoke:
	@mkdir -p targets/build
	@for s in $(SMOKE_SKETCHES); do \
		$(GXX) -std=c++20 -O3 -Wall -Wextra -Iinc -Itargets -Ilibs/etl/include \
			-Iexamples/control/$$s -o targets/build/$${s}_sketch_host.exe \
			examples/control/$$s/$${s}_sketch.cpp \
		&& ./targets/build/$${s}_sketch_host.exe || exit 1; \
	done
	@echo "target-smoke: host compile + run of $(SMOKE_SKETCHES) sketches (ETL profile) OK"

# PlatformIO board builds (needs `pio` on PATH). Default: every env in
# targets/platformio.ini. Subset: make targets ENVS="teensy41 nucleo_g474re"
ENVS ?=
targets:
	@cd targets && pio run -s $(foreach e,$(ENVS),-e $(e))
