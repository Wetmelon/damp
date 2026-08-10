# Target compile-smoke (PlatformIO®)

Design Is Deploy means the same C++ types run on MCU. These envs compile a
product sketch from `examples/` against real Arduino® cores (C++20 required).

There is no parallel `targets/sketches/` tree — board smoke uses the same
`*_sketch.cpp` files you fork for product work.

## Requirements

- [PlatformIO® Core](https://docs.platformio.org/en/latest/core/installation.html) (`pio` on PATH)
- C++20 — do not drop to C++17 (concepts, `constinit`)
- ETL submodule at `libs/etl` (default target profile uses `DAMP_BACKEND_ETL`)
- Runtime-oriented compile flags (see `platformio.ini`): `-O3`, `-ffast-math`,
  `-fno-rtti`, `-fno-exceptions`, `-fomit-frame-pointer`, section GC, `-DNDEBUG`.
  Framework `-Os`/`-O2` are unflagged so our `-O3` wins. No LTO (header-only
  library + single sketch TU; LTO is framework/toolchain friction for little gain).
- First build per platform downloads toolchains (can take several minutes)

## Default profile

[`damp_profile.hpp`](damp_profile.hpp) (on the PIO include path via `-I.`):

| Macro | Effect |
| ----- | ------ |
| `DAMP_BACKEND_ETL` | `damp::array` / `optional` / … → `etl::` |
| `DAMP_MATH_BACKEND_DAMP` | runtime float math via `damp/math/trig.hpp` |

Host SILs under `examples/` keep stdlib + `std::` math (their own `damp_profile.hpp`).

## Board matrix (verified 2026-07)

| Env | Board | Status | Notes |
|-----|--------|--------|--------|
| `teensy41` | Teensy 4.1 | Green | Default. GCC 15 (Teensy package) |
| `teensy40` | Teensy 4.0 | Green | Same family as 4.1 |
| `nucleo_g474re` | Nucleo G474RE | Green | STM32G4, GCC 12.3 (`ststm32`) |
| `blackpill_f411ce` | Black Pill F411 | Green | Cheap STM32F4 lab board |
| `pico` | Raspberry Pi Pico | Green | earlephilhower core via [maxgerhardt/platform-raspberrypi](https://github.com/maxgerhardt/platform-raspberrypi) (stock PIO RP2040 is GCC 9 → no `gnu++20`) |
| `esp32s3` | ESP32-S3-DevKitC-1 | Green | [pioarduino](https://github.com/pioarduino/platform-espressif32) 53.03.13 + GCC 13 Xtensa (stock `espressif32` is GCC 8) |
| Uno R4 Minima | — | Blocked | PIO `renesas-ra` still GCC 7.2 |

Primary pair for CI-ish smoke: `teensy41` + `nucleo_g474re`.

## Which sketch is built?

`platformio.ini` points `src_dir` at a product folder and builds only the
`*_sketch.cpp` (not the SIL):

| Default | Product folder | Source |
| ------- | -------------- | ------ |
| PID | [`examples/control/pid/`](../examples/control/pid/) | `pid_sketch.cpp` |
| LQR (optional) | [`examples/control/cart_pole/`](../examples/control/cart_pole/) | set `src_dir` + `build_src_filter` as commented in `platformio.ini` |

## Commands

```bash
cd targets

# Full green matrix (same as make targets) — PID product sketch
pio run

# Primary smoke pair only
pio run -e teensy41 -e nucleo_g474re

# Flash (board plugged in)
pio run -e teensy41 -t upload
```

From repo root:

```bash
make targets                                    # all green envs (pio -s; quiet compile lines)
make targets ENVS="teensy41 nucleo_g474re"      # primary pair only
make target-smoke                               # host g++ of pid + cart_pole (ETL profile)
# Verbose compile lines: cd targets && pio run -e teensy41
```

## Toolchain notes

- Stock platforms are often too old. Pico and ESP32 need the forks/pins in `platformio.ini` — do not “simplify” back to bare `platform = raspberrypi` / `espressif32` without re-verifying C++20.
- ESP32 first install can fail mid-postinstall on some Python setups; re-run `pio run -e esp32s3` once packages are cached (often succeeds on the second try).
- TI AM263x / modern F28 are product SDK targets (`tiarmclang`), not PIO envs here.

## Relation to host build

| Host (`make` / tup) | Target (`pio` / `target-smoke`) |
|---------------------|----------------------------------|
| tests, SIL, REFERENCE | real MCU cores or host stub of the same `*_sketch.cpp` |
| `examples/` damp_profile (stdlib) | `targets/damp_profile.hpp` (ETL + Damp math) |
| `embedded-check` / `freestanding-check` | flashable artifact |

Both matter; neither replaces the other. Domain spines: [docs/spines/](../docs/spines/).
