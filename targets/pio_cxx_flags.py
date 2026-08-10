# PlatformIO pre-script: C++-only flags must not land on C framework TUs.
# (cc1 warns: "-fno-rtti is valid for C++ but not for C" for every .c file.)
Import("env")  # type: ignore[name-defined]  # noqa: F821 — PlatformIO injects Import

env.Append(  # noqa: F821
    CXXFLAGS=[
        "-std=gnu++20",
        "-fno-rtti",
        "-fno-exceptions",
        "-fno-threadsafe-statics",
        "-fno-use-cxa-atexit",
        # C++20 deprecates volatile++/--; Teensy/Pico/ESP Arduino cores still use them.
        "-Wno-volatile",
    ]
)
