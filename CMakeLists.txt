cmake_minimum_required(VERSION 3.20)
project(damp
    VERSION 0.1.0
    LANGUAGES CXX
    DESCRIPTION "Header-only C++20 control systems library for embedded systems"
)

# Packaging only: INTERFACE include path + C++20 for find_package / FetchContent
# consumers. Day-to-day develop/test/examples use make + tup (see README).
# Target MCU smoke uses PlatformIO® under targets/. Do not grow a second test suite
# here unless CI explicitly needs it.

# --- Header-only library target ---
add_library(damp INTERFACE)

target_include_directories(damp INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(damp INTERFACE cxx_std_20)

# --- Install rules ---
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS damp
    EXPORT dampTargets
)

install(DIRECTORY inc/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)

install(EXPORT dampTargets
    FILE dampTargets.cmake
    NAMESPACE wetmelon::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/damp
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/dampConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/dampConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/damp
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/dampConfigVersion.cmake
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/dampConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/dampConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/damp
)

# --- Optional: tests & examples ---
option(DAMP_BUILD_TESTS    "Build unit tests"  OFF)
option(DAMP_BUILD_EXAMPLES "Build examples"    OFF)

if(DAMP_BUILD_TESTS OR DAMP_BUILD_EXAMPLES)
    add_library(fmt STATIC
        libs/fmt/src/format.cc
        libs/fmt/src/os.cc
    )
    target_include_directories(fmt PUBLIC libs/fmt/include)
    target_compile_features(fmt PUBLIC cxx_std_20)
endif()

# Shared warning/optimization flags (mirrors Tuprules.lua WARNINGS + CXXFLAGS).
# FAST_MATH adds -ffast-math, which the tup build applies to test TUs only --
# constexpr static_asserts act as a tripwire that compile-time eval stays
# IEEE-strict even under -ffast-math.
function(damp_apply_flags target fast_math)
    target_compile_features(${target} PRIVATE cxx_std_20)
    # Mirrors Tuprules.lua WARNINGS (Jason Turner / cppbestpractices set).
    target_compile_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:
            -O3 -march=native
            -ffunction-sections -fdata-sections
            -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align
            -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion
            -Wnull-dereference -Wdouble-promotion -Wfloat-conversion -Wformat=2
            -Wimplicit-fallthrough
        >
        $<$<CXX_COMPILER_ID:GNU>:
            -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches
            -Wlogical-op -Wuseless-cast -Wno-psabi
        >
        $<$<CXX_COMPILER_ID:Clang>:-Wno-nan-infinity-disabled>
        $<$<CXX_COMPILER_ID:MSVC>:/W4 /O2>
    )
    target_link_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wl,--gc-sections>)
    # -ffast-math + -static apply to the test runners; examples build IEEE-strict
    # and link shared, matching Tuprules.lua (TEST_CXXFLAGS / LDFLAGS).
    if(fast_math)
        target_compile_options(${target} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-ffast-math>)
        target_link_options(${target} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-static>)
    endif()
endfunction()

# --- Tests ---
if(DAMP_BUILD_TESTS)
    enable_testing()

    # Default profile: every tests/*.cpp linked into one runner (matches tup).
    file(GLOB TEST_SOURCES tests/*.cpp)
    add_executable(test_runner ${TEST_SOURCES})
    target_link_libraries(test_runner PRIVATE damp fmt)
    target_include_directories(test_runner PRIVATE
        tests
        libs
        libs/plotlypp/include
        libs/json/single_include
    )
    damp_apply_flags(test_runner TRUE)
    add_test(NAME damp_tests COMMAND test_runner)

    # The Damp-math and ETL container backends are ODR-incompatible with the
    # default std backend the runner links against, so each gets its own
    # executable (mirrors tests/{damp,etl}_backend/Tupfile.lua). No fmt; the ETL
    # variant's own include dir must precede tests/ so its damp_profile.hpp wins.
    file(GLOB DAMP_BACKEND_SOURCES tests/damp_backend/*.cpp)
    add_executable(test_damp_backend ${DAMP_BACKEND_SOURCES})
    target_link_libraries(test_damp_backend PRIVATE damp)
    target_include_directories(test_damp_backend PRIVATE tests/damp_backend tests)
    damp_apply_flags(test_damp_backend TRUE)
    add_test(NAME damp_tests_damp_backend COMMAND test_damp_backend)

    file(GLOB ETL_BACKEND_SOURCES tests/etl_backend/*.cpp)
    add_executable(test_etl_backend ${ETL_BACKEND_SOURCES})
    target_link_libraries(test_etl_backend PRIVATE damp)
    target_include_directories(test_etl_backend PRIVATE
        tests/etl_backend tests libs/etl/include)
    damp_apply_flags(test_etl_backend TRUE)
    add_test(NAME damp_tests_etl_backend COMMAND test_etl_backend)
endif()

# --- Examples ---
if(DAMP_BUILD_EXAMPLES)
    file(GLOB EXAMPLE_SOURCES examples/*.cpp)
    foreach(src ${EXAMPLE_SOURCES})
        get_filename_component(name ${src} NAME_WE)
        add_executable(${name} ${src})
        target_link_libraries(${name} PRIVATE damp fmt)
        target_include_directories(${name} PRIVATE
            examples
            libs
            libs/plotlypp/include
            libs/json/single_include
        )
        damp_apply_flags(${name} FALSE)  # examples build IEEE-strict (no -ffast-math)
    endforeach()
endif()

# --- Docs (mirrors `make docs`) ---
find_package(Doxygen)
if(DOXYGEN_FOUND AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
    add_custom_target(docs
        COMMAND ${CMAKE_COMMAND} -E make_directory docs/html
        COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()
