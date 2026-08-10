// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file npy_export.hpp
 * @brief Host-only NumPy .npy writer for dense simulation traces
 *
 * Dense hybrid / multi-rate runs (resonant tanks, high-rate IMU dumps) should
 * not go to Plotly HTML. Write planar or (N,C) float64 arrays here; load with
 * `numpy.load` in Python. Optional later: Python Plotly rendering of those files.
 *
 * Format: NumPy 1.0 `.npy` little-endian C-order, dtype float64 only (v0).
 *
 * @note Not part of `control.hpp`. No HDF5/MF4 dependency.
 */

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include "damp/simulation/hybrid.hpp"

namespace damp::sim {

/**
 * @brief Write a C-order float64 array as a .npy file
 *
 * @param path   Output path (e.g. "buck_trace.npy")
 * @param shape  Dimension sizes (e.g. {N, C} for N samples × C channels)
 * @param data   Row-major (C-order) length = product(shape)
 * @return true on success
 */
[[nodiscard]] inline bool write_npy_f64(
    std::string_view                path,
    const std::vector<std::size_t>& shape,
    const std::vector<double>&      data
) {
    if (shape.empty()) {
        return false;
    }
    std::size_t n = 1;
    for (std::size_t d : shape) {
        n *= d;
    }
    if (data.size() != n) {
        return false;
    }

    // Header dict, e.g. "{'descr': '<f8', 'fortran_order': False, 'shape': (N, C), }"
    std::string dict = "{'descr': '<f8', 'fortran_order': False, 'shape': (";
    for (std::size_t i = 0; i < shape.size(); ++i) {
        dict += std::to_string(shape[i]);
        if (i + 1 < shape.size()) {
            dict += ", ";
        } else if (shape.size() == 1) {
            dict += ","; // (N,) singleton tuple
        }
        // multi-dim: (N, C) no trailing comma after last — numpy accepts both
    }
    if (shape.size() > 1) {
        // remove trailing style: we produced (N, C) without trailing comma — good
    }
    dict += "), }";

    // Pad header so magic+version+hdr_len+dict is 64-byte aligned (npy 1.0)
    constexpr std::size_t preamble = 10; // magic(6)+ver(2)+hdrlen(2)
    std::size_t           pad = 16 - ((preamble + dict.size() + 1) % 16);
    if (pad == 16) {
        pad = 0;
    }
    dict.append(pad, ' ');
    dict.push_back('\n');

    const auto hdr_len = static_cast<std::uint16_t>(dict.size());

    std::ofstream out{std::string(path), std::ios::binary};
    if (!out) {
        return false;
    }
    const char magic[] = {'\x93', 'N', 'U', 'M', 'P', 'Y', 0x01, 0x00};
    out.write(magic, 8);
    out.put(static_cast<char>(hdr_len & 0xff));
    out.put(static_cast<char>((hdr_len >> 8) & 0xff));
    out.write(dict.data(), static_cast<std::streamsize>(dict.size()));
    out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(data.size() * sizeof(double)));
    return static_cast<bool>(out);
}

/**
 * @brief Pack HybridSimulationResult into (N, 2+NX+NU+NY) float64: t, mode, x..., u..., y...
 */
template<size_t NX, size_t NU, size_t NY, typename T>
[[nodiscard]] bool write_hybrid_trace_npy(
    std::string_view                             path,
    const HybridSimulationResult<NX, NU, NY, T>& tr
) {
    const std::size_t   N = tr.t.size();
    const std::size_t   C = 2 + NX + NU + NY;
    std::vector<double> data(N * C);
    for (std::size_t i = 0; i < N; ++i) {
        std::size_t c = 0;
        data[i * C + c++] = static_cast<double>(tr.t[i]);
        data[i * C + c++] = static_cast<double>(tr.mode[i]);
        for (std::size_t k = 0; k < NX; ++k) {
            data[i * C + c++] = static_cast<double>(tr.x[i](k));
        }
        for (std::size_t k = 0; k < NU; ++k) {
            data[i * C + c++] = static_cast<double>(tr.u[i](k));
        }
        for (std::size_t k = 0; k < NY; ++k) {
            data[i * C + c++] = static_cast<double>(tr.y[i](k));
        }
    }
    return write_npy_f64(path, {N, C}, data);
}

} // namespace damp::sim
