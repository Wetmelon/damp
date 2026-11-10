// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file core.hpp
 * @brief Fixed-size stack-allocated matrix type
 */

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "matrix_traits.hpp"

namespace damp {

/* Concrete (owning) types */
template<size_t Rows, size_t Cols, typename T>
struct Matrix;

template<size_t N, typename T>
struct ColVec;

template<size_t N, typename T>
struct RowVec;

/* Non-owning Views */
template<size_t N, typename T>
struct Diagonal;

template<size_t N, typename T>
struct UpperTriangle;

template<size_t N, typename T>
struct LowerTriangle;

template<size_t Rows, size_t Cols, size_t ParentCols, typename T>
struct Block;

template<size_t Rows, size_t Cols, typename T>
struct RowView;

template<size_t Rows, size_t Cols, typename T>
struct ColView;

template<size_t Rows, size_t Cols, typename T>
struct TransposeView;

namespace detail {

/**
 * @brief Inclusive address ranges. Relational `<` is defined inside one array, which
 * is the only case where two views can share an element. Different objects
 * never alias; a spurious overlap there only costs a snapshot.
 */
template<typename T>
[[nodiscard]] constexpr bool address_ranges_overlap(const T* a_first, const T* a_last, const T* b_first, const T* b_last) {
    const T* const a_lo = (a_last < a_first) ? a_last : a_first;
    const T* const a_hi = (a_last < a_first) ? a_first : a_last;
    const T* const b_lo = (b_last < b_first) ? b_last : b_first;
    const T* const b_hi = (b_last < b_first) ? b_first : b_last;
    return !(a_hi < b_lo) && !(b_hi < a_lo);
}

/**
 * @brief Element overlap of two row-major rectangles. Bounding spans reject the
 * disjoint case; the address walk rejects a strided gap that the span still covers.
 */
template<typename T>
[[nodiscard]] constexpr bool rects_share_element(
    const T* a, size_t a_rows, size_t a_cols, size_t a_stride,
    const T* b, size_t b_rows, size_t b_cols, size_t b_stride
) {
    if (a_rows == 0 || a_cols == 0 || b_rows == 0 || b_cols == 0) {
        return false;
    }
    if (damp::is_constant_evaluated()) {
        return true; // `<` across objects is not a constant expression; snapshot instead
    }
    const T* const a_last = a + ((a_rows - 1) * a_stride) + (a_cols - 1);
    const T* const b_last = b + ((b_rows - 1) * b_stride) + (b_cols - 1);
    if (!address_ranges_overlap(a, a_last, b, b_last)) {
        return false;
    }
    for (size_t i = 0; i < a_rows; ++i) {
        for (size_t j = 0; j < a_cols; ++j) {
            const T* const p = a + (i * a_stride) + j;
            for (size_t u = 0; u < b_rows; ++u) {
                const T* const row = b + (u * b_stride);
                for (size_t v = 0; v < b_cols; ++v) {
                    if (p == (row + v)) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

/**
 * @brief True when @p src's backing elements intersect [dest_first, dest_last].
 * Types without storage_first/storage_last are treated as aliasing so the
 * caller snapshots. A strided view's span includes the gaps between its rows,
 * so a non-overlapping column of the same parent can still look like an alias.
 */
template<typename T, typename Src>
[[nodiscard]] constexpr bool storage_overlaps(const T* dest_first, const T* dest_last, const Src& src) {
    if constexpr (requires { src.storage_first(); src.storage_last(); }) {
        using SrcElem = damp::remove_cv_t<damp::remove_pointer_t<decltype(src.storage_first())>>;
        if constexpr (!damp::is_same_v<damp::remove_cv_t<T>, SrcElem>) {
            return false;
        } else if (damp::is_constant_evaluated()) {
            return true; // `<` across objects is not a constant expression; snapshot instead
        } else {
            const T* const src_first = static_cast<const T*>(src.storage_first());
            const T* const src_last = static_cast<const T*>(src.storage_last());
            return address_ranges_overlap(dest_first, dest_last, src_first, src_last);
        }
    } else {
        return true;
    }
}

/**
 * @brief True when writing [dest_first, dest_last] element-wise from @p src may read
 * an element already overwritten. An owning (Matrix-derived) source can only
 * share storage by being the destination itself, element (r,c) onto (r,c), so
 * it never needs a snapshot.
 */
template<typename T, typename Src>
[[nodiscard]] constexpr bool may_alias(const T* dest_first, const T* dest_last, const Src& src) {
    if constexpr (damp::is_base_of_v<Matrix<Src::rows(), Src::cols(), typename Src::value_type>, Src>) {
        return false;
    } else {
        return storage_overlaps(dest_first, dest_last, src);
    }
}

/**
 * @brief Euclidean norm of n magnitudes mag(0) … mag(n−1), none of whose squares can overflow or underflow.
 *
 * A plain sum of squares overflows past √max (1.8e19 in float) and flushes below √min
 * (1.1e-19 in float, under flush-to-zero). Each magnitude is divided by d, the largest
 * magnitude clamped to [min, 1/min], so the scaled squares sum to less than 16·n and the
 * largest of them is not negligible. Two passes: the largest magnitude, then the scaled sum.
 *
 * The clamp keeps 1/d normal: under -freciprocal-math x/d becomes x·(1/d), and 1/amax is
 * subnormal past 1/min (flushed to zero under flush-to-zero) and overflows for subnormal
 * amax. The divisor is a run-time value because a constant scale s lets -fassociative-math
 * rewrite (x·s)² as x²·s², and g++ folds an s² below the subnormal range to 0. The result
 * is accurate to rounding for normal input; subnormal input keeps the precision it carries.
 *
 * A magnitude that is not finite (damp::isfinite, on its bits) is returned as the norm, so a
 * NaN or ±Inf entry gives a NaN or +Inf norm, as a plain sum of squares would and as LAPACK's
 * xNRM2 has since 3.10. damp::max drops a NaN against the running maximum, so without the
 * test an all-NaN vector, or a NaN among zeros, read as norm 0 and passed every ‖r‖ ≤ tol
 * certificate.
 *
 * @see "A portable Fortran program to find the Euclidean norm of a vector" (Blue, 1978), ACM TOMS 4(1)
 * @see "Algorithm 978: Safe Scaling in the Level 1 BLAS" (Anderson, 2017), ACM TOMS 44(1)
 *
 * @param n    Number of magnitudes
 * @param mag  Callable size_t → R returning |x_k| (≥ 0, or NaN / +Inf for a non-finite x_k)
 * @return √(Σ mag(k)²), or the first magnitude that is not finite
 */
template<typename R, typename Mag>
[[nodiscard]] constexpr R scaled_norm2(size_t n, Mag mag) {
    R amax = R{0};
    for (size_t k = 0; k < n; ++k) {
        const R a = mag(k);
        if (!damp::isfinite(a)) {
            return a;
        }
        amax = damp::max(amax, a);
    }
    if (!(amax > R{0})) {
        return R{0};
    }
    constexpr R lo = damp::numeric_limits<R>::min();
    const R     d = damp::clamp(amax, lo, R{1} / lo);
    R           ssq = R{0};
    for (size_t k = 0; k < n; ++k) {
        const R a = mag(k) / d;
        ssq += a * a;
    }
    return d * damp::sqrt(ssq);
}

/** @brief e with |x| = 2^e·f, f ∈ [1, 2); 0 for x = 0 or non-finite x. */
template<typename T>
[[nodiscard]] constexpr int binary_exponent(T x) {
    T v = damp::abs(x);
    if (!(v > T{0}) || !damp::isfinite(v)) {
        return 0;
    }
    constexpr int cap = damp::numeric_limits<T>::max_exponent - damp::numeric_limits<T>::min_exponent + damp::numeric_limits<T>::digits;
    int           e = 0;
    while (v >= T{2} && e < cap) {
        v *= static_cast<T>(0.5);
        ++e;
    }
    while (v < T{1} && e > -cap) {
        v *= T{2};
        --e;
    }
    return e;
}

/**
 * @brief Tag for a Matrix whose every entry is written before it is read: no zero-fill at run time
 *
 * Not API. The fill of a large matrix is a memset call the compiler cannot remove when a loop with
 * run-time indices does the writing, since it cannot prove the loop covers every entry. Constant
 * evaluation still starts from zeros, so nothing there reads an indeterminate value.
 */
struct Unfilled {
    explicit constexpr Unfilled() = default;
};

/** @brief The Unfilled tag. */
inline constexpr Unfilled unfilled{};

/** @brief Largest array any supported ABI may return in registers: four doubles on AArch64. */
inline constexpr size_t kReturnedInRegistersBytes = 32;

/** @brief Largest shape for_each_entry folds over: a local this size is worth keeping in registers. */
inline constexpr size_t kConstantIndexEntries = 64;

/** @brief One row of for_index_pairs: @p f<R, C>() for each column. */
template<size_t R, typename F, size_t... C>
[[gnu::always_inline]] constexpr void for_index_pairs_in_row(F& f, damp::index_sequence<C...> /*columns*/) {
    (f.template operator()<R, C>(), ...);
}

/**
 * @brief Call @p f<R, C>() for every pair of a Rows × Cols shape, the indices compile-time constants
 *
 * One Cols-wide fold per row, never a Rows·Cols pack: the arity stays within the compilers' default
 * limits and the frontend's cost stays low, as in Matrix's copy. Everything is forced inline, so @p f
 * reads the same data through no closure pointer and the compiler can keep it in registers.
 */
template<size_t Rows, size_t Cols, typename F>
[[gnu::always_inline]] constexpr void for_index_pairs(F&& f) {
    [&]<size_t... R>(damp::index_sequence<R...>) __attribute__((always_inline)) {
        (for_index_pairs_in_row<R>(f, damp::make_index_sequence<Cols>{}), ...);
    }(damp::make_index_sequence<Rows>{});
}

/**
 * @brief Call @p f(r, c) for every entry of a Rows × Cols shape: on a small shape at run time with
 *        constant indices, otherwise in a plain loop.
 *
 * GCC scalar-replaces a local only if every access has a constant index, and it stops unrolling a
 * loop past 16 iterations, so a loop over a ColVec<18> leaves that vector, and everything assigned
 * to it, in memory. for_index_pairs keeps the indices constant. Constant evaluation keeps the
 * loop: it costs less there.
 */
template<size_t Rows, size_t Cols, typename F>
[[gnu::always_inline]] constexpr void for_each_entry(F&& f) {
    if constexpr (Rows * Cols > 0 && Rows * Cols <= kConstantIndexEntries) {
        if (!damp::is_constant_evaluated()) {
            // Not forced as a whole: a dense operation inlines into its caller only where the compiler judges it worthwhile.
            [&] { for_index_pairs<Rows, Cols>([&]<size_t R, size_t C>() __attribute__((always_inline)) { f(R, C); }); }();
            return;
        }
    }
    for (size_t r = 0; r < Rows; ++r) {
        for (size_t c = 0; c < Cols; ++c) {
            f(r, c);
        }
    }
}

} // namespace detail

/**
 * @ingroup linear_algebra
 * @brief Fixed-size, stack-allocated matrix for linear algebra operations
 *
 * Provides efficient matrix arithmetic, multiplication, transposition, trace, and inversion.
 * Supports constexpr evaluation for compile-time matrix computations.
 *
 * @tparam Rows Number of rows
 * @tparam Cols Number of columns
 * @tparam T Element type — only @c float, @c double, or @c damp::complex of those
 *         (see @ref is_matrix_element_v). Nested @c Matrix / block types are
 *         not elements: use a flat matrix + @c .block<>() or multi-axis
 *         @c Lut1D/@c Lut2D/@c Lut3D for multi-parameter tables.
 */
template<size_t Rows, size_t Cols, typename T = double>
struct Matrix {
protected:
    damp::array<T, Rows * Cols> data_{};

    template<size_t, size_t, typename>
    friend struct Matrix;

    /** @brief One row of constant-index stores (fold arity = Cols, not Rows*Cols). */
    template<size_t R>
    constexpr void copy_row_from(const Matrix& other) {
        [&]<size_t... C>(damp::index_sequence<C...>) {
            ((data_[(R * Cols) + C] = other.data_[(R * Cols) + C]), ...);
        }(damp::make_index_sequence<Cols>{});
    }

    /**
     * @brief Element-wise copy via per-row folds
     *
     * Each store uses a compile-time constant index so GCC can SROA small local
     * matrices. Fold arity is Cols (one row), with an outer pack over Rows —
     * not a single pack of length Rows*Cols.
     */
    constexpr void copy_from(const Matrix& other) {
        [&]<size_t... R>(damp::index_sequence<R...>) {
            (copy_row_from<R>(other), ...);
        }(damp::make_index_sequence<Rows>{});
    }

    /*
     * Constructors build the elements in a local and hand it to the data_
     * mem-initializer. GCC rejects a constructor body that stores into data_
     * when the Matrix is a member of a const aggregate under constant
     * evaluation ("modifying a const object").
     */

    /** @brief Elements of any (r, c)-indexable source, cast to T. */
    template<typename Src>
    [[nodiscard]] static constexpr damp::array<T, Rows * Cols> elements_of(const Src& src) {
        damp::array<T, Rows * Cols> out = seed_elements();
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                out[(r * Cols) + c] = static_cast<T>(src(r, c));
            }
        }
        return out;
    }

    /**
     * @brief Elements of a nested initializer list; missing entries are zero. More rows or columns
     * than the matrix has do not compile when constant-evaluated and trap at run time.
     */
    [[nodiscard]] static constexpr damp::array<T, Rows * Cols> elements_of(
        std::initializer_list<std::initializer_list<T>> init
    ) {
        damp::array<T, Rows * Cols> out{};
        if (init.size() > Rows) {
            ::damp::detail::refuse_initializer_overflow();
        }
        size_t r = 0;
        for (const auto& row : init) {
            if (row.size() > Cols) {
                ::damp::detail::refuse_initializer_overflow();
            }
            size_t c = 0;
            for (const auto& val : row) {
                out[(r * Cols) + c] = val;
                ++c;
            }
            ++r;
        }
        return out;
    }

    /**
     * @brief Elements of a flat row-major list (vectors); missing entries are zero. More elements
     * than the vector has do not compile when constant-evaluated and trap at run time.
     */
    [[nodiscard]] static constexpr damp::array<T, Rows * Cols> flat_elements_of(std::initializer_list<T> init) {
        damp::array<T, Rows * Cols> out{};
        if (init.size() > Rows * Cols) {
            ::damp::detail::refuse_initializer_overflow();
        }
        size_t i = 0;
        for (const auto& val : init) {
            out[i] = val;
            ++i;
        }
        return out;
    }

    /** @brief Elements of a sized C array (vectors). */
    template<size_t M>
        requires(M == Rows * Cols)
    [[nodiscard]] static constexpr damp::array<T, Rows * Cols> flat_elements_of(const T (&arr)[M]) {
        damp::array<T, Rows * Cols> out{};
        for (size_t i = 0; i < M; ++i) {
            out[i] = arr[i];
        }
        return out;
    }

    /**
     * @brief Elements with no values yet: every entry must be written before it is read. Not constexpr,
     * so constant evaluation cannot reach it and GCC does not fold away the constructor that calls it
     * (on ARM a constructor returns `this`, and an empty one warns "value computed is not used").
     *
     * An array small enough to come back in registers would be read by the return itself, so it gets
     * zeros: the compiler removes that small a fill anyway. A larger one is built in the caller's storage.
     */
    [[nodiscard]] static damp::array<T, Rows * Cols> unset_elements() {
        damp::array<T, Rows * Cols> out; // NOLINT(cppcoreguidelines-pro-type-member-init): see detail::Unfilled
        if constexpr (sizeof(damp::array<T, Rows * Cols>) <= detail::kReturnedInRegistersBytes) {
            out = {};
        }
        return out; // one variable at function scope: GCC before 14 elides no other
    }

    /** @brief What a matrix about to be overwritten starts from: zeros under constant evaluation, else nothing. */
    [[nodiscard]] static constexpr damp::array<T, Rows * Cols> seed_elements() {
        if (damp::is_constant_evaluated()) {
            return {};
        }
        return unset_elements();
    }

    /**
     * @brief Copy-constructor seed: the whole array under constant evaluation; at run
     * time no values yet, for the constant-index folds that GCC can SROA to write.
     *
     * The unset array is built here rather than by unset_elements: GCC 14 crashes (an internal
     * compiler error in build_data_member_initialization) when a copy constructor reached from a
     * member's default initializer, in a local of a `constexpr` constructor, calls a function that is
     * not `constexpr`.
     */
    [[nodiscard]] static constexpr damp::array<T, Rows * Cols> copy_seed(const Matrix& other) {
        damp::array<T, Rows * Cols> out; // NOLINT(cppcoreguidelines-pro-type-member-init): copy_from writes every entry
        if (damp::is_constant_evaluated()) {
            out = other.data_;
        } else if constexpr (sizeof(damp::array<T, Rows * Cols>) <= detail::kReturnedInRegistersBytes) {
            out = {};
        }
        return out; // the only return, so built in the caller's storage: no copy reads the unset entries
    }

public:
    using value_type = T;

    static_assert(
        is_matrix_element_v<T>,
        "Matrix element type must be float, double, or damp::complex<float|double> "
        "(not nested Matrix or other types — use flat Matrix + block views, "
        "or Lut1D/Lut2D/Lut3D for multi-parameter lookup tables)"
    );

    /**
     * @brief Default constructor, initializes all elements to zero
     */
    constexpr Matrix() = default;

    /** @brief No zero-fill at run time. Not API: see detail::Unfilled. */
    constexpr explicit Matrix(const detail::Unfilled& /*tag*/) : data_(seed_elements()) {}

    /**
     * @brief Copy constructor and assignment operator
     *
     * Copies per row with a Cols-wide index_sequence fold (constant indices for
     * GCC SROA of small locals). Arity is Cols rather than Rows*Cols, which
     * keeps g++ frontend cost reasonable for typical matrix shapes.
     */
#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ < 12)
    // GCC 11 crashes (an internal compiler error in build_data_member_initialization) on a copy
    // constructor that branches on is_constant_evaluated, reached from a member's default
    // initializer in a local of a `constexpr` constructor, as RlsConfig's Bounds is. It copies the
    // array whole instead.
    constexpr Matrix(const Matrix& other) : data_(other.data_) {}
#else
    constexpr Matrix(const Matrix& other) : data_(copy_seed(other)) {
        if (!damp::is_constant_evaluated()) {
            copy_from(other);
        }
    }
#endif

    constexpr Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            copy_from(other);
        }
        return *this;
    }

    /**
     * @brief Move constructor and assignment operator
     */
    constexpr Matrix(Matrix&&) = default;
    constexpr Matrix& operator=(Matrix&&) = default;

    /**
     * @brief Destructor
     */
    constexpr ~Matrix() = default;

    /**
     * @brief Assignment from a matrix whose elements convert exactly (real into complex); a change
     * of precision is written `.as<U>()`
     */
    template<typename U>
        requires ConvertsExactly<U, T>
    constexpr Matrix& operator=(const Matrix<Rows, Cols, U>& other) {
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                data_[(r * Cols) + c] = static_cast<T>(other.data_[(r * Cols) + c]);
            }
        }
        return *this;
    }

    /**
     * @brief Explicit conversion from another element type; with `.as<U>()`, the only spelling of a
     * change in precision
     */
    template<typename U>
    explicit constexpr Matrix(const Matrix<Rows, Cols, U>& other) : data_(elements_of(other)) {}

    /**
     * @brief Constructor from nested initializer list; missing entries are zero. `A = {{…}}` assigns
     * through this constructor, so it zero-fills and refuses overflow the same way.
     */
    constexpr Matrix(std::initializer_list<std::initializer_list<T>> init) : data_(elements_of(init)) {}

    /**
     * @brief Constructor from damp::array, filled in row-major order
     */
    constexpr explicit Matrix(const damp::array<T, Rows * Cols>& arr) : data_(arr) {}

    /**
     * @brief A 1×1 matrix from its one element. Implicit, since a 1×1 matrix is its scalar: a SISO
     * weight or noise reads `.R = 1.0`.
     */
    constexpr Matrix(T value) // NOLINT(google-explicit-constructor): a 1×1 matrix is its scalar
        requires(Rows == 1 && Cols == 1)
        : data_{value} {}

    /**
     * @brief Returns a view over a sub-block of the matrix
     *
     * @param start_row Starting row index of the block
     * @param start_col Starting column index of the block
     * @tparam Brows Number of rows in the block
     * @tparam Bcols Number of columns in the block
     *
     * @return a Block view of the specified submatrix
     */
    template<size_t Brows, size_t Bcols>
    [[nodiscard]] constexpr Block<Brows, Bcols, Cols, T> block(size_t start_row, size_t start_col) & {
        static_assert(Brows <= Rows && Bcols <= Cols, "block larger than its parent");
        if constexpr (Brows == 0 || Bcols == 0) {
            return Block<Brows, Bcols, Cols, T>{nullptr, 0};
        } else {
            return Block<Brows, Bcols, Cols, T>{data_.data(), (start_row * Cols) + start_col};
        }
    }

    template<size_t Brows, size_t Bcols>
    [[nodiscard]] constexpr Block<Brows, Bcols, Cols, const T> block(size_t start_row, size_t start_col) const& {
        static_assert(Brows <= Rows && Bcols <= Cols, "block larger than its parent");
        if constexpr (Brows == 0 || Bcols == 0) {
            return Block<Brows, Bcols, Cols, const T>{nullptr, 0};
        } else {
            return Block<Brows, Bcols, Cols, const T>{data_.data(), (start_row * Cols) + start_col};
        }
    }

    /** @brief Owning copy. An rvalue's view would dangle at the end of the full expression. */
    template<size_t Brows, size_t Bcols>
    [[nodiscard]] constexpr Matrix<Brows, Bcols, T> block(size_t start_row, size_t start_col) && {
        return static_cast<const Matrix&>(*this).template block<Brows, Bcols>(start_row, start_col);
    }

    template<size_t Brows, size_t Bcols>
    [[nodiscard]] constexpr Matrix<Brows, Bcols, T> block(size_t start_row, size_t start_col) const&& {
        const auto view = static_cast<const Matrix&>(*this).template block<Brows, Bcols>(start_row, start_col);
        return Matrix<Brows, Bcols, T>{view};
    }

    /**
     * @brief Returns a view of the diagonal elements
     *
     * @return Diagonal view of the matrix
     */
    [[nodiscard]] constexpr Diagonal<Rows, T> diagonal() &
        requires(Rows == Cols)
    {
        return Diagonal<Rows, T>{*this};
    }

    [[nodiscard]] constexpr Diagonal<Rows, const T> diagonal() const&
        requires(Rows == Cols)
    {
        return Diagonal<Rows, const T>{*this};
    }

    /** @brief Owning diagonal matrix (off-diagonal zero). An rvalue's view would dangle at the end of the full expression. */
    [[nodiscard]] constexpr Matrix diagonal() const&&
        requires(Rows == Cols)
    {
        return Matrix{Diagonal<Rows, const T>{*this}};
    }

    /**
     * @brief Returns a view of the upper triangular elements
     *
     * @return UpperTriangle view of the matrix
     */
    [[nodiscard]] constexpr UpperTriangle<Rows, T> upper_triangle() const&
        requires(Rows == Cols)
    {
        return UpperTriangle<Rows, T>{*this};
    }

    /**
     * @brief Owning copy with the strict lower triangle zeroed. An rvalue's view would dangle at the
     * end of the full expression; mat::solve treats the copy as a general matrix.
     */
    [[nodiscard]] constexpr Matrix upper_triangle() const&&
        requires(Rows == Cols)
    {
        return Matrix{UpperTriangle<Rows, T>{*this}};
    }

    /**
     * @brief Returns a view of the lower triangular elements
     *
     * @return LowerTriangle view of the matrix
     */
    [[nodiscard]] constexpr LowerTriangle<Rows, T> lower_triangle() const&
        requires(Rows == Cols)
    {
        return LowerTriangle<Rows, T>{*this};
    }

    /**
     * @brief Owning copy with the strict upper triangle zeroed. An rvalue's view would dangle at the
     * end of the full expression; mat::solve treats the copy as a general matrix.
     */
    [[nodiscard]] constexpr Matrix lower_triangle() const&&
        requires(Rows == Cols)
    {
        return Matrix{LowerTriangle<Rows, T>{*this}};
    }

    /**
     * @brief Returns a view of the specified row
     *
     * @param row_index The index of the row to view
     * @return RowView of the specified row
     */
    [[nodiscard]] constexpr RowView<Rows, Cols, T> row(size_t row_index) & {
        return RowView<Rows, Cols, T>{*this, row_index};
    }

    [[nodiscard]] constexpr RowView<Rows, Cols, const T> row(size_t row_index) const& {
        return RowView<Rows, Cols, const T>{*this, row_index};
    }

    [[nodiscard]] constexpr Matrix<1, Cols, T> row(size_t row_index) && {
        return row_vector(row_index);
    }

    [[nodiscard]] constexpr Matrix<1, Cols, T> row(size_t row_index) const&& {
        return row_vector(row_index);
    }

    /**
     * @brief Returns a view of the specified column
     *
     * @param col_index The index of the column to view
     * @return ColView of the specified column
     */
    [[nodiscard]] constexpr ColView<Rows, Cols, T> col(size_t col_index) & {
        return ColView<Rows, Cols, T>{*this, col_index};
    }

    [[nodiscard]] constexpr ColView<Rows, Cols, const T> col(size_t col_index) const& {
        return ColView<Rows, Cols, const T>{*this, col_index};
    }

    [[nodiscard]] constexpr Matrix<Rows, 1, T> col(size_t col_index) && {
        return col_vector(col_index);
    }

    [[nodiscard]] constexpr Matrix<Rows, 1, T> col(size_t col_index) const&& {
        return col_vector(col_index);
    }

    /**
     * @brief Constructor from a compatible Block (copies data); explicit when the elements change
     * precision
     */
    template<size_t ParentCols, typename U>
    explicit(!ConvertsExactly<U, T>) constexpr Matrix(const Block<Rows, Cols, ParentCols, U>& block) : data_(elements_of(block)) {}

    /**
     * @brief Constructor from any MatrixLike type with matching dimensions (copies data); explicit
     * when the elements change precision
     */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == Cols && !damp::is_same_v<damp::remove_cvref_t<M>, Matrix<Rows, Cols, T>>)
    explicit(!ConvertsExactly<typename M::value_type, T>) constexpr Matrix(const M& other) : data_(elements_of(other)) {}

    /**
     * @brief Assignment from any MatrixLike type with matching dimensions whose elements convert exactly
     */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == Cols && ConvertsExactly<typename M::value_type, T>)
    constexpr Matrix& operator=(const M& other) {
        if constexpr (Rows == 0 || Cols == 0) {
            return *this;
        }
        // Owning the elements is one copy. A second copy is only for A = A.t()
        // and any other source that reads this storage while we write it.
        if (aliases(other)) {
            Matrix tmp{};
            for (size_t r = 0; r < Rows; ++r) {
                for (size_t c = 0; c < Cols; ++c) {
                    tmp.data_[(r * Cols) + c] = static_cast<T>(other(r, c));
                }
            }
            copy_from(tmp);
        } else {
            for (size_t r = 0; r < Rows; ++r) {
                for (size_t c = 0; c < Cols; ++c) {
                    data_[(r * Cols) + c] = static_cast<T>(other(r, c));
                }
            }
        }
        return *this;
    }

    /**
     * @brief Create a zero matrix
     */
    [[nodiscard]] static constexpr Matrix zeros() { return Matrix{}; }

    /**
     * @brief Create a matrix with all elements set to a constant value
     */
    [[nodiscard]] static constexpr Matrix constant(const T& value) {
        Matrix result;
        result.data_.fill(value);
        return result;
    }

    /**
     * @brief  Create an identity matrix (square matrix with ones on the diagonal)
     */
    [[nodiscard]] static constexpr Matrix identity()
        requires(Rows == Cols)
    {
        Matrix result{};
        for (size_t r = 0; r < Rows; ++r) {
            result.data_[(r * Cols) + r] = T{1};
        }
        return result;
    }

    /**
     * @brief Create a diagonal matrix from an array of diagonal elements
     *
     * @param diag Array containing the diagonal elements
     * @return Diagonal matrix with the specified diagonal elements
     */
    [[nodiscard]] static constexpr Matrix diagonal(const damp::array<T, Rows>& diag)
        requires(Rows == Cols)
    {
        Matrix result{};
        for (size_t r = 0; r < Rows; ++r) {
            result.data_[(r * Cols) + r] = diag[r];
        }
        return result;
    }

    /**
     * @brief Fill all elements with a constant value
     */
    constexpr void fill(const T& value) {
        data_.fill(value);
    }

    /**
     * @brief Type conversion method for converting to a different element type
     */
    template<typename U>
    [[nodiscard]] constexpr Matrix<Rows, Cols, U> as() const {
        return Matrix<Rows, Cols, U>(*this);
    }

    /**
     * @brief Get number of rows
     */
    [[nodiscard]] static constexpr size_t rows() { return Rows; }

    /**
     * @brief Get number of columns
     */
    [[nodiscard]] static constexpr size_t cols() { return Cols; }

    /**
     * @brief Get const pointer to data in row-major order
     */
    [[nodiscard]] constexpr const T* data() const { return data_.data(); }

    /** @brief First backing element. Pair with @ref storage_last for alias checks. */
    [[nodiscard]] constexpr const T* storage_first() const { return data_.data(); }

    /** @brief Last backing element (inclusive). Equal to @ref storage_first when empty. */
    [[nodiscard]] constexpr const T* storage_last() const {
        if constexpr (Rows == 0 || Cols == 0) {
            return data_.data();
        }
        return data_.data() + ((Rows * Cols) - 1);
    }

    /** @brief True when @p other is not this object but reads this storage (a view of it). */
    template<typename M>
    [[nodiscard]] constexpr bool aliases(const M& other) const {
        if constexpr (Rows == 0 || Cols == 0) {
            return false;
        } else {
            return detail::may_alias(storage_first(), storage_last(), other);
        }
    }

    /**
     * @brief Get pointer to data in row-major order
     */
    [[nodiscard]] constexpr T* data() { return data_.data(); }

    /**
     * @brief Iterator begin for range-based access
     */
    [[nodiscard]] constexpr const T* begin() const { return data(); }

    /**
     * @brief Iterator end for range-based access
     */
    [[nodiscard]] constexpr const T* end() const { return data() + (Rows * Cols); }

    /**
     * @brief Iterator begin for range-based access
     */
    [[nodiscard]] constexpr T* begin() { return data(); }

    /**
     * @brief Iterator end for range-based access
     */
    [[nodiscard]] constexpr T* end() { return data() + (Rows * Cols); }

    /**
     * @brief Element access operator
     * @param row Row index
     * @param col Column index
     * @return Value of element at (row, col)
     */
    constexpr T operator()(size_t row, size_t col) const { return data_[(row * Cols) + col]; }

    /**
     * @brief Element access operator
     * @param row Row index
     * @param col Column index
     * @return Reference to element at (row, col)
     */
    constexpr T& operator()(size_t row, size_t col) { return data_[(row * Cols) + col]; }

    /**
     * @brief Equality comparison operator
     */
    [[nodiscard]] constexpr bool operator==(const Matrix& other) const {
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                if (data_[(r * Cols) + c] != other.data_[(r * Cols) + c]) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Inequality comparison operator
     */
    [[nodiscard]] constexpr bool operator!=(const Matrix& other) const { return !(*this == other); }

    /**
     * @brief In-place addition with any MatrixLike type whose elements convert exactly
     */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == Cols && ConvertsExactly<typename M::value_type, T>)
    constexpr Matrix& operator+=(const M& other) {
        if (aliases(other)) {
            return *this += Matrix(other); // A += A.t(): read before writing
        }
        detail::for_each_entry<Rows, Cols>([&](size_t r, size_t c) { data_[(r * Cols) + c] += static_cast<T>(other(r, c)); });
        return *this;
    }

    /**
     * @brief In-place subtraction with any MatrixLike type whose elements convert exactly
     */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == Cols && ConvertsExactly<typename M::value_type, T>)
    constexpr Matrix& operator-=(const M& other) {
        if (aliases(other)) {
            return *this -= Matrix(other);
        }
        detail::for_each_entry<Rows, Cols>([&](size_t r, size_t c) { data_[(r * Cols) + c] -= static_cast<T>(other(r, c)); });
        return *this;
    }

    /**
     * @brief In-place addition with scalar
     */
    constexpr Matrix& operator+=(T scalar) {
        for (size_t i = 0; i < Rows * Cols; ++i) {
            data_[i] += scalar;
        }
        return *this;
    }

    /**
     * @brief In-place subtraction with scalar
     */
    constexpr Matrix& operator-=(T scalar) {
        for (size_t i = 0; i < Rows * Cols; ++i) {
            data_[i] -= scalar;
        }
        return *this;
    }

    /**
     * @brief In-place multiplication with scalar
     */
    constexpr Matrix& operator*=(T scalar) {
        for (size_t i = 0; i < Rows * Cols; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    /**
     * @brief In-place division with scalar
     */
    constexpr Matrix& operator/=(T scalar) {
        for (size_t i = 0; i < Rows * Cols; ++i) {
            data_[i] /= scalar;
        }
        return *this;
    }

    /**
     * @brief Matrix transpose (concrete copy)
     * @return Transposed matrix (Cols × Rows)
     * @see t() for a zero-copy TransposeView
     */
    [[nodiscard]] constexpr Matrix<Cols, Rows, T> transpose() const {
        Matrix<Cols, Rows, T> result{detail::unfilled};
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                result.data_[(c * Rows) + r] = data_[(r * Cols) + c];
            }
        }
        return result;
    }

    /**
     * @brief Non-owning transpose view (zero-copy)
     * @return TransposeView that swaps row/column indexing
     */
    [[nodiscard]] constexpr TransposeView<Rows, Cols, T> t() & {
        return TransposeView<Rows, Cols, T>{*this};
    }

    [[nodiscard]] constexpr TransposeView<Rows, Cols, const T> t() const& {
        return TransposeView<Rows, Cols, const T>{*this};
    }

    [[nodiscard]] constexpr Matrix<Cols, Rows, T> t() && { return transpose(); }

    [[nodiscard]] constexpr Matrix<Cols, Rows, T> t() const&& { return transpose(); }

    /**
     * @brief Conjugate transpose (Hermitian adjoint, concrete copy)
     *
     * Returns Aᴴ where Aᴴ(i,j) = conj(A(j,i)). For real types, identical
     * to transpose().
     *
     * @return Conjugate-transposed matrix (Cols × Rows)
     */
    [[nodiscard]] constexpr Matrix<Cols, Rows, T> conjugate_transpose() const {
        Matrix<Cols, Rows, T> result{detail::unfilled};
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                result.data_[(c * Rows) + r] = damp::conj(data_[(r * Cols) + c]);
            }
        }
        return result;
    }

    /**
     * @brief Reshape to new dimensions (total elements must match), taking elements in row-major order (NumPy C order)
     * @note Differs from MATLAB®'s reshape, which is column-major.
     * @tparam NewRows New number of rows
     * @tparam NewCols New number of columns
     * @return Reshaped matrix
     */
    template<size_t NewRows, size_t NewCols>
    [[nodiscard]] constexpr Matrix<NewRows, NewCols, T> reshape() const
        requires(NewRows * NewCols == Rows * Cols)
    {
        Matrix<NewRows, NewCols, T> result{detail::unfilled};
        for (size_t i = 0; i < Rows * Cols; ++i) {
            size_t r_new = i / NewCols;
            size_t c_new = i % NewCols;
            size_t r_old = i / Cols;
            size_t c_old = i % Cols;
            result(r_new, c_new) = (*this)(r_old, c_old);
        }
        return result;
    }

    /**
     * @brief Sum of diagonal elements (square matrices only)
     *
     * @return Trace of the matrix
     */
    [[nodiscard]] constexpr T trace() const
        requires(Rows == Cols)
    {
        T accum = T{0};
        for (size_t i = 0; i < Rows; ++i) {
            accum += data_[(i * Cols) + i];
        }
        return accum;
    }

    /**
     * @brief Frobenius norm (generalized "length") of the matrix
     *
     * The Frobenius norm is the square root of the sum of the absolute squares
     * of every element: ||A||_F = sqrt( sum |a_ij|^2 ).
     *
     * For a column vector this is the Euclidean length. For a matrix it is
     * a measure of overall magnitude, analogous to the Euclidean norm but
     * applied to all entries.
     *
     * Always returns a real value, even for complex-valued matrices. Entries are
     * divided by the largest magnitude before squaring (detail::scaled_norm2), so
     * no square overflows or underflows, under -ffast-math as well; the result
     * overflows only when the norm itself exceeds the largest finite value. A NaN
     * or ±Inf entry (in either part of a complex one) gives a norm that is not
     * finite, so a residual certificate cannot read an all-NaN residual as 0.
     *
     * @return Frobenius norm (real scalar)
     */
    [[nodiscard]] constexpr scalar_type_t<T> norm() const {
        if constexpr (is_complex_v<T>) {
            // Over the real and imaginary parts: |z| = hypot(re, im) reads (NaN, 0) as 0.
            return detail::scaled_norm2<scalar_type_t<T>>(2 * Rows * Cols, [&](size_t k) {
                const T& z = data_[k / 2];
                return damp::abs((k % 2 == 0) ? z.real() : z.imag());
            });
        } else {
            return detail::scaled_norm2<scalar_type_t<T>>(Rows * Cols, [&](size_t k) { return damp::abs(data_[k]); });
        }
    }

    /**
     * @brief Frobenius inner product with another matrix
     *
     * For real matrices:   dot(A,B) = sum( a_ij * b_ij )
     * For complex matrices: dot(A,B) = sum( conj(a_ij) * b_ij )  (=trace(A^H B))
     *
     * @param other Matrix to compute inner product with
     * @return Scalar inner product
     */
    [[nodiscard]] constexpr T dot(const Matrix& other) const {
        T result = T{0};
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                result += damp::conj(data_[(r * Cols) + c]) * other.data_[(r * Cols) + c];
            }
        }
        return result;
    }

    /**
     * @brief Sum of all elements in the matrix
     * @return Sum of all elements
     */
    [[nodiscard]] constexpr T sum() const {
        T result = T{0};
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                result += data_[(r * Cols) + c];
            }
        }
        return result;
    }

    /**
     * @brief Mean (average) of all elements in the matrix
     * @return Mean of all elements
     */
    [[nodiscard]] constexpr T mean() const {
        return sum() / static_cast<T>(Rows * Cols);
    }

    /**
     * @brief Matrix inverse A⁻¹ via mat::solve(A, I)
     *
     * Solves @f$ A X = I @f$ (Cholesky when Hermitian PD, else LU with partial
     * pivoting). Prefer mat::solve / @ref mat::cholesky_solve when the
     * inverse is only an intermediate (@f$ A^{-1} b @f$) — this method is for
     * when the inverse matrix itself is the deliverable.
     *
     * @return Inverse if nonsingular, damp::nullopt if singular
     * @note Compare with MATLAB®'s inv(A).
     * @see mat::solve()
     */
    [[nodiscard]] constexpr damp::optional<Matrix> inverse() const
        requires(Rows == Cols);

    /**
     * @brief Extract a row as a 1×Cols matrix
     * @param r_idx Row index
     * @return Row vector
     */
    [[nodiscard]] constexpr Matrix<1, Cols, T> row_vector(size_t r_idx) const {
        Matrix<1, Cols, T> result{detail::unfilled};
        for (size_t c = 0; c < Cols; ++c) {
            result.data_[c] = data_[(r_idx * Cols) + c];
        }
        return result;
    }

    /**
     * @brief Extract a column as a Rows×1 matrix
     * @param c_idx Column index
     * @return Column vector
     */
    [[nodiscard]] constexpr Matrix<Rows, 1, T> col_vector(size_t c_idx) const {
        Matrix<Rows, 1, T> result{detail::unfilled};
        for (size_t r = 0; r < Rows; ++r) {
            result.data_[r] = data_[(r * Cols) + c_idx];
        }
        return result;
    }

    /**
     * @brief Extract the first NewRows rows as a thin view (Block)
     * @tparam NewRows Number of rows to extract (must be <= Rows)
     * @return Block view with the first NewRows rows
     */
    template<size_t NewRows>
    [[nodiscard]] constexpr Block<NewRows, Cols, Cols, T> head() & {
        static_assert(NewRows <= Rows, "head<NewRows> called with NewRows > Rows");
        return Block<NewRows, Cols, Cols, T>(data_.data(), 0);
    }

    template<size_t NewRows>
    [[nodiscard]] constexpr Matrix<NewRows, Cols, T> head() && {
        return static_cast<const Matrix&&>(*this).template head<NewRows>();
    }

    /**
     * @brief Extract the first NewRows rows as a thin const view (Block)
     * @tparam NewRows Number of rows to extract (must be <= Rows)
     * @return Const Block view with the first NewRows rows
     */
    template<size_t NewRows>
    [[nodiscard]] constexpr Block<NewRows, Cols, Cols, const T> head() const& {
        static_assert(NewRows <= Rows, "head<NewRows> called with NewRows > Rows");
        return Block<NewRows, Cols, Cols, const T>(data_.data(), 0);
    }

    template<size_t NewRows>
    [[nodiscard]] constexpr Matrix<NewRows, Cols, T> head() const&& {
        static_assert(NewRows <= Rows, "head<NewRows> called with NewRows > Rows");
        return Matrix<NewRows, Cols, T>{Block<NewRows, Cols, Cols, const T>(data_.data(), 0)};
    }

    /**
     * @brief Extract the last NewRows rows as a thin view (Block)
     * @tparam NewRows Number of rows to extract (must be <= Rows)
     * @return Block view with the last NewRows rows
     */
    template<size_t NewRows>
    [[nodiscard]] constexpr Block<NewRows, Cols, Cols, T> tail() & {
        static_assert(NewRows <= Rows, "tail<NewRows> called with NewRows > Rows");
        return Block<NewRows, Cols, Cols, T>(data_.data() + ((Rows - NewRows) * Cols));
    }

    template<size_t NewRows>
    [[nodiscard]] constexpr Matrix<NewRows, Cols, T> tail() && {
        return static_cast<const Matrix&&>(*this).template tail<NewRows>();
    }

    /**
     * @brief Extract the last NewRows rows as a thin const view (Block)
     * @tparam NewRows Number of rows to extract (must be <= Rows)
     * @return Const Block view with the last NewRows rows
     */
    template<size_t NewRows>
    [[nodiscard]] constexpr Block<NewRows, Cols, Cols, const T> tail() const& {
        static_assert(NewRows <= Rows, "tail<NewRows> called with NewRows > Rows");
        return Block<NewRows, Cols, Cols, const T>(data_.data() + ((Rows - NewRows) * Cols));
    }

    template<size_t NewRows>
    [[nodiscard]] constexpr Matrix<NewRows, Cols, T> tail() const&& {
        static_assert(NewRows <= Rows, "tail<NewRows> called with NewRows > Rows");
        return Matrix<NewRows, Cols, T>{Block<NewRows, Cols, Cols, const T>(data_.data() + ((Rows - NewRows) * Cols))};
    }
};

// Common small-matrix aliases
template<typename T>
using Mat2 = Matrix<2, 2, T>;

template<typename T>
using Mat3 = Matrix<3, 3, T>;

template<typename T>
using Mat4 = Matrix<4, 4, T>;

template<typename T>
using Mat3x4 = Matrix<3, 4, T>;

template<typename T>
using Mat4x3 = Matrix<4, 3, T>;

// =============================================================================
// Free operators for MatrixLike types (Matrix, Block, RowView, ColView, etc.)
// =============================================================================

namespace detail {
/**
 * @brief Whether an element-wise free operator returns its operand's own type (ColVec, RowVec, …)
 * rather than Matrix. Only owning types qualify: a view built over the operator's local result
 * would dangle. Matrix itself does not: its result is already the type, and returns with no copy.
 */
template<typename A>
inline constexpr bool keeps_operand_type = !damp::is_same_v<A, Matrix<A::rows(), A::cols(), typename A::value_type>>
                                        && damp::is_base_of_v<Matrix<A::rows(), A::cols(), typename A::value_type>, A>
                                        && damp::is_constructible_v<A, Matrix<A::rows(), A::cols(), typename A::value_type>>;
} // namespace detail

/**
 * @brief Addition of two MatrixLike types of the same precision (with broadcasting support)
 */
template<MatrixLike A, MatrixLike B>
    requires SamePrecision<typename A::value_type, typename B::value_type>
[[nodiscard]] constexpr auto operator+(const A& a, const B& b) {
    constexpr size_t ResRows = damp::max(A::rows(), B::rows());
    constexpr size_t ResCols = damp::max(A::cols(), B::cols());
    static_assert(
        (A::rows() == 1 || A::rows() == ResRows) && (B::rows() == 1 || B::rows() == ResRows),
        "Incompatible row dimensions for broadcasting"
    );
    static_assert(
        (A::cols() == 1 || A::cols() == ResCols) && (B::cols() == 1 || B::cols() == ResCols),
        "Incompatible column dimensions for broadcasting"
    );
    using T = detail::common_element_t<typename A::value_type, typename B::value_type>;
    Matrix<ResRows, ResCols, T> result{detail::unfilled};
    detail::for_each_entry<ResRows, ResCols>([&](size_t i, size_t j) {
        const size_t ai = (A::rows() == 1) ? 0 : i;
        const size_t aj = (A::cols() == 1) ? 0 : j;
        const size_t bi = (B::rows() == 1) ? 0 : i;
        const size_t bj = (B::cols() == 1) ? 0 : j;
        result(i, j) = static_cast<T>(a(ai, aj)) + static_cast<T>(b(bi, bj));
    });
    if constexpr (damp::is_same_v<A, B> && detail::keeps_operand_type<A>) {
        return A(result);
    } else {
        return result;
    }
}

/**
 * @brief Subtraction of two MatrixLike types of the same precision (with broadcasting support)
 */
template<MatrixLike A, MatrixLike B>
    requires SamePrecision<typename A::value_type, typename B::value_type>
[[nodiscard]] constexpr auto operator-(const A& a, const B& b) {
    constexpr size_t ResRows = damp::max(A::rows(), B::rows());
    constexpr size_t ResCols = damp::max(A::cols(), B::cols());
    static_assert(
        (A::rows() == 1 || A::rows() == ResRows) && (B::rows() == 1 || B::rows() == ResRows),
        "Incompatible row dimensions for broadcasting"
    );
    static_assert(
        (A::cols() == 1 || A::cols() == ResCols) && (B::cols() == 1 || B::cols() == ResCols),
        "Incompatible column dimensions for broadcasting"
    );
    using T = detail::common_element_t<typename A::value_type, typename B::value_type>;
    Matrix<ResRows, ResCols, T> result{detail::unfilled};
    detail::for_each_entry<ResRows, ResCols>([&](size_t i, size_t j) {
        const size_t ai = (A::rows() == 1) ? 0 : i;
        const size_t aj = (A::cols() == 1) ? 0 : j;
        const size_t bi = (B::rows() == 1) ? 0 : i;
        const size_t bj = (B::cols() == 1) ? 0 : j;
        result(i, j) = static_cast<T>(a(ai, aj)) - static_cast<T>(b(bi, bj));
    });
    if constexpr (damp::is_same_v<A, B> && detail::keeps_operand_type<A>) {
        return A(result);
    } else {
        return result;
    }
}

/**
 * @brief Matrix multiplication of two MatrixLike types of the same precision
 */
template<MatrixLike A, MatrixLike B>
    requires(A::cols() == B::rows() && SamePrecision<typename A::value_type, typename B::value_type>)
[[nodiscard]] constexpr auto operator*(const A& a, const B& b) {
    using T = detail::common_element_t<typename A::value_type, typename B::value_type>;
    using Result = Matrix<A::rows(), B::cols(), T>;
    //! Each row is written from its first term, added to +0 as a zero-filled sum would be, so a sum
    //! of -0 terms is +0 as it always was. One result at function scope: GCC before 14 elides no other.
    Result result{detail::unfilled};
    if constexpr (A::cols() == 0) {
        result = Result{}; // an empty sum on every entry
    } else if constexpr (damp::is_base_of_v<Matrix<A::rows(), A::cols(), T>, A> && damp::is_base_of_v<Matrix<B::rows(), B::cols(), T>, B> && A::rows() * A::cols() * B::cols() > 0) {
        // Owning operands: walk the storage directly. Same arithmetic as the view
        // path; it only spares constant evaluation the per-element accessor calls,
        // and there skips the zero entries that block-structured design matrices
        // are full of, adding to the zeros that evaluation starts from.
        const T* const ap = a.data();
        const T* const bp = b.data();
        T* const       rp = result.data();
        for (size_t i = 0; i < A::rows(); ++i) {
            T* const rrow = rp + (i * B::cols());
            size_t   first = 0;
            if (!damp::is_constant_evaluated()) {
                const T ai0 = ap[i * A::cols()];
                for (size_t j = 0; j < B::cols(); ++j) {
                    rrow[j] = T{} + (ai0 * bp[j]);
                }
                first = 1;
            }
            for (size_t k = first; k < A::cols(); ++k) {
                const T aik = ap[(i * A::cols()) + k];
                if (damp::is_constant_evaluated() && aik == T{0}) {
                    continue;
                }
                const T* const brow = bp + (k * B::cols());
                for (size_t j = 0; j < B::cols(); ++j) {
                    rrow[j] += aik * brow[j];
                }
            }
        }
    } else {
        for (size_t i = 0; i < A::rows(); ++i) {
            const T ai0 = static_cast<T>(a(i, 0));
            for (size_t j = 0; j < B::cols(); ++j) {
                result(i, j) = T{} + (ai0 * static_cast<T>(b(0, j)));
            }
            for (size_t k = 1; k < A::cols(); ++k) {
                const T aik = static_cast<T>(a(i, k));
                for (size_t j = 0; j < B::cols(); ++j) {
                    result(i, j) += aik * static_cast<T>(b(k, j));
                }
            }
        }
    }
    return result;
}

/**
 * @brief Symmetric congruence (quadratic) form  S = M X Mᵀ
 *
 * Computes the congruence transform of a symmetric matrix X by M. Because X is
 * symmetric the result is symmetric, so only the lower triangle is evaluated and
 * mirrored: the output is symmetric to the last bit (no rounding-induced
 * asymmetry) and the outer product costs ~half a general triple product. This is
 * the standard building block for covariance propagation — A P Aᵀ, C P Cᵀ, the
 * Joseph-form update — where keeping the covariance exactly symmetric is what
 * preserves its positive-definiteness across iterations.
 *
 * Precondition: X is symmetric. If it is not, the result is built from the lower
 * triangle of M X Mᵀ and is therefore not the true product. (For complex T this
 * computes M X Mᵀ with the plain transpose, not the Hermitian adjoint.)
 *
 * @tparam Mat MatrixLike, p×n
 * @tparam Sym MatrixLike, n×n (symmetric)
 * @return p×p symmetric Matrix  M X Mᵀ
 */
template<MatrixLike Mat, MatrixLike Sym>
    requires(Mat::cols() == Sym::rows() && Sym::rows() == Sym::cols() && SamePrecision<typename Mat::value_type, typename Sym::value_type>)
[[nodiscard]] constexpr auto quadratic_form(const Mat& m, const Sym& x) {
    using T = detail::common_element_t<typename Mat::value_type, typename Sym::value_type>;
    constexpr size_t P = Mat::rows();
    constexpr size_t N = Mat::cols();
    const auto       mx = m * x; // p×n
    Matrix<P, P, T>  result{detail::unfilled};
    for (size_t i = 0; i < P; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            T acc{};
            for (size_t k = 0; k < N; ++k) {
                acc += static_cast<T>(mx(i, k)) * static_cast<T>(m(j, k));
            }
            result(i, j) = acc;
            result(j, i) = acc; // mirror → exact symmetry
        }
    }
    return result;
}

/**
 * @brief Unary negation for any MatrixLike type
 */
template<MatrixLike A>
[[nodiscard]] constexpr auto operator-(const A& a) {
    using T = A::value_type;
    Matrix<A::rows(), A::cols(), T> result{detail::unfilled};
    detail::for_each_entry<A::rows(), A::cols()>([&](size_t i, size_t j) { result(i, j) = -static_cast<T>(a(i, j)); });
    if constexpr (detail::keeps_operand_type<A>) {
        return A(result);
    } else {
        return result;
    }
}

// MatrixLike + Scalar
template<MatrixLike A, typename Scalar>
    requires(damp::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator+(const A& a, Scalar scalar) {
    using T = A::value_type;
    Matrix<A::rows(), A::cols(), T> result{detail::unfilled};
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t j = 0; j < A::cols(); ++j) {
            result(i, j) = static_cast<T>(a(i, j)) + static_cast<T>(scalar);
        }
    }
    if constexpr (detail::keeps_operand_type<A>) {
        return A(result);
    } else {
        return result;
    }
}

// MatrixLike - Scalar
template<MatrixLike A, typename Scalar>
    requires(damp::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator-(const A& a, Scalar scalar) {
    using T = A::value_type;
    Matrix<A::rows(), A::cols(), T> result{detail::unfilled};
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t j = 0; j < A::cols(); ++j) {
            result(i, j) = static_cast<T>(a(i, j)) - static_cast<T>(scalar);
        }
    }
    if constexpr (detail::keeps_operand_type<A>) {
        return A(result);
    } else {
        return result;
    }
}

// Scalar - MatrixLike
template<typename Scalar, MatrixLike A>
    requires(damp::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator-(Scalar scalar, const A& a) {
    using T = A::value_type;
    Matrix<A::rows(), A::cols(), T> result{detail::unfilled};
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t j = 0; j < A::cols(); ++j) {
            result(i, j) = static_cast<T>(scalar) - static_cast<T>(a(i, j));
        }
    }
    if constexpr (detail::keeps_operand_type<A>) {
        return A(result);
    } else {
        return result;
    }
}

// MatrixLike * Scalar
template<MatrixLike A, typename Scalar>
    requires(damp::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator*(const A& a, Scalar scalar) {
    using T = A::value_type;
    Matrix<A::rows(), A::cols(), T> result{detail::unfilled};
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t j = 0; j < A::cols(); ++j) {
            result(i, j) = static_cast<T>(a(i, j)) * static_cast<T>(scalar);
        }
    }
    if constexpr (detail::keeps_operand_type<A>) {
        return A(result);
    } else {
        return result;
    }
}

// Scalar * MatrixLike
template<typename Scalar, MatrixLike A>
    requires(damp::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator*(Scalar scalar, const A& a) {
    return a * scalar;
}

// MatrixLike / Scalar
template<MatrixLike A, typename Scalar>
    requires(damp::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator/(const A& a, Scalar scalar) {
    using T = A::value_type;
    Matrix<A::rows(), A::cols(), T> result{detail::unfilled};
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t j = 0; j < A::cols(); ++j) {
            result(i, j) = static_cast<T>(a(i, j)) / static_cast<T>(scalar);
        }
    }
    if constexpr (detail::keeps_operand_type<A>) {
        return A(result);
    } else {
        return result;
    }
}

namespace mat {

/**
 * @brief Apply @p f to every entry, B(i, j) = f(A(i, j)): the matrix counterpart of std::transform
 *
 * The result's element type is what @p f returns, so the magnitudes of a complex matrix
 * come back real; it must be a Matrix element type (float, double or damp::complex).
 * @code
 * const auto absA = mat::transform(A, [](double v) { return damp::abs(v); }); // |A|, entrywise
 * @endcode
 *
 * @param A Matrix
 * @param f Callable on one entry
 * @return f applied to each entry of A
 */
template<size_t R, size_t C, typename T, typename F>
[[nodiscard]] constexpr auto transform(const Matrix<R, C, T>& A, F f) {
    using U = damp::remove_cvref_t<decltype(f(A(0, 0)))>;
    Matrix<R, C, U> out{};
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            out(i, j) = f(A(i, j));
        }
    }
    return out;
}

/**
 * @brief Whether @p p holds for every entry of @p A: the matrix counterpart of std::all_of
 *
 * Row by row, stopping at the first entry that fails.
 * @code
 * const bool finite = mat::all_of(x, [](double v) { return damp::isfinite(v); });
 * @endcode
 *
 * @param A Matrix
 * @param p Predicate on one entry
 * @return true when p(A(i, j)) holds for every i, j
 */
template<size_t R, size_t C, typename T, typename P>
[[nodiscard]] constexpr bool all_of(const Matrix<R, C, T>& A, P p) {
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            if (!p(A(i, j))) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Whether @p p holds for some entry of @p A: the matrix counterpart of std::any_of
 *
 * Row by row, stopping at the first entry that passes.
 * @code
 * const bool saturated = mat::any_of(u, [&](double v) { return damp::abs(v) >= u_max; });
 * @endcode
 *
 * @param A Matrix
 * @param p Predicate on one entry
 * @return true when p(A(i, j)) holds for at least one i, j
 */
template<size_t R, size_t C, typename T, typename P>
[[nodiscard]] constexpr bool any_of(const Matrix<R, C, T>& A, P p) {
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            if (p(A(i, j))) {
                return true;
            }
        }
    }
    return false;
}

} // namespace mat

} // namespace damp
