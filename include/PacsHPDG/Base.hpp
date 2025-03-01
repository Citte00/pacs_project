/**
 * @file Base.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-11-28
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef INCLUDE_PACSHPDG_BASE_HPP
#define INCLUDE_PACSHPDG_BASE_HPP

#include <concepts>
#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <type_traits>
#include <omp.h>

// TOLERANCES.

// Zero tolerance.
#ifndef TOLERANCE
#define TOLERANCE 1E-20
#endif

// Geometry tolerance.
#ifndef GEOMETRY_TOLERANCE
#define GEOMETRY_TOLERANCE 1E-10
#endif

// Quadrature tolerance.
#ifndef QUADRATURE_TOLERANCE
#define QUADRATURE_TOLERANCE 5E-16
#endif

// Collapse (mesh_diagram) tolerance.
#ifndef COLLAPSE_TOLERANCE
#define COLLAPSE_TOLERANCE 1E-1
#endif

// Lloyd tolerance.
#ifndef LLOYD_TOLERANCE
#define LLOYD_TOLERANCE 1E-4
#endif


// CONSTANTS.

// Quadrature order.
#ifndef GAUSS_ORDER
#define GAUSS_ORDER 15
#endif

// Algebra iterations limit.
#ifndef ALGEBRA_ITER_MAX
#define ALGEBRA_ITER_MAX 25E3
#endif

#ifndef LLOYD_ITER_MAX
#define LLOYD_ITER_MAX 2E2
#endif

// Algebra m limit.
#ifndef ALGEBRA_M_MAX
#define ALGEBRA_M_MAX 25E1
#endif

// Tests limits.
#ifndef TESTS_MAX
#define TESTS_MAX 1E2
#endif

#ifndef DOFS_MAX
#define DOFS_MAX 3E4
#endif


// PARALLELISM.

// STL Parallelism.
#ifdef PARALLEL
#include <execution>
#include <tbb/tbb.h>
#define POLICY std::execution::par_unseq
#endif

// OpenMP Parallelism.
#ifdef _OPENMP
#include <omp.h>
#endif


// BASE TYPES, ALIASES AND CONCEPTS.

namespace pacs {

    /**
     * @brief Real alias.
     * 
     */
    using Real = long double;

    /**
     * @brief Function variadic template alias.
     * 
     * @tparam ReturnType Type to return from the function.
     * @tparam Params Parameters to pass to the function.
     */
    template <typename ReturnType, typename... Params>
    using Function = std::function<ReturnType(Params...)>;

    /**
     * @brief Mask alias.
     * 
     */
    using Mask = std::vector<bool>;

    // Numeric type concept for Matrix and Vector.

    /**
     * @brief Addable types.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Addable = requires(T first, T second) {
        {first += second} -> std::convertible_to<T>;
        {first -= second} -> std::convertible_to<T>;
    };

    /**
     * @brief Multipliable types.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Multipliable = requires(T first, T second) {
        {first *= second} -> std::convertible_to<T>;
        {first /= second} -> std::convertible_to<T>;
    };

    /**
     * @brief Absolute value supported types.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Absolute = requires(T value) {
        {std::abs(value)} -> std::convertible_to<Real>;
    };

    /**
     * @brief Numeric type.
     * 
     * @tparam T 
     */
    template<typename T>
    concept NumericType = Addable<T> && Multipliable<T> && Absolute<T>;

    /**
     * @brief Conjugable types.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Conjugable = requires(T value) {
        {std::conj(value)} -> std::convertible_to<T>;
    };
}

#endif