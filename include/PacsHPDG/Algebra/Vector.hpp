/**
 * @file Vector.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Vector class and methods.
 * @date 2025-01-18
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#ifndef INCLUDE_PACSHPDG_ALGEBRA_VECTOR_HPP
#define INCLUDE_PACSHPDG_ALGEBRA_VECTOR_HPP

#include "../Base.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace pacs {

/**
 * @brief Vector structure.
 *
 * @tparam T
 */
template <NumericType T> struct Vector {

  // Shape.
  std::size_t length;

  // Elements.
  std::vector<T> elements;

  // CONSTRUCTORS.

  /**
   * @brief Constructs a new empty Vector.
   *
   * @param length Length.
   */
  Vector(const std::size_t &length)
      : length{length}, elements(length, static_cast<T>(0)) {
#ifndef NDEBUG // Integrity check.
    assert(length > 0);
#endif
  }

  /**
   * @brief Constructs a new homogeneous Vector.
   *
   * @param length Length.
   * @param value Scalar.
   */
  Vector(const std::size_t &length, const T &value)
      : length{length}, elements(length, value) {
#ifndef NDEBUG // Integrity check.
    assert(length > 0);
#endif
  }

  /**
   * @brief Constructs a new Vector from a given std::vector.
   *
   * @param length Length.
   * @param elements Elements.
   */
  Vector(const std::size_t &length, const std::vector<T> &elements)
      : length{length}, elements(elements.begin(), elements.end()) {
#ifndef NDEBUG // Integrity check.
    assert(length > 0);
    assert(elements.size() == length);
#endif
  }

  /**
   * @brief Copy constructor.
   *
   * @param vector Vector.
   */
  Vector(const Vector &vector)
      : length{vector.length},
        elements(vector.elements.begin(), vector.elements.end()) {}

  /**
   * @brief Copy operator.
   *
   * @param vector Vector.
   * @return Vector&
   */
  Vector &operator=(const Vector &vector) {
#ifndef NDEBUG // Integrity check.
    assert(this->length == vector.length);
#endif

#ifdef PARALLEL
    std::copy(POLICY, vector.elements.begin(), vector.elements.end(),
              this->elements.begin());
#else
    std::copy(vector.elements.begin(), vector.elements.end(),
              this->elements.begin());
#endif

    return *this;
  }

  /**
   * @brief Scalar copy operator.
   *
   * @param scalar Scalar.
   * @return Vector&
   */
  Vector &operator=(const T &scalar) {
#ifdef PARALLEL
    std::for_each(POLICY, this->elements.begin(), this->elements.end(),
                  [scalar](auto &element) { element = scalar; });
#else
    std::for_each(this->elements.begin(), this->elements.end(),
                  [scalar](auto &element) { element = scalar; });
#endif

    return *this;
  }

  // CONVERSION.

  /**
   * @brief Converts the Vector into a std::vector<T>.
   *
   * @return std::vector<T>
   */
  operator std::vector<T>() const { return this->elements; }

  // RESIZE.

  /**
   * @brief Resizing the vector length keeping the elements.
   *
   * @param new_len New length.
   * @return Vector
   */
  void resize(const std::size_t &new_len) {
    this->elements.resize(new_len, 0.0);
    this->length = new_len;
  }

  // READ AND WRITE.

  /**
   * @brief Const subscript operator, returns the j-th element.
   *
   * @param j Index.
   * @return T
   */
  inline T operator[](const std::size_t &j) const {
#ifndef NDEBUG // Integrity check.
    assert(j < this->length);
#endif

    return this->elements[j];
  }

  /**
   * @brief Subscript operator, returns a reference to the j-th element.
   *
   * @param j Index.
   * @return T&
   */
  inline T &operator[](const std::size_t &j) {
#ifndef NDEBUG // Integrity check.
    assert(j < this->length);
#endif

    return this->elements[j];
  }

  /**
   * @brief Returns the [j, k) range.
   *
   * @param j Index.
   * @param k
   * @return Vector
   */
  Vector operator()(const std::size_t &j, const std::size_t &k) const {
#ifndef NDEBUG // Integrity check.
    assert(j != k);
    assert((j < this->length) && (k < this->length + 1));
#endif

    Vector result{(j < k) ? k - j : j - k};

    if (j < k) {
      for (std::size_t h = j; h < k; ++h)
        result[h - j] = this->elements[h];
    } else {
      for (std::size_t h = k; h > j; h--)
        result[h - k] = this->elements[h];
    }

    return result;
  }

  /**
   * @brief Returns the [j, end) range.
   *
   * @param j Index.
   * @return Vector
   */
  Vector operator()(const std::size_t &j) const {
#ifndef NDEBUG // Integrity check.
    assert(j < this->length);
#endif

    Vector result{this->length - j};

    for (std::size_t h = j; h < this->length; ++h)
      result[h - j] = this->elements[h];

    return result;
  }

  /**
   * @brief Returns a sub-Vector given a std::vector of indices.
   *
   * @param indices Indices.
   * @return Vector
   */
  Vector operator()(const std::vector<std::size_t> &indices) const {
#ifndef NDEBUG // Integrity check.
    for (const auto &index : indices)
      assert(index < this->length);
#endif

    Vector result{indices.size()};

    for (std::size_t j = 0; j < indices.size(); ++j)
      result[j] = this->elements[indices[j]];

    return result;
  }

  /**
   * @brief Returns a sub-Vector given a Mask.
   *
   * @param mask Mask.
   * @return Vector
   */
  Vector operator()(const Mask &mask) const {
#ifndef NDEBUG // Integrity check.
    assert(this->length == mask.size());
#endif

    std::vector<T> result;

    for (std::size_t j = 0; j < this->length; ++j)
      if (mask[j])
        result.emplace_back(this->elements[j]);

    return {result.size(), result};
  }

  /**
   * @brief Sets a sub-Vector given a std::vector of indices and a Vector of
   * values.
   *
   * @param indices Indices.
   * @param values Values.
   * @return Vector
   */
  void operator()(const std::vector<std::size_t> &indices,
                  const Vector &values) {
#ifndef NDEBUG // Integrity check.
    assert(indices.size() == values.length);
    for (const auto &index : indices)
      assert(index < this->length);
#endif

    for (std::size_t j = 0; j < indices.size(); ++j)
      this->elements[indices[j]] = values[j];
  }

  // COMPARISONS.

  /**
   * @brief Vector < Scalar.
   *
   * scalar Scalar.
   */
  Mask operator<(const T &scalar) const {
    Mask mask(this->length, false);

#ifdef PARALLEL
    std::transform(POLICY, this->elements.begin(), this->elements.end(),
                   mask.begin(),
                   [scalar](const auto &element) { return element < scalar; });
#else
    std::transform(this->elements.begin(), this->elements.end(), mask.begin(),
                   [scalar](const auto &element) { return element < scalar; });
#endif

    return mask;
  }

  /**
   * @brief Vector > Scalar.
   *
   * @param scalar Scalar.
   */
  Mask operator>(const T &scalar) const {
    Mask mask(this->length, false);

#ifdef PARALLEL
    std::transform(POLICY, this->elements.begin(), this->elements.end(),
                   mask.begin(),
                   [scalar](const auto &element) { return element > scalar; });
#else
    std::transform(this->elements.begin(), this->elements.end(), mask.begin(),
                   [scalar](const auto &element) { return element > scalar; });
#endif

    return mask;
  }

  /**
   * @brief Vector < Vector.
   *
   * @param vector Vector.
   * @return Mask
   */
  Mask operator<(const Vector &vector) const {
#ifndef NDEBUG // Integrity check.
    assert(this->length == vector.length);
#endif

    Mask mask(this->length, false);

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        vector.elements.begin(), mask.begin(),
        [](const auto &first, const auto &second) { return first < second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), vector.elements.begin(),
        mask.begin(),
        [](const auto &first, const auto &second) { return first < second; });
#endif

    return mask;
  }

  /**
   * @brief Vector > Vector.
   *
   * @param vector Vector.
   * @return Mask
   */
  Mask operator>(const Vector &vector) const {
#ifndef NDEBUG // Integrity check.
    assert(this->length == vector.length);
#endif

    Mask mask(this->length, false);

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        vector.elements.begin(), mask.begin(),
        [](const auto &first, const auto &second) { return first > second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), vector.elements.begin(),
        mask.begin(),
        [](const auto &first, const auto &second) { return first > second; });
#endif

    return mask;
  }

  // OPERATIONS.

  /**
   * @brief Vector unary +.
   *
   * @return Vector
   */
  Vector operator+() const { return *this; }

  /**
   * @brief Vector unary -.
   *
   * @return Vector
   */
  Vector operator-() const {
    Vector result{this->length};

#ifdef PARALLEL
    std::transform(POLICY, this->elements.begin(), this->elements.end(),
                   result.elements.begin(),
                   [](const auto &element) { return -element; });
#else
    std::transform(this->elements.begin(), this->elements.end(),
                   result.elements.begin(),
                   [](const auto &element) { return -element; });
#endif

    return result;
  }

  /**
   * @brief Scalar product.
   *
   * @param scalar Scalar.
   * @return Vector
   */
  Vector operator*(const T &scalar) const {
    Vector result{this->length};

#ifdef PARALLEL
    std::transform(POLICY, this->elements.begin(), this->elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element * scalar; });
#else
    std::transform(this->elements.begin(), this->elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element * scalar; });
#endif

    return result;
  }

  /**
   * @brief Friend scalar product.
   *
   * @param scalar Scalar.
   * @param vector Vector.
   * @return Vector
   */
  friend Vector operator*(const T &scalar, const Vector &vector) {
    Vector result{vector.length};

#ifdef PARALLEL
    std::transform(POLICY, vector.elements.begin(), vector.elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element * scalar; });
#else
    std::transform(vector.elements.begin(), vector.elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element * scalar; });
#endif

    return result;
  }

  /**
   * @brief Scalar product and assignation.
   *
   * @param scalar Scalar.
   * @return Vector&
   */
  Vector &operator*=(const T &scalar) {
#ifdef PARALLEL
    std::transform(POLICY, this->elements.begin(), this->elements.end(),
                   this->elements.begin(),
                   [scalar](const auto &element) { return element * scalar; });
#else
    std::transform(this->elements.begin(), this->elements.end(),
                   this->elements.begin(),
                   [scalar](const auto &element) { return element * scalar; });
#endif

    return *this;
  }

  /**
   * @brief Scalar division.
   *
   * @param scalar Scalar.
   * @return Vector
   */
  Vector operator/(const T &scalar) const {
    Vector result{this->length};

#ifdef PARALLEL
    std::transform(POLICY, this->elements.begin(), this->elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element / scalar; });
#else
    std::transform(this->elements.begin(), this->elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element / scalar; });
#endif

    return result;
  }

  /**
   * @brief Friend scalar division.
   *
   * @param scalar Scalar.
   * @param vector Vector.
   * @return Vector
   */
  friend Vector operator/(const T &scalar, const Vector &vector) {
    Vector result{vector.length};

#ifdef PARALLEL
    std::transform(POLICY, vector.elements.begin(), vector.elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return scalar / element; });
#else
    std::transform(vector.elements.begin(), vector.elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return scalar / element; });
#endif

    return result;
  }

  /**
   * @brief Scalar division and assignation.
   *
   * @param scalar Scalar.
   * @return Vector&
   */
  Vector &operator/=(const T &scalar) {
#ifdef PARALLEL
    std::transform(POLICY, this->elements.begin(), this->elements.end(),
                   this->elements.begin(),
                   [scalar](const auto &element) { return element / scalar; });
#else
    std::transform(this->elements.begin(), this->elements.end(),
                   this->elements.begin(),
                   [scalar](const auto &element) { return element / scalar; });
#endif

    return *this;
  }

  /**
   * @brief Vector sum.
   *
   * @param vector Vector.
   * @return Vector
   */
  Vector operator+(const Vector &vector) const {
#ifndef NDEBUG // Integrity check.
    assert(this->length == vector.length);
#endif

    Vector result{this->length};

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        vector.elements.begin(), result.elements.begin(),
        [](const auto &first, const auto &second) { return first + second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), vector.elements.begin(),
        result.elements.begin(),
        [](const auto &first, const auto &second) { return first + second; });
#endif

    return result;
  }

  /**
   * @brief Vector sum and assignation.
   *
   * @param vector Vector.
   * @return Vector&
   */
  Vector &operator+=(const Vector &vector) {
#ifndef NDEBUG // Integrity check.
    assert(this->length == vector.length);
#endif

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        vector.elements.begin(), this->elements.begin(),
        [](const auto &first, const auto &second) { return first + second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), vector.elements.begin(),
        this->elements.begin(),
        [](const auto &first, const auto &second) { return first + second; });
#endif

    return *this;
  }

  /**
   * @brief Scalar "sum".
   *
   * @param scalar Scalar.
   * @return Vector
   */
  Vector operator+(const T &scalar) const {
    Vector result{this->length};

#ifdef PARALLEL
    std::transform(POLICY, this->elements.begin(), this->elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element + scalar; });
#else
    std::transform(this->elements.begin(), this->elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element + scalar; });
#endif

    return result;
  }

  /**
   * @brief Friend scalar "sum".
   *
   * @param scalar Scalar.
   * @param vector Vector.
   * @return Vector
   */
  friend Vector operator+(const T &scalar, const Vector &vector) {
    Vector result{vector.length};

#ifdef PARALLEL
    std::transform(POLICY, vector.elements.begin(), vector.elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element + scalar; });
#else
    std::transform(vector.elements.begin(), vector.elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element + scalar; });
#endif

    return result;
  }

  /**
   * @brief Scalar "sum" and assignation.
   *
   * @param scalar Scalar.
   * @return Vector&
   */
  Vector &operator+=(const T &scalar) {
#ifdef PARALLEL
    std::transform(POLICY, this->elements.begin(), this->elements.end(),
                   this->elements.begin(),
                   [scalar](const auto &element) { return element + scalar; });
#else
    std::transform(this->elements.begin(), this->elements.end(),
                   this->elements.begin(),
                   [scalar](const auto &element) { return element + scalar; });
#endif

    return *this;
  }

  /**
   * @brief Vector difference.
   *
   * @param vector Vector.
   * @return Vector
   */
  Vector operator-(const Vector &vector) const {
#ifndef NDEBUG // Integrity check.
    assert(this->length == vector.length);
#endif

    Vector result{this->length};

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        vector.elements.begin(), result.elements.begin(),
        [](const auto &first, const auto &second) { return first - second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), vector.elements.begin(),
        result.elements.begin(),
        [](const auto &first, const auto &second) { return first - second; });
#endif

    return result;
  }

  /**
   * @brief Vector difference and assignation.
   *
   * @param vector Vector.
   * @return Vector&
   */
  Vector &operator-=(const Vector &vector) {
#ifndef NDEBUG // Integrity check.
    assert(this->length == vector.length);
#endif

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        vector.elements.begin(), this->elements.begin(),
        [](const auto &first, const auto &second) { return first - second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), vector.elements.begin(),
        this->elements.begin(),
        [](const auto &first, const auto &second) { return first - second; });
#endif

    return *this;
  }

  /**
   * @brief Scalar "difference".
   *
   * @param scalar Scalar.
   * @return Vector
   */
  Vector operator-(const T &scalar) const {
    Vector result{this->length};

#ifdef PARALLEL
    std::transform(POLICY, this->elements.begin(), this->elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element - scalar; });
#else
    std::transform(this->elements.begin(), this->elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element - scalar; });
#endif

    return result;
  }

  /**
   * @brief Friend scalar "difference".
   *
   * @param scalar Scalar.
   * @param vector Vector.
   * @return Vector
   */
  friend Vector operator-(const T &scalar, const Vector &vector) {
    Vector result{vector.length};

#ifdef PARALLEL
    std::transform(POLICY, vector.elements.begin(), vector.elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return scalar - element; });
#else
    std::transform(vector.elements.begin(), vector.elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return scalar - element; });
#endif

    return result;
  }

  /**
   * @brief Scalar "difference" and assignation.
   *
   * @param scalar Scalar.
   * @return Vector&
   */
  Vector &operator-=(const T &scalar) {
#ifdef PARALLEL
    std::transform(POLICY, this->elements.begin(), this->elements.end(),
                   this->elements.begin(),
                   [scalar](const auto &element) { return element - scalar; });
#else
    std::transform(this->elements.begin(), this->elements.end(),
                   this->elements.begin(),
                   [scalar](const auto &element) { return element - scalar; });
#endif

    return *this;
  }

  /**
   * @brief Vector element-wise product.
   *
   * @param vector Vector.
   * @return T
   */
  Vector operator*(const Vector &vector) const {
#ifndef NDEBUG // Integrity check.
    assert(this->length == vector.length);
#endif

    Vector result{this->length};

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        vector.elements.begin(), result.elements.begin(),
        [](const auto &first, const auto &second) { return first * second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), vector.elements.begin(),
        result.elements.begin(),
        [](const auto &first, const auto &second) { return first * second; });
#endif

    return result;
  }

  /**
   * @brief Vector element-wise product and assignation.
   *
   * @param vector Vector.
   * @return T
   */
  Vector &operator*=(const Vector &vector) {
#ifndef NDEBUG // Integrity check.
    assert(this->length == vector.length);
#endif

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        vector.elements.begin(), this->elements.begin(),
        [](const auto &first, const auto &second) { return first * second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), vector.elements.begin(),
        this->elements.begin(),
        [](const auto &first, const auto &second) { return first * second; });
#endif

    return *this;
  }

  /**
   * @brief Vector element-wise division.
   *
   * @param vector Vector.
   * @return T
   */
  Vector operator/(const Vector &vector) const {
#ifndef NDEBUG // Integrity check.
    assert(this->length == vector.length);
#endif

    Vector result{this->length};

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        vector.elements.begin(), result.elements.begin(),
        [](const auto &first, const auto &second) { return first / second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), vector.elements.begin(),
        result.elements.begin(),
        [](const auto &first, const auto &second) { return first / second; });
#endif

    return result;
  }

  /**
   * @brief Vector element-wise division and assignation.
   *
   * @param vector Vector.
   * @return T
   */
  Vector &operator/=(const Vector &vector) {
#ifndef NDEBUG // Integrity check.
    assert(this->length == vector.length);
#endif

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        vector.elements.begin(), this->elements.begin(),
        [](const auto &first, const auto &second) { return first / second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), vector.elements.begin(),
        this->elements.begin(),
        [](const auto &first, const auto &second) { return first / second; });
#endif

    return *this;
  }

  // SIZE operator

  std::size_t size() const { return this->length; }

  // OUTPUT.

  /**
   * @brief Vector output.
   *
   * @param ost
   * @param vector Vector.
   * @return std::ostream&
   */
  friend std::ostream &operator<<(std::ostream &ost, const Vector &vector) {
    for (const auto &element : vector.elements)
      ost << element << " ";

    return ost;
  }
};

} // namespace pacs

#endif