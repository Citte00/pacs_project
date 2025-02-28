/**
 * @file Matrix.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Matrix class and methods.
 * @date 2025-01-18
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_ALGEBRA_MATRIX_HPP
#define INCLUDE_PACSHPDG_ALGEBRA_MATRIX_HPP

// Type.
#include "../Base.hpp"

#include "./Vector.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace pacs {

/**
 * @brief Matrix structure.
 *
 * @tparam T
 */
template <NumericType T> struct Matrix {

  // Shape.
  std::size_t m_rows;
  std::size_t m_columns;

  // Elements.
  std::vector<T> elements;

  // CONSTRUCTORS.

  /**
   * @brief Constructs a new empty Matrix.
   *
   * @param m_rows Rows.
   * @param m_columns Columns.
   */
  inline Matrix(const std::size_t &m_rows, const std::size_t &m_columns)
      : m_rows{m_rows}, m_columns{m_columns},
        elements(m_rows * m_columns, static_cast<T>(0)) {
#ifndef NDEBUG // Integrity check.
    assert((m_rows > 0) && (m_columns > 0));
#endif
  }

  /**
   * @brief Constructs a new Matrix from a given std::vector.
   *
   * @param m_rows Rows.
   * @param m_columns Columns.
   * @param elements
   */
  inline Matrix(const std::size_t &m_rows, const std::size_t &m_columns,
                const std::vector<T> &elements)
      : m_rows{m_rows}, m_columns{m_columns},
        elements(elements.begin(), elements.end()) {
#ifndef NDEBUG // Integrity check.
    assert((m_rows > 0) && (m_columns > 0));
    assert(elements.size() == m_rows * m_columns);
#endif
  }

  /**
   * @brief Copy constructor.
   *
   * @param matrix Matrix.
   */
  inline Matrix(const Matrix &matrix)
      : m_rows{matrix.m_rows}, m_columns{matrix.m_columns},
        elements(matrix.elements.begin(), matrix.elements.end()) {}

  /**
   * @brief Copy operator.
   *
   * @param matrix Matrix.
   * @return Matrix&
   */
  inline Matrix &operator=(const Matrix &matrix) {
#ifndef NDEBUG
    assert((this->m_rows == matrix.m_rows) &&
           (this->m_columns == matrix.m_columns));
#endif

#ifdef PARALLEL
    std::copy(POLICY, matrix.elements.begin(), matrix.elements.end(),
              this->elements.begin());
#else
    std::copy(matrix.elements.begin(), matrix.elements.end(),
              this->elements.begin());
#endif

    return *this;
  }

  // READ AND WRITE.

  /**
   * @brief Const call operator, returns the (i, j)-th element.
   *
   * @param j Index.
   * @param k Index.
   * @return T
   */
  inline T operator()(const std::size_t &j, const std::size_t &k) const {
#ifndef NDEBUG // Integrity check.
    assert((j < this->m_rows) && (k < this->m_columns));
#endif

    return this->elements[j * this->m_columns + k];
  }

  /**
   * @brief Call operator, returns a reference to the (i, j)-th element.
   *
   * @param j Index.
   * @param k Index.
   * @return T&
   */
  inline T &operator()(const std::size_t &j, const std::size_t &k) {
#ifndef NDEBUG // Integrity check.
    assert((j < this->m_rows) && (k < this->m_columns));
#endif

    return this->elements[j * this->m_columns + k];
  }

  /**
   * @brief Returns the j-th row as a Vector.
   *
   * @param j Index.
   * @return Vector<T>
   */
  Vector<T> row(const std::size_t j) const {
#ifndef NDEBUG // Integrity check.
    assert(j < this->m_rows);
#endif

    Vector<T> row{this->m_columns};

#ifdef PARALLEL
    std::copy(POLICY, this->elements.begin() + j * this->m_columns,
              this->elements.begin() + (j + 1) * this->m_columns,
              row.elements.begin());
#else
    std::copy(this->elements.begin() + j * this->m_columns,
              this->elements.begin() + (j + 1) * this->m_columns,
              row.elements.begin());
#endif

    return row;
  }

  /**
   * @brief Sets the j-th row to the given scalar.
   *
   * @param j Index.
   * @param scalar Scalar.
   */
  void row(const std::size_t j, const T &scalar) {
#ifndef NDEBUG // Integrity check.
    assert(j < this->m_rows);
#endif

#ifdef PARALLEL
    std::for_each_n(POLICY, this->elements.begin() + j * this->m_columns,
                    this->m_columns,
                    [scalar](auto &element) { element = scalar; });
#else
    std::for_each_n(this->elements.begin() + j * this->m_columns, this->m_columns,
                    [scalar](auto &element) { element = scalar; });
#endif
  }

  /**
   * @brief Sets the j-th row to the given Vector.
   *
   * @param j Index.
   * @param vector Vector.
   */
  void row(const std::size_t j, const Vector<T> &vector) {
#ifndef NDEBUG // Integrity check.
    assert(j < this->m_rows);
    assert(vector.length == this->m_columns);
#endif

#ifdef PARALLEL
    std::copy(POLICY, vector.elements.begin(), vector.elements.end(),
              this->elements.begin() + j * this->m_columns);
#else
    std::copy(vector.elements.begin(), vector.elements.end(),
              this->elements.begin() + j * this->m_columns);
#endif
  }

  /**
   * @brief Returns the k-th column as a Vector.
   *
   * @param k Index.
   * @return Vector<T>
   */
  Vector<T> column(const std::size_t &k) const {
#ifndef NDEBUG // Integrity check.
    assert(k < this->m_columns);
#endif

    Vector<T> column{this->m_rows};

    for (std::size_t j = 0; j < this->m_rows; ++j)
      column[j] = this->elements[j * this->m_columns + k];

    return column;
  }

  /**
   * @brief Sets the k-th column to the given scalar.
   *
   * @param k Index.
   * @param scalar Scalar.
   */
  void column(const std::size_t &k, const T &scalar) {
#ifndef NDEBUG // Integrity check.
    assert(k < this->m_columns);
#endif

    for (std::size_t j = 0; j < this->m_rows; ++j)
      this->elements[j * this->m_columns + k] = scalar;
  }

  /**
   * @brief Sets the k-th column to the given vector.
   *
   * @param k Index.
   * @param vector Vector.
   */
  void column(const std::size_t &k, const Vector<T> &vector) {
#ifndef NDEBUG // Integrity check.
    assert(k < this->m_columns);
    assert(vector.length == this->m_rows);
#endif

    for (std::size_t j = 0; j < this->m_rows; ++j)
      this->elements[j * this->m_columns + k] = vector.elements[j];
  }

  // SIZE.

  /**
   * @brief Returns the Matrix size.
   *
   * @return std::size_t
   */
  inline std::size_t size() const { return this->m_rows * this->m_columns; }

  // SHAPE.

  /**
   * @brief Returns the reshaped Matrix.
   *
   * @param m_rows Rows.
   * @param m_columns Columns.
   */
  inline void reshape(const std::size_t &m_rows, const std::size_t &m_columns) {
    this->m_rows = m_rows;
    this->m_columns = m_columns;
    this->elements.resize(m_rows * m_columns, 0.0);
  };

  /**
   * @brief Returns the transpose matrix.
   *
   * @return Matrix
   */
  Matrix transpose() const {
    Matrix transpose{this->m_columns, this->m_rows};

    for (std::size_t j = 0; j < this->m_rows; ++j)
      for (std::size_t k = 0; k < this->m_columns; ++k)
        transpose.elements[k * this->m_rows + j] =
            this->elements[j * this->m_columns + k];

    return transpose;
  }

  /**
   * @brief Returns the diagonal.
   *
   * @return Matrix
   */
  Matrix diagonal() const {
#ifndef NDEBUG // Integrity check.
    assert(this->m_rows == this->m_columns);
#endif

    Matrix diagonal{this->m_rows, this->m_columns};

    for (std::size_t j = 0; j < this->m_rows; ++j)
      diagonal.elements[j * (this->m_columns + 1)] =
          this->elements[j * (this->m_columns + 1)];

    return diagonal;
  }

  /**
   * @brief Returns the lower triangular part of the Matrix.
   *
   * @return Matrix
   */
  Matrix lower() const {
#ifndef NDEBUG // Integrity check.
    assert(this->m_rows == this->m_columns);
#endif

    Matrix lower{this->m_rows, this->m_columns};

    for (std::size_t j = 0; j < this->m_rows; ++j)
#ifdef PARALLEL
      std::copy(POLICY, this->elements.begin() + j * this->m_columns,
                this->elements.begin() + j * this->m_columns + j,
                lower.elements.begin() + j * this->m_columns);
#else
      std::copy(this->elements.begin() + j * this->m_columns,
                this->elements.begin() + j * this->m_columns + j,
                lower.elements.begin() + j * this->m_columns);
#endif

    return lower;
  }

  /**
   * @brief Returns the upper triangular part of the Matrix.
   *
   * @return Matrix
   */
  Matrix upper() const {
#ifndef NDEBUG // Integrity check.
    assert(this->m_rows == this->m_columns);
#endif

    Matrix upper{this->m_rows, this->m_columns};

    for (std::size_t j = 0; j < this->m_rows; ++j)
#ifdef PARALLEL
      std::copy(POLICY, this->elements.begin() + j * this->m_columns + j + 1,
                this->elements.begin() + (j + 1) * this->m_columns,
                upper.elements.begin() + j * this->m_columns + j + 1);
#else
      std::copy(this->elements.begin() + j * this->m_columns + j + 1,
                this->elements.begin() + (j + 1) * this->m_columns,
                upper.elements.begin() + j * this->m_columns + j + 1);
#endif

    return upper;
  }

  // OPERATIONS.

  /**
   * @brief Matrix unary +.
   *
   * @return Matrix
   */
  Matrix operator+() const { return *this; }

  /**
   * @brief Matrix unary -.
   *
   * @return Matrix
   */
  Matrix operator-() const {
    Matrix result{this->m_rows, this->m_columns};

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
   * @brief Matrix scalar product.
   *
   * @param scalar Scalar.
   * @return Matrix
   */
  Matrix operator*(const T &scalar) const {
    Matrix result{this->m_rows, this->m_columns};

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
   * @brief Friend Matrix scalar product.
   *
   * @param scalar Scalar.
   * @param matrix Matrix.
   * @return Matrix
   */
  friend Matrix operator*(const T &scalar, const Matrix &matrix) {
    Matrix result{matrix.m_rows, matrix.m_columns};

#ifdef PARALLEL
    std::transform(POLICY, matrix.elements.begin(), matrix.elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element * scalar; });
#else
    std::transform(matrix.elements.begin(), matrix.elements.end(),
                   result.elements.begin(),
                   [scalar](const auto &element) { return element * scalar; });
#endif

    return result;
  }

  /**
   * @brief Matrix scalar product and assignation.
   *
   * @param scalar Scalar.
   * @return Matrix&
   */
  Matrix &operator*=(const T scalar) {
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
   * @brief Matrix scalar division.
   *
   * @param scalar Scalar.
   * @return Matrix
   */
  Matrix operator/(const T &scalar) const {
    Matrix result{this->m_rows, this->m_columns};

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
   * @brief Matrix scalar division and assignation.
   *
   * @param scalar Scalar.
   * @return Matrix&
   */
  Matrix &operator/=(const T scalar) {
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
   * @brief Matrix sum.
   *
   * @param matrix Matrix.
   * @return Matrix
   */
  Matrix operator+(const Matrix &matrix) const {
#ifndef NDEBUG // Integrity check.
    assert((this->m_rows == matrix.m_rows) &&
           (this->m_columns == matrix.m_columns));
#endif

    Matrix result{this->m_rows, this->m_columns};

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        matrix.elements.begin(), result.elements.begin(),
        [](const auto &first, const auto &second) { return first + second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), matrix.elements.begin(),
        result.elements.begin(),
        [](const auto &first, const auto &second) { return first + second; });
#endif

    return result;
  }

  /**
   * @brief Matrix sum and assignation.
   *
   * @param matrix Matrix.
   * @return Matrix&
   */
  Matrix &operator+=(const Matrix &matrix) {
#ifndef NDEBUG // Integrity check.
    assert((this->m_rows == matrix.m_rows) &&
           (this->m_columns == matrix.m_columns));
#endif

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        matrix.elements.begin(), this->elements.begin(),
        [](const auto &first, const auto &second) { return first + second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), matrix.elements.begin(),
        this->elements.begin(),
        [](const auto &first, const auto &second) { return first + second; });
#endif

    return *this;
  }

  /**
   * @brief Matrix difference.
   *
   * @param matrix Matrix.
   * @return Matrix
   */
  Matrix operator-(const Matrix &matrix) const {
#ifndef NDEBUG // Integrity check.
    assert((this->m_rows == matrix.m_rows) &&
           (this->m_columns == matrix.m_columns));
#endif

    Matrix result{this->m_rows, this->m_columns};

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        matrix.elements.begin(), result.elements.begin(),
        [](const auto &first, const auto &second) { return first - second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), matrix.elements.begin(),
        result.elements.begin(),
        [](const auto &first, const auto &second) { return first - second; });
#endif

    return result;
  }

  /**
   * @brief Matrix difference and assignation.
   *
   * @param matrix Matrix.
   * @return Matrix&
   */
  Matrix &operator-=(const Matrix &matrix) {
#ifndef NDEBUG // Integrity check.
    assert((this->m_rows == matrix.m_rows) &&
           (this->m_columns == matrix.m_columns));
#endif

#ifdef PARALLEL
    std::transform(
        POLICY, this->elements.begin(), this->elements.end(),
        matrix.elements.begin(), this->elements.begin(),
        [](const auto &first, const auto &second) { return first - second; });
#else
    std::transform(
        this->elements.begin(), this->elements.end(), matrix.elements.begin(),
        this->elements.begin(),
        [](const auto &first, const auto &second) { return first - second; });
#endif

    return *this;
  }

  /**
   * @brief Matrix * Vector product.
   *
   * @param vector Vector.
   * @return Vector<T>
   */
  Vector<T> operator*(const Vector<T> &vector) const {
#ifndef NDEBUG // Integrity check.
    assert(this->m_columns == vector.length);
#endif

    Vector<T> result{this->m_rows};

    for (std::size_t j = 0; j < this->m_rows; ++j)
#ifdef PARALLEL
      result.elements[j] = std::transform_reduce(
          POLICY, vector.elements.begin(), vector.elements.end(),
          this->elements.begin() + j * this->m_columns, static_cast<T>(0),
          std::plus{},
          [](const auto &first, const auto &second) { return first * second; });
#else
      result.elements[j] = std::inner_product(
          vector.elements.begin(), vector.elements.end(),
          this->elements.begin() + j * this->m_columns, static_cast<T>(0));
#endif

    return result;
  }

  /**
   * @brief Friend Vector * Matrix product.
   *
   * @param vector Vector.
   * @param matrix Matrix.
   * @return Vector<T>
   */
  friend Vector<T> operator*(const Vector<T> &vector, const Matrix &matrix) {
#ifndef NDEBUG // Integrity check.
    assert(vector.length == matrix.m_rows);
#endif

    Vector<T> result{matrix.m_columns};

    for (std::size_t j = 0; j < matrix.m_columns; ++j)
      for (std::size_t k = 0; k < matrix.m_rows; ++k)
        result.elements[j] +=
            vector.elements[k] * matrix.elements[k * matrix.m_rows + j];

    return result;
  }

  /**
   * @brief Matrix * Matrix product.
   *
   * @param matrix Matrix.
   * @return Matrix
   */
  Matrix operator*(const Matrix &matrix) const {
#ifndef NDEBUG // Integrity check.
    assert(this->m_columns == matrix.m_rows);
#endif

    Matrix result{this->m_rows, matrix.m_columns};

    for (std::size_t j = 0; j < this->m_rows; ++j)
      for (std::size_t k = 0; k < matrix.m_columns; ++k)
        for (std::size_t h = 0; h < this->m_columns; ++h)
          result.elements[j * result.m_columns + k] +=
              this->elements[j * this->m_columns + h] *
              matrix.elements[h * matrix.m_columns + k];

    return result;
  }

  // OUTPUT.

  /**
   * @brief Matrix output.
   *
   * @param ost
   * @param matrix Matrix.
   * @return std::ostream&
   */
  friend std::ostream &operator<<(std::ostream &ost, const Matrix &matrix) {
    for (std::size_t j = 0; j < matrix.m_rows; ++j) {
      for (std::size_t k = 0; k < matrix.m_columns; ++k)
        ost << matrix.elements[j * matrix.m_columns + k] << " ";

      if (j < matrix.m_rows - 1)
        ost << std::endl;
    }
    return ost;
  }
};

} // namespace pacs

#endif