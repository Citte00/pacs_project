/**
 * @file Errors.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief
 * @date 2024-12-22
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef INCLUDE_PACSHPDG_ERRORS_ERRORS_HPP
#define INCLUDE_PACSHPDG_ERRORS_ERRORS_HPP

#include <array>
#include <iostream>

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"

namespace pacs {

/**
 * @brief Polymorphism base class error.
 *
 */
class Error {
protected:
  std::size_t elements;
  std::size_t dofs;
  std::size_t degree; // p.
  Real size;          // h.

  Real dg_error;
  Real l2_error;

  Vector<Real> l2_errors;
  Vector<Real> h1_errors;

public:
  Error(std::size_t num_elements = 0)
      : elements{num_elements}, l2_errors{num_elements},
        h1_errors{num_elements} {};

  // Virtual constructor for polymorphic behaviour.
  virtual ~Error() = default;

  // Virtual print method for polymorphic behavior.
  virtual void print(std::ostream &os) const = 0;

  // Virtual getters.
  virtual std::size_t getDofs() const = 0;
  virtual Real getElements() const = 0;
  virtual Vector<Real> getDegree() const = 0;
  virtual Vector<Real> getSize() const = 0;
  virtual Real getL2Error() const = 0;
  virtual Real getDGError() const = 0;
  virtual Vector<Real> getl2Error() const = 0;
  virtual Vector<Real> getH1Error() const = 0;

  // Virtual setters.
  virtual void setDofs(const size_t &) = 0;
  virtual void setElements(const size_t &) = 0;
  virtual void setDegree(const size_t &) = 0;
  virtual void setSize(const Real &) = 0;
  virtual void setL2Error(const Real &) = 0;
  virtual void setDGError(const Real &) = 0;
  virtual void setl2Error(const Vector<Real> &) = 0;
  virtual void setH1Error(const Vector<Real> &) = 0;

  // Friend operator<< for polymorphic printing.
  friend std::ostream &operator<<(std::ostream &os, const Error &error) {
    error.print(os);
    return os;
  };
};

/**
 * @brief Laplace equation error class.
 *
 */
class LapError : public Error {
public:
  // CONSTRUCTORS.
  LapError(const Mesh &, const std::array<Sparse<Real>, 2> &,
           const Vector<Real> &, const Functor &, const TwoFunctor &);

  // Getters.
  std::size_t getDofs() const { return dofs; };
  Real getElements() const { return elements; };
  Vector<Real> getDegree() const { return degree; };
  Vector<Real> getSize() const { return size; };
  Real getL2Error() const { return l2_error; };
  Real getDGError() const { return dg_error; };
  Vector<Real> getl2Error() const { return l2_errors; };
  Vector<Real> getH1Error() const { return h1_errors; };

  // Setters.
  void setDofs(const size_t &dofs_) { dofs = dofs_; };
  void setElements(const size_t &elements_) { elements = elements_; };
  void setDegree(const size_t &degree_) { degree = degree_; };
  void setSize(const Real &size_) { size = size_; };
  void setL2Error(const Real &error_) { l2_error = error_; };
  void setDGError(const Real &error_) { dg_error = error_; };
  void setl2Error(const Vector<Real> &errors_) { l2_errors = errors_; };
  void setH1Error(const Vector<Real> &errors_) { h1_errors = errors_; };

  // Output
  void print(std::ostream &os) const {
    os << "Elements: " << elements << "\n";
    os << "Dofs: " << dofs << "\n";
    os << "Degree (p): " << degree << "\n";
    os << "Size (h): " << size << "\n";
    os << "L2 Error: " << l2_error << "\n";
    os << "DG Error: " << dg_error << std::endl;
  };
};

/**
 * @brief Heat equation error class.
 *
 */
class HeatError : public Error {
public:
  // CONSTRUCTORS.
  HeatError(const DataHeat &, const Mesh &, const std::array<Sparse<Real>, 2> &,
            const Vector<Real> &, const Real &);

  // Getters.
  std::size_t getDofs() const { return dofs; };
  Real getElements() const { return elements; };
  Vector<Real> getDegree() const { return degree; };
  Vector<Real> getSize() const { return size; };
  Real getL2Error() const { return l2_error; };
  Real getDGError() const { return dg_error; };
  Vector<Real> getl2Error() const { return l2_errors; };
  Vector<Real> getH1Error() const { return h1_errors; };

  // Setters.
  void setDofs(const size_t &dofs_) { dofs = dofs_; };
  void setElements(const size_t &elements_) { elements = elements_; };
  void setDegree(const size_t &degree_) { degree = degree_; };
  void setSize(const Real &size_) { size = size_; };
  void setL2Error(const Real &error_) { l2_error = error_; };
  void setDGError(const Real &error_) { dg_error = error_; };
  void setl2Error(const Vector<Real> &errors_) { l2_errors = errors_; };
  void setH1Error(const Vector<Real> &errors_) { h1_errors = errors_; };

  // Output
  void print(std::ostream &os) const {
    os << "Elements: " << elements << "\n";
    os << "Dofs: " << dofs << "\n";
    os << "Degree (p): " << degree << "\n";
    os << "Size (h): " << size << "\n";
    os << "L2 Error: " << l2_error << "\n";
    os << "DG Error: " << dg_error << std::endl;
  };
};
}

#endif