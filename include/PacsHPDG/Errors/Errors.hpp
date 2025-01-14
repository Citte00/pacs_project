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
#include "../Laplacian.hpp"

namespace pacs {

/**
 * @brief Polymorphism base class error.
 *
 */
class LaplaceError {
protected:
  std::size_t m_degree, m_dofs;
  Real m_size;

  Real m_dg_error;
  Real m_l2_error;

  Vector<Real> m_l2_errors;
  Vector<Real> m_h1_errors;

public:
  LaplaceError(const Mesh &mesh_)
      : m_l2_errors{mesh_.elements.size()}, m_h1_errors{mesh_.elements.size()} {

    this->m_dofs = mesh_.dofs();

    this->m_degree = 0;
    for (const auto &element : mesh_.elements)
      this->m_degree =
          (element.degree > this->m_degree) ? element.degree : this->m_degree;

    this->m_size = 0.0;
    for (const auto &element : mesh_.elements)
      for (const auto &p : element.element.points)
        for (const auto &q : element.element.points)
          this->m_size =
              (distance(p, q) > this->m_size) ? distance(p, q) : this->m_size;
  };

  // Getter.
  std::size_t dofs() const { return this->m_dofs; };
  std::size_t p() const { return this->m_degree; };
  Real h() const { return this->m_size; };
  Real L2error() const { return this->m_l2_error; };
  Real DGerror() const { return this->m_dg_error; };
  Vector<Real> L2errors() const { return this->m_l2_errors; };
  Vector<Real> H1errors() const { return this->m_h1_errors; };

  void computeErrors(const DataLaplace &, const Mesh &, const Laplace &, const Vector<Real> &);

  // Friend operator<< for polymorphic printing.
  friend std::ostream &operator<<(std::ostream &ost, const LaplaceError &error) {
    ost << "Elements: " << error.L2errors().size() << "\n";
    ost << "Dofs: " << error.dofs() << "\n";
    ost << "Degree (p): " << error.p() << "\n";
    ost << "Size (h): " << error.h() << "\n";
    ost << "L2 Error: " << error.L2error() << "\n";
    return ost << "DG Error: " << error.DGerror() << std::endl;
  };
};

} // namespace pacs

#endif