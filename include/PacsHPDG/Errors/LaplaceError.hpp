/**
 * @file LaplaceErrors.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Laplace equation error computation class object.
 * @date 2025-01-11
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_ERRORS_LAPLACE_ERROR_HPP
#define INCLUDE_PACSHPDG_ERRORS_LAPLACE_ERROR_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"
#include "../Solvers.hpp"
  
#include <array>
#include <iostream>

namespace pacs {

class LaplaceError {
protected:
  std::size_t m_degree; // p.
  Real m_size;          // h.

  Real m_DG_error;
  Real m_L2_error;

  Vector<Real> m_L2_errors;
  Vector<Real> m_H1_errors;

public:
  // CONSTRUCTORS.
  LaplaceError(const Mesh &mesh_)
      : m_L2_errors{mesh_.elements.size()}, m_H1_errors{mesh_.elements.size()} {

    m_degree = 0;
    for (const auto &element : mesh_.elements)
      m_degree = (element.degree > m_degree) ? element.degree : m_degree;

    m_size = 0.0;
    for (const auto &element : mesh_.elements)
      for (const auto &p : element.element.points)
        for (const auto &q : element.element.points)
          m_size = (distance(p, q) > m_size) ? distance(p, q) : m_size;

    m_DG_error = 0.0;
    m_L2_error = 0.0;
  };

  // GETTERS.
  std::size_t p() const { return m_degree; };
  Real h() const { return m_size; };
  Real DGerror() const { return m_DG_error; };
  Real L2error() const { return m_L2_error; };
  Vector<Real> L2errors() const { return m_L2_errors; };
  Vector<Real> H1errors() const { return m_H1_errors; };

  // METHODS.
  // Compute L2 and DG errors.
  template <typename Functor>
  void computeError(const Mesh &, const Laplace &, const Functor &);

  // Compute errors.
  void computeErrors(const DataLap &, const Mesh &, const Laplace &);

  // Print to stream file for plotting.
  virtual void print(std::ostream &os, const Mesh &mesh_) const {
    os << "Elements: " << mesh_.elements.size() << "\n";
    os << "Dofs: " << mesh_.dofs() << "\n";
    os << "Degree (p): " << this->p() << "\n";
    os << "Size (h): " << this->h() << "\n";
    os << "L2 Error: " << this->L2error() << "\n";
    os << "DG Error: " << this->DGerror() << std::endl;
  };
};

} // namespace pacs

#endif