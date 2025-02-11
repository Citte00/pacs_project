/**
 * @file Heat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief
 * @date 2024-12-16
 *
 * @copyright Copyright (c) 2024
 *
 */
#ifndef INCLUDE_PACSHPDG_SOLVERS_HEAT_HPP
#define INCLUDE_PACSHPDG_SOLVERS_HEAT_HPP

#include "./Laplacian.hpp"

namespace pacs {

class Heat : public Laplace {
protected:
  // Forcing term and numerical solution.
  Vector<Real> m_forcing;

  // Time step.
  Real m_t;

public:
  // CONSTRUCTOR.
  Heat(const Mesh &mesh_)
      : Laplace(mesh_), m_forcing{mesh_.dofs()}, m_t{0.0} {
    this->m_forcing.elements.reserve(DOFS_MAX);
  };

  // GETTERS.
  Vector<Real> forcing() const { return this->m_forcing; };
  Vector<Real> &forcing() { return this->m_forcing; };
  Real t() const { return this->m_t; };
  Real &t() { return this->m_t; };

  // METHODS.
  // Update matrices and forcing term in adaptive framework.
  void update(const DataHeat &, const Mesh &);

  // Assembly the heat equation system matrices.
  void assembly(const DataHeat &, const Mesh &);
  
  // Assembly the forcing term.
  void assembly_force(const DataHeat &, const Mesh &);
  
  // Solver of the Heat equation.
  Vector<Real> solver(const DataHeat &, const Mesh &, const Vector<Real> &,
                      const Vector<Real> &, const Real &TOL = 1E-15);
  
  // Get functions modal coefficients.
  Vector<Real> modal(const Mesh &, const TriFunctor &) const;
  
  // Get source function modal coefficient.
  Vector<Real> modal_source(const DataHeat &, const Mesh &) const;
  
  // hp-adaptive methods.
  
  // Construct matrix with base indeces for each degree.
  Matrix<int> transition(const std::size_t &) const;
  
  // prolong solution for p.
  Vector<Real> prolong_solution_p(const Mesh &, const Mesh &,
                                  const Vector<Real> &, const Mask &) const;
  
  // prolong solution for h.
  Vector<Real> prolong_solution_h(const Mesh &, const Mesh &,
                                  const Vector<Real> &, const Mask &) const;
};
} // namespace pacs

#endif