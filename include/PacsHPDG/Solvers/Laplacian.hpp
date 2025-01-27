/**
 * @file Laplacian.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Laplace equation class.
 * @date 2025-01-14
 *
 * @copyright Copyright (c) 2025
 *
 */

#ifndef INCLUDE_PACSHPDG_SOLVERS_LAPLACIAN_HPP
#define INCLUDE_PACSHPDG_SOLVERS_LAPLACIAN_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"

namespace pacs {

class Laplace {
protected:
  Sparse<Real> m_mass, m_stiff, m_dg_stiff;

public:
  // CONSTRUCTOR.
  Laplace(const Mesh &mesh_)
      : m_mass{mesh_.dofs(), mesh_.dofs()}, m_stiff{mesh_.dofs(), mesh_.dofs()},
        m_dg_stiff{mesh_.dofs(), mesh_.dofs()} {};

  // GETTERS.
  Sparse<Real> M() const { return this->m_mass; };
  Sparse<Real> A() const { return this->m_stiff; };
  Sparse<Real> DG() const { return this->m_dg_stiff; };

  // METHODS.
  // Blocks.
  std::vector<std::array<std::vector<std::size_t>, 2>>
  block_mass(const Mesh &) const;
  // Assembly the laplace system matrix.
  void assembly(const DataLaplace &, const Mesh &);
  // Assembly the forcing term.
  Vector<Real> assembly_force(const DataLaplace &, const Mesh &) const;
  // Solver of the Laplace equation.
  Vector<Real> solver(const Mesh &, const Vector<Real> &,
                      const Real &TOL = 1E-15) const;
  // Get functions modal coefficients.
  Vector<Real> modal(const Mesh &, const BiFunctor &) const;
};
} // namespace pacs

#endif