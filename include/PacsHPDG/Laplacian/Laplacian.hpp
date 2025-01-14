/**
 * @file Laplacian.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief
 * @date 2024-05-08
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef LAPLACIAN_MATRIX__PACS
#define LAPLACIAN_MATRIX__PACS

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"

namespace pacs {

class Laplace {
private:
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
  std::vector<std::array<std::vector<std::size_t>, 2>> block_mass(const Mesh &) const;
  // Assembly the laplace system matrix.
  void assembly(const DataLaplace &, const Mesh &);
  // Assembly the forcing term.
  Vector<Real> forcing(const DataLaplace &, const Mesh &);
  // Solver of the Laplace equation.
  Vector<Real> lapsolver(const Mesh &, const Vector<Real> &,
                         const Real &TOL = 1E-15);
  // Get functions modal coefficients.
  Vector<Real> modal(const Mesh &, const BiFunctor &);
};
} // namespace pacs

#endif