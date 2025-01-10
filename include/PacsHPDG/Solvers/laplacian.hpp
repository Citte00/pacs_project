/**
 * @file laplacian.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Class for laplacian solver.
 * @date 2025-01-10
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
  // Matrices.
  Sparse<Real> m_mass;
  Sparse<Real> m_stiff;
  Sparse<Real> m_DG_stiff;

  // Forcing term.
  Vector<Real> m_forcing;

public:
  // CONSTRUCTOR.
  Laplace(const Mesh &mesh_)
      : m_mass{mesh_.dofs(), mesh_.dofs()}, m_stiff{mesh_.dofs(), mesh_.dofs()},
        m_DG_stiff{mesh_.dofs(), mesh_.dofs()}, m_forcing{mesh_.dofs()} {};

  // GETTER.
  std::array<Sparse<Real>, 3> matrices() const { return {m_mass, m_stiff, m_DG_stiff}; };
  Vector<Real> forcing() const { return m_forcing; };

  // METHODS.
  // Blocks.
  std::vector<std::array<std::vector<std::size_t>, 2>> block_mass(const Mesh &);

  // Assembly of laplacian matrix.
  void assembly(const DataLap &, const Mesh &);

  // Assembly of forcing term.
  void assemblyforce(const DataLap &, const Mesh &);

  // Laplacian solver.
  Vector<Real> solver(const Mesh &, const Real &TOL = 1E-15);

  // Evaluate coefficient.
  Vector<Real> evaluateCoeff(const Mesh &, const BiFunctor &);

  // Get coefficients of source function.
  Vector<Real> evaluateSource(const DataLap &, const Mesh &);
};

} // namespace pacs

#endif