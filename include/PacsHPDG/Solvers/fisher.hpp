/**
 * @file fisher.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Class for Fisher-KPP equation solver.
 * @date 2025-01-10
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_SOLVERS_FISHER_HPP
#define INCLUDE_PACSHPDG_SOLVERS_FISHER_HPP

#include "./heat.hpp"

namespace pacs {

class Fisher : public Heat {
private:
  Sparse<Real> m_nl_mass;
  Vector<Real> m_ch_old;

public:
    // CONSTRUCTOR.
  Fisher(const Mesh &mesh_)
      : Heat(mesh_), m_nl_mass{mesh_.dofs(), mesh_.dofs()},
        m_ch_old{mesh_.dofs()} {};

  // GETTER.
  Vector<Real> ch_old() const { return m_ch_old; };
  Sparse<Real> nl_mass() const { return m_nl_mass; };

  // METHODS.
  // Assembly of the system matrices.
  void assembly(const DataFKPP &, const Mesh &);

  // Assembly of the forcing term.
  void assemblyforce(const DataFKPP &, const Mesh &);

  // Assembly of the non-linear matrix.
  Sparse<Real> assemblyNL(const DataFKPP &, const Mesh &, const TriFunctor &);

  // Solve the heat equation system.
  void solver(const DataFKPP &, const Mesh &, const Real &TOL = 1E-15);

  // Get coefficients of source function.
  Vector<Real> evaluateSource(const DataFKPP &, const Mesh &, const Real &) const;

  // Get initial condition.
  void evaluateIC(const DataFKPP &, const Mesh &);
};

} // namespace pacs

#endif