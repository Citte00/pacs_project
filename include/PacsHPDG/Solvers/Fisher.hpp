/**
 * @file Fisher.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief
 * @date 2024-11-28
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef INCLUDE_PACSHPDG_FISHER_FISHER_HPP
#define INCLUDE_PACSHPDG_FISHER_FISHER_HPP

#include "./Heat.hpp"

namespace pacs {

class Fisher : public Heat {
protected:
  // Non-linear matrix.
  Sparse<Real> m_nl_mass;
  Vector<Real> m_ch_old;

public:
  // CONSTRUCTOR.
  Fisher(const DataFKPP &data_, const Mesh &mesh_)
      : Heat(mesh_), m_nl_mass{mesh_.dofs(), mesh_.dofs()},
        m_ch_old{mesh_.dofs()} {
    this->m_t = -data_.dt;

    // Reserve space for nl_mass matrix.
    this->m_nl_mass.inner.reserve(DOFS_MAX);
    this->m_nl_mass.outer.reserve(DOFS_MAX);
    this->m_nl_mass.values.reserve(DOFS_MAX);

    // Reserve space for ch_old vector.
    this->m_ch_old.elements.reserve(DOFS_MAX);
  };

  // GETTER.
  Sparse<Real> M_alpha() const { return this->m_nl_mass; };
  Vector<Real> ch_old() const { return this->m_ch_old; };
  Vector<Real> &ch_old() { return this->m_ch_old; };

  // METHODS.
  // Assembly the fisher equation system matrices.
  void assembly(const DataFKPP &, const Mesh &);
  // Assembly the non-linear term.
  Sparse<Real> assembly_nl(const DataFKPP &, const Mesh &,
                           const Vector<Real> &);
  // Assembly the forcing term.
  void assembly_force(const DataFKPP &, const Mesh &);
  // Solver of the Heat equation.
  Vector<Real> solver(const DataFKPP &, const Mesh &, const Vector<Real> &,
                      const Vector<Real> &, const Real &TOL = 1E-15);
  // Get source function modal coefficient.
  Vector<Real> modal_source(const DataFKPP &, const Mesh &) const;
};
} // namespace pacs

#endif