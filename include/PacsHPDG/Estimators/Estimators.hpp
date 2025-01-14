/**
 * @file Estimators.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Laplace error estimator class.
 * @date 2024-12-21
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef INCLUDE_PACSHPDG_ESTIMATORS_LAPLACE_ESTIMATORS_HPP
#define INCLUDE_PACSHPDG_ESTIMATORS_LAPLACE_ESTIMATORS_HPP

#include <array>
#include <iostream>

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"
#include "../Laplacian.hpp"

namespace pacs {

class LaplaceEstimator {
protected:
  // DOFs.
  std::size_t m_dofs;
  // Estimates.
  Real m_estimate;
  Vector<Real> m_estimates, m_fits;

public:
  // CONSTRUCTOR.
  LaplaceEstimator(const Mesh &mesh_)
      : m_estimates{mesh_.elements.size()}, m_fits{mesh_.elements.size()} {
    m_dofs = mesh_.dofs();
    m_estimate = 0.0;
  };

  // GETTERS.
  std::size_t dofs() const { return this->m_dofs; };
  Real estimate() const { return this->m_estimate; };
  Vector<Real> estimates() const { return this->m_estimates; };
  Vector<Real> fits() const { return this->m_fits; };

  // METHODS.
  // Compute error estimates.
  void computeEstimates(const DataLaplace &, const Mesh &, const Laplace &,
                        const Vector<Real> &numerical);
  // Polynomial fit.
  Vector<Real> polyfit(const Vector<Real> &, const Vector<Real> &,
                       const std::size_t &);

  // Refinement.
  void mesh_refine_size(Mesh &, const Mask &);
  void mesh_refine_degree(Mesh &, const Mask &);
  
  // Mesh refinement.
  void mesh_refine(Mesh &, const LaplaceEstimator &, const Real &refine = 0.75,
                   const Real &speed = 1.0);

  // Friend operator<< for output printing.
  friend std::ostream &operator<<(std::ostream &ost,
                                  const LaplaceEstimator &estimator) {
    ost << "Dofs: " << estimator.dofs() << std::endl;
    ost << "Estimate: " << estimator.estimate() << std::endl;
  }
};

} // namespace pacs

#endif