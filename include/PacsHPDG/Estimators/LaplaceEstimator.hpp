/**
 * @file EstimatorLaplace.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Laplace equation error estimator computation class object.
 * @date 2025-01-11
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_ESTIMATORS_LAPLACE_ESTIMATOR_HPP
#define INCLUDE_PACSHPDG_ESTIMATORS_LAPLACE_ESTIMATOR_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"
#include "../Solvers.hpp"

#include <array>
#include <iostream>

namespace pacs
{
class LaplaceEstimator {
protected:
  // Estimates.
  Real m_estimate;
  Vector<Real> m_estimates;

  // Fits.
  Vector<Real> m_fits;

public:
  // CONTRUCTOR.
  LaplaceEstimator(const Mesh &mesh_)
      : m_estimates{mesh_.elements.size()}, m_fits{mesh_.elements.size()} {
    this->m_estimate = 0.0;
  };

  // GETTERS.
  Real estimate() const { return this->m_estimate; };
  Vector<Real> estimates() const { return this->m_estimates; };
  Vector<Real> fits() const { return this->m_fits; };

  // METHODS.
  // Polynomial fit.
  Vector<Real> polyfit(const Vector<Real> &, const Vector<Real> &,
                       const std::size_t &);

  // Compute error estimate.
  void computeEstimate(const DataLap &, const Mesh &,
                       const Laplace &);

  // Adaptive refinement.
  void mesh_refine(Mesh &, const Real &refine = 0.75,
                   const Real &speed = 1.0);

  // Print to stream file for plotting.
  void print(std::ostream &ost, const Mesh &mesh_) const {
    ost << "Dofs: " << mesh_.dofs() << std::endl;
    ost << "Estimate: " << this->estimate() << std::endl;
  };
};
} // namespace pacs

#endif