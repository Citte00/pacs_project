/**
 * @file Estimators.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Error estimates computation classes.
 * @date 2025-01-20
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef INCLUDE_PACSHPDG_ERRORS_ESTIMATORS_HPP
#define INCLUDE_PACSHPDG_ERRORS_ESTIMATORS_HPP

#include <array>
#include <iostream>

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"
#include "../Solvers.hpp"

namespace pacs {

/**
 * @brief Laplace equation error estimator class.
 * 
 */
class LaplaceEstimator {
protected:
  // DOFs.
  std::size_t m_dofs;
  // Estimates.
  Real m_estimate;
  Vector<Real> m_estimates, m_fits;
  // Refined mesh.
  Mesh m_mesh;

public:
  // CONSTRUCTOR.
  LaplaceEstimator(const Mesh &mesh_)
      : m_estimates{mesh_.elements.size()}, m_fits{mesh_.elements.size()},
        m_mesh{mesh_} {
    m_dofs = mesh_.dofs();
    m_estimate = 0.0;
  };

  // GETTERS.
  std::size_t dofs() const { return this->m_dofs; };
  Real estimate() const { return this->m_estimate; };
  Vector<Real> estimates() const { return this->m_estimates; };
  Vector<Real> fits() const { return this->m_fits; };
  Mesh mesh() const { return this->m_mesh; };

  // METHODS.
  // Compute error estimates.
  void computeEstimates(const DataLaplace &, const Laplace &,
                        const Vector<Real> &numerical);
  // Polynomial fit.
  Vector<Real> polyfit(const Vector<Real> &, const Vector<Real> &,
                       const std::size_t &) const;

  // Refinement.
  void mesh_refine_size(const Mask &);
  void mesh_refine_degree(const Mask &);

  // Elements to refine refinement.
  std::array<Mask, 2> find_elem_to_refine(const Real &refine = 0.75, const Real &speed = 1.0);

  // Friend operator<< for output printing.
  friend std::ostream &operator<<(std::ostream &ost,
                                  const LaplaceEstimator &estimator) {
    ost << "Dofs: " << estimator.dofs() << std::endl;
    ost << "Estimate: " << estimator.estimate() << std::endl;
    return ost;
  };

  // Outputs the error estimate distribution.
  void write(const std::string &, const bool &estimates = false);
};

/**
 * @brief Heat equation error estimator class.
 * 
 */
class HeatEstimator : public LaplaceEstimator {
public:
  // CONSTRUCTOR.
  HeatEstimator(const Mesh &mesh_) : LaplaceEstimator(mesh_) {};

  // METHODS.
  // Compute error estimates.
  void computeEstimates(const DataHeat &, const Heat &, const Vector<Real> &, const Vector<Real> &);
};

/**
 * @brief Fisher equation error estimator class.
 *
 */
class FisherEstimator : public HeatEstimator {
public:
  // CONSTRUCTOR.
  FisherEstimator(const Mesh &mesh_) : HeatEstimator(mesh_) {};

  // METHODS.
  // Compute error estimates.
  void computeEstimates(const DataFKPP &, const Fisher &,
                        const Vector<Real> &);
};

} // namespace pacs

#endif