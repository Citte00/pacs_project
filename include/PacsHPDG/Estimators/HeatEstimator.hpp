/**
 * @file HeatEstimator.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Heat equation error estimator computation class object.
 * @date 2025-01-11
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_ESTIMATORS_HEAT_ESTIMATOR_HPP
#define INCLUDE_PACSHPDG_ESTIMATORS_HEAT_ESTIMATOR_HPP

#include "./LaplaceEstimator.hpp"

namespace pacs {
class HeatEstimator : public LaplaceEstimator {
public:
  // CONTRUCTOR.
  HeatEstimator(const Mesh &mesh_) : LaplaceEstimator(mesh_) {};

  // METHODS.
  // Compute error estimate.
  void computeEstimate(const DataHeat &, const Mesh &, const Heat &);
};
} // namespace pacs

#endif