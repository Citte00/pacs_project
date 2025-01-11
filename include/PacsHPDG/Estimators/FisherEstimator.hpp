/**
 * @file FisherEstimator.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Fisher-KPP equation error estimator computation class object.
 * @date 2025-01-11
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_ESTIMATORS_FISHER_ESTIMATOR_HPP
#define INCLUDE_PACSHPDG_ESTIMATORS_FISHER_ESTIMATOR_HPP

#include "./HeatEstimator.hpp"

namespace pacs {
class FisherEstimator : public HeatEstimator {
public:
  // CONTRUCTOR.
  FisherEstimator(const Mesh &mesh_) : HeatEstimator(mesh_) {};

  // METHODS.
  // Compute error estimate.
  void computeEstimate(const DataFKPP &, const Mesh &, const Fisher &);
};
} // namespace pacs

#endif