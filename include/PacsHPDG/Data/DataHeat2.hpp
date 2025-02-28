/**
 * @file DataHeat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Data for the Heat equation problem.
 * @date 2024-12-19
 *
 * @copyright Copyright (c) 2024
 *
 */
#ifndef INCLUDE_PACSHPDG_DATA_DATAHEAT2_HPP
#define INCLUDE_PACSHPDG_DATA_DATAHEAT2_HPP

#include <array>
#include <iostream>

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"

namespace pacs {

struct DataHeat2 : public DataHeat {

  std::vector<Point> domain = {
      {-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};

  // Forcing Term
  SourceFunc source_f = [](Vector<Real> x, Vector<Real> y, Real t,
                           Vector<Real> D) {
    Vector<Real> result{x.length};

    for (size_t i = 0; i < x.length; i++)
      result[i] =
          y[i] / 10.0 *
              std::exp((1.0L - std::exp((x[i] - t * y[i]) / 10.0)) /
                           (1.0L - std::exp(-1 / 10.0)) +
                       (x[i] - t * y[i]) / 10.0) /
              (1.0L - std::exp(-1 / 10.0)) +
          D[i] * (-1.0 / (10.0 * (1.0L - std::exp(-1 / 10.0))) *
                      std::exp((1.0L - std::exp((x[i] - t * y[i]) / 10.0)) /
                                   (1.0L - std::exp(-1 / 10.0)) +
                               (x[i] - t * y[i]) / 10.0) *
                      (1.0 / 10.0 - std::exp((x[i] - t * y[i]) / 10.0)) /
                      (10.0 * (1.0L - std::exp(-1 / 10.0))) +
                  t * t / (10.0 * (1.0L - std::exp(-1 / 10.0))) *
                      std::exp((1.0L - std::exp((x[i] - t * y[i]) / 10.0)) /
                                   (1.0L - std::exp(-1 / 10.0)) +
                               (x[i] - t * y[i]) / 10.0) *
                      (-1.0 / 10.0 + std::exp((x[i] - t * y[i]) / 10.0)) /
                      (10.0 * (1.0L - std::exp(-1 / 10.0))));

    return result;
  };

  // Boundary Conditions
  SpatialTimeFunc DirBC = [](Vector<Real> x, Vector<Real> y, Real t) {
    Vector<Real> result{x.length};

    for (size_t i = 0; i < x.length; i++)
      result[i] = std::exp((1.0L - std::exp((x[i] - t * y[i]) / 10.0)) /
                           (1.0L - std::exp(-1 / 10.0)));

    return result;
  };

  // Gradients of the Boundary Conditions
  SpatialTimeFunc DirBC_dx = [](Vector<Real> x, Vector<Real> y, Real t) {
    Vector<Real> result{x.size()};

    for (size_t i = 0; i < x.length; i++)
      result[i] = -1.0 / 10.0 *
                  std::exp((1.0L - std::exp((x[i] - t * y[i]) / 10.0)) /
                               (1.0L - std::exp(-1 / 10.0)) +
                           (x[i] - t * y[i]) / 10.0) /
                  (1.0L - std::exp(-1 / 10.0));

    return result;
  };

  SpatialTimeFunc DirBC_dy = [](Vector<Real> x, Vector<Real> y, Real t) {
    Vector<Real> result{x.size()};

    for (size_t i = 0; i < x.length; i++)
      result[i] = t / 10.0 *
                  std::exp((1.0L - std::exp((x[i] - t * y[i]) / 10.0)) /
                               (1.0L - std::exp(-1 / 10.0)) +
                           (x[i] - t * y[i]) / 10.0) /
                  (1.0L - std::exp(-1 / 10.0));

    return result;
  };

  SpatialTimeFunc DirBC_dt = [](Vector<Real> x, Vector<Real> y, Real t) {
    Vector<Real> result{x.length};

    for (size_t i = 0; i < x.length; i++)
      result[i] = y[i] / 10.0 *
                  std::exp((1.0L - std::exp((x[i] - t * y[i]) / 10.0)) /
                               (1.0L - std::exp(-1 / 10.0)) +
                           (x[i] - t * y[i]) / 10.0) /
                  (1.0L - std::exp(-1 / 10.0));

    return result;
  };

  // Exact Solution
  SpatialTimeFunc c_ex = DirBC;

  // Gradients of the Exact Solution
  SpatialTimeFunc dc_dx_ex = DirBC_dx;
  SpatialTimeFunc dc_dy_ex = DirBC_dy;
  SpatialTimeFunc dc_dt_ex = DirBC_dt;
};

} // namespace pacs

#endif