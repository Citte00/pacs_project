/**
 * @file DataFKPP.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Data struct for the Fisher-KPP equation problem.
 * @date 2024-12-27
 *
 * @copyright Copyright (c) 2024
 *
 */
#ifndef INCLUDE_PACSHPDG_DATA_DATAFKPP_HPP
#define INCLUDE_PACSHPDG_DATA_DATAFKPP_HPP

#include <array>
#include <iostream>

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"

namespace pacs {

struct DataFKPP {

  using SpatialFunc = Functor<Vector<Real>, Vector<Real>, Vector<Real>>;
  using SpatialTimeFunc =
      Functor<Vector<Real>, Vector<Real>, Vector<Real>, Real>;
  using SourceFunc =
      Functor<Vector<Real>, Vector<Real>, Vector<Real>, Real, Vector<Real>, Vector<Real>>;

  // Geometrical properties
  std::vector<Point> domain = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
  int elements = 125;

  // Material properties
  SpatialTimeFunc D_ext = [](Vector<Real> x, Vector<Real> y, Real t) {
    return 1.0 + 0.0 * x;
  };
  SpatialTimeFunc alpha = [](Vector<Real> x, Vector<Real> y, Real t) {
    return 1.0 + 0.0 * x;
  };

  // Forcing Term
  bool homog_source_f = false;
  SourceFunc source_f =
      [](Vector<Real> x, Vector<Real> y, Real t, Vector<Real> D,
         Vector<Real> alpha) {
        Vector<Real> result{x.size()};

        for (size_t i = 0; i < x.length; i++)
          result[i] =
              -(1.0L - std::exp(-100.0L * x[i])) / (1.0L - std::exp(-100.0)) *
                  std::sin(M_PI * y[i]) * (1.0L - x[i]) +
              D[i] *
                  ((10000.0L * std::exp(-100.0L * x[i])) /
                       (1.0L - std::exp(-100.0)) * std::sin(M_PI * y[i]) *
                       (1.0L - x[i]) +
                   (200.0L * std::exp(-100.0L * x[i])) /
                       (1.0L - std::exp(-100.0)) * std::sin(M_PI * y[i]) +
                   M_PI * M_PI * (1.0L - std::exp(-100.0L * x[i])) /
                       (1.0L - std::exp(-100.0)) * std::sin(M_PI * y[i]) *
                       (1.0L - x[i])) *
                  (1 - t) -
              alpha[i] * (1 - t) * (1.0L - std::exp(-100.0L * x[i])) /
                  (1.0L - std::exp(-100.0)) * std::sin(M_PI * y[i]) *
                  (1.0L - x[i]) *
                  (1 - (1 - t) * (1.0L - std::exp(-100.0L * x[i])) /
                           (1.0L - std::exp(-100.0)) * std::sin(M_PI * y[i]) *
                           (1.0L - x[i]));

        return result;
      };

  // Boundary Conditions
  SpatialTimeFunc DirBC =
      [](Vector<Real> x, Vector<Real> y, Real t) {
        Vector<Real> result{x.length};

        for (size_t i = 0; i < x.length; i++)
          result[i] = (1.0L - std::exp(-100.0L * x[i])) /
                      (1.0L - std::exp(-100.0)) * std::sin(M_PI * y[i]) *
                      (1.0L - x[i]) * (1 - t);

        return result;
      };

  // Gradients of the Boundary Conditions
  SpatialTimeFunc DirBC_dx =
      [](Vector<Real> x, Vector<Real> y, Real t) {
        Vector<Real> result{x.size()};

        for (size_t i = 0; i < x.length; i++)
          result[i] = std::sin(M_PI * y[i]) * (1 - t) /
                      (1.0L - std::exp(-100.0)) *
                      ((100.0L * std::exp(-100.0L * x[i])) * (1.0L - x[i]) -
                       (1.0L - std::exp(-100.0L * x[i])));

        return result;
      };

  SpatialTimeFunc DirBC_dy =
      [](Vector<Real> x, Vector<Real> y, Real t) {
        Vector<Real> result{x.size()};

        for (size_t i = 0; i < x.length; i++)
          result[i] = M_PI * (1.0L - std::exp(-100.0L * x[i])) /
                      (1.0L - std::exp(-100.0)) * std::cos(M_PI * y[i]) *
                      (1.0L - x[i]) * (1 - t);

        return result;
      };

  SpatialFunc DirBC_dt =
      [](Vector<Real> x, Vector<Real> y) {
        Vector<Real> result{x.length};

        for (size_t i = 0; i < x.length; i++)
          result[i] = -(1.0L - std::exp(-100.0L * x[i])) /
                      (1.0L - std::exp(-100.0)) * std::sin(M_PI * y[i]) *
                      (1.0L - x[i]);

        return result;
      };

  // Exact Solution
  SpatialTimeFunc c_ex = DirBC;

  // Gradients of the Exact Solution
  SpatialTimeFunc dc_dx_ex = DirBC_dx;
  SpatialTimeFunc dc_dy_ex = DirBC_dy;
  SpatialFunc dc_dt_ex = DirBC_dt;

  // Time discretization
  Real t_0 = 0.0;
  Real t_f = 1e-3;
  Real dt = 1e-5;
  Real theta = 0.5;

  // Space discretization
  size_t degree = 3;
  Real penalty_coeff = 10.0;

  // Visualization settings
  int VisualizationStep = 10;
  int NqnVisualization = 5;
};

} // namespace pacs

#endif