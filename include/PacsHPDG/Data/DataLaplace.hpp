/**
 * @file DataHeat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Data for the Heat equation problem.
 * @date 2024-12-19
 *
 * @copyright Copyright (c) 2024
 *
 */
#ifndef INCLUDE_PACSHPDG_DATA_DATALAPLACE_HPP
#define INCLUDE_PACSHPDG_DATA_DATALAPLACE_HPP

#include <array>
#include <iostream>

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"

namespace pacs {

struct DataLaplace {

  // Defining some aliases.
  using SpatialFunc = Functor<Vector<Real>, Vector<Real>, Vector<Real>>;

  // Geometrical properties
  std::vector<Point> domain = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
  int elements = 125;

  // Forcing Term
  SpatialFunc source_f =
      [](Vector<Real> x, Vector<Real> y) {
        Vector<Real> result{x.size()};

        for (size_t i = 0; i < x.size(); i++)
          result[i] = 8 * M_PI * M_PI * std::cos(2 * M_PI * x[i]) *
                  std::cos(2 * M_PI * y[i]);

        return result;
      };

  // Boundary Conditions
  SpatialFunc DirBC = [](Vector<Real> x, Vector<Real> y) {
    Vector<Real> result{x.size()};

    for (size_t i = 0; i < x.size(); i++)
      result[i] = std::cos(2.0 * M_PI * x[i]) * std::cos(2.0 * M_PI * y[i]);

    return result;
  };

  // Gradients of the Boundary Conditions
  SpatialFunc DirBC_dx = [](Vector<Real> x, Vector<Real> y) {
    Vector<Real> result{x.size()};

    for (size_t i = 0; i < x.size(); i++)
      result[i] = -2.0L * M_PI * std::sin(2.0 * M_PI * x[i]) *
                  std::cos(2.0 * M_PI * y[i]);

    return result;
  };

  SpatialFunc DirBC_dy = [](Vector<Real> x, Vector<Real> y) {
    Vector<Real> result{x.size()};

    for (size_t i = 0; i < x.size(); i++)
      result[i] = -2.0L * M_PI * std::cos(2.0 * M_PI * x[i]) *
                  std::sin(2.0 * M_PI * y[i]);

    return result;
  };

  // Exact Solution
  SpatialFunc c_ex = DirBC;

  // Gradients of the Exact Solution
  SpatialFunc dc_dx_ex = DirBC_dx;
  SpatialFunc dc_dy_ex = DirBC_dy;
  
  // Space discretization
  size_t degree = 2;
  Real penalty_coeff = 10.0;

  int VisualizationStep = 10;
  int NqnVisualization = 5;
};

} // namespace pacs

#endif