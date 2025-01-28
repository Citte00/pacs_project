/**
 * @file DataHeat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Data for the Heat equation problem.
 * @date 2024-12-19
 *
 * @copyright Copyright (c) 2024
 *
 */
#ifndef INCLUDE_PACSHPDG_DATA_DATAHEAT_HPP
#define INCLUDE_PACSHPDG_DATA_DATAHEAT_HPP

#include <array>
#include <iostream>

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Geometry.hpp"

namespace pacs {

struct DataHeat {

  // Definition.
  using SpatialTemporalFunction = GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Real>;
  using SpatialFunction =
      GenFunc<Vector<Real>, Vector<Real>, Vector<Real>>;

  // Geometrical properties
  std::vector<Point> domain = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
  int elements = 125;
  bool meshFromFile = true;
  std::string VTKMeshFileName = "Mesh.vtk";
  std::string meshFileSeq = "meshes/square/square_300.poly";

  // Material properties
  Function<Vector<Real>, Vector<Real>, Vector<Real>, Real> D_ext =
      [](Vector<Real> x, Vector<Real> y, Real t) { return 1.0 + 0.0 * x; };

  // Forcing Term
  bool homog_source_f = false;
  Function<Vector<Real>, Vector<Real>, Vector<Real>, Real, Vector<Real>>
      source_f = [](Vector<Real> x, Vector<Real> y, Real t, Vector<Real> D) {
        Vector<Real> result{x.length};

        for (size_t i = 0; i < x.length; i++)
          result[i] =
              -(std::cos(2.0 * M_PI * x[i]) * std::cos(2 * M_PI * y[i]) + 2) +
              8 * M_PI * M_PI * D[i] * std::cos(2 * M_PI * x[i]) *
                  std::cos(2 * M_PI * y[i]) * (1 - t);

        return result;
      };

  // Boundary Conditions
  Function<Vector<Real>, Vector<Real>, Vector<Real>, Real> DirBC =
      [](Vector<Real> x, Vector<Real> y, Real t) {
        Vector<Real> result{x.length};

        for (size_t i = 0; i < x.length; i++)
          result[i] =
              (std::cos(2.0 * M_PI * x[i]) * std::cos(2.0 * M_PI * y[i]) +
               2.0) *
              (1 - t);

        return result;
      };

  // Gradients of the Boundary Conditions
  Function<Vector<Real>, Vector<Real>, Vector<Real>, Real> DirBC_dx =
      [](Vector<Real> x, Vector<Real> y, Real t) {
        Vector<Real> result{x.size()};

    for (size_t i = 0; i < x.length; i++)
      result[i] = -2.0L * M_PI * std::sin(2.0 * M_PI * x[i]) *
                  std::cos(2.0 * M_PI * y[i]) * (1 - t);

    return result;
  };

  Function<Vector<Real>, Vector<Real>, Vector<Real>, Real> DirBC_dy =
      [](Vector<Real> x, Vector<Real> y, Real t) {
        Vector<Real> result{x.size()};

    for (size_t i = 0; i < x.length; i++)
      result[i] = -2.0L * M_PI * std::cos(2.0 * M_PI * x[i]) *
                  std::sin(2.0 * M_PI * y[i]) * (1 - t);

    return result;
  };

  Function<Vector<Real>, Vector<Real>, Vector<Real>> DirBC_dt =
      [](Vector<Real> x, Vector<Real> y) {
        Vector<Real> result{x.length};

        for (size_t i = 0; i < x.length; i++)
          result[i] = -(
              std::cos(2.0 * M_PI * x[i]) * std::cos(2.0 * M_PI * y[i]) + 2.0);

        return result;
      };

  // Exact Solution
  Function<Vector<Real>, Vector<Real>, Vector<Real>, Real> c_ex = DirBC;

  // Gradients of the Exact Solution
  Function<Vector<Real>, Vector<Real>, Vector<Real>, Real> dc_dx_ex = DirBC_dx;
  Function<Vector<Real>, Vector<Real>, Vector<Real>, Real> dc_dy_ex = DirBC_dy;
  Function<Vector<Real>, Vector<Real>, Vector<Real>> dc_dt_ex = DirBC_dt;

  // Time discretization
  Real t_0 = 0.0;
  Real t_f = 2e-3;
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