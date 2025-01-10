/**
 * @file DataLap.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Data for the Laplace equation problem.
 * @date 2025-01-10
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef INCLUDE_PACSHPDG_DATA_DATALAP_HPP
#define INCLUDE_PACSHPDG_DATA_DATALAP_HPP

#include <array>
#include <iostream>

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Fem/GeneralFunctor.hpp"
#include "../Geometry.hpp"

namespace pacs {

struct DataLap {

  // Definition.
  using SpatialTemporalFunction =
      GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Real>;
  using SpatialFunction = GenFunc<Vector<Real>, Vector<Real>, Vector<Real>>;

  // Geometrical properties
  std::vector<Point> domain = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
  int N = 300;
  bool meshFromFile = true;
  std::string VTKMeshFileName = "Mesh.vtk";
  std::string meshFileSeq = "meshes/square/square_300.poly";

  // Material properties
  SpatialFunction D_ext = [](Vector<Real> x, Vector<Real> y) {
    return 0.0 + 0.0 * x;
  };

  // Forcing Term
  bool homog_source_f = false;
  GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Vector<Real>>
      source_f = [](Vector<Real> x, Vector<Real> y, Vector<Real> D) {
        Vector<Real> result{x.length};

        for (size_t i = 0; i < x.length; i++)
          result[i] = 8 * M_PI * M_PI * D[i] * std::cos(2 * M_PI * x[i]) *
                      std::cos(2 * M_PI * y[i]);

        return result;
      };

  // Boundary Conditions
  SpatialFunction DirBC = [](Vector<Real> x, Vector<Real> y) {
    Vector<Real> result{x.length};

    for (size_t i = 0; i < x.length; i++)
      result[i] =
          (std::cos(2.0 * M_PI * x[i]) * std::cos(2.0 * M_PI * y[i]) + 2.0);

    return result;
  };

  // Gradients of the Boundary Conditions
  SpatialFunction DirBC_dx = [](Vector<Real> x, Vector<Real> y) {
    Vector<Real> result{x.length};

    for (size_t i = 0; i < x.length; i++)
      result[i] = -2.0L * M_PI * std::sin(2.0 * M_PI * x[i]) *
                  std::cos(2.0 * M_PI * y[i]);

    return result;
  };

  SpatialFunction DirBC_dy = [](Vector<Real> x, Vector<Real> y) {
    Vector<Real> result{x.length};

    for (size_t i = 0; i < x.length; i++)
      result[i] = -2.0L * M_PI * std::cos(2.0 * M_PI * x[i]) *
                  std::sin(2.0 * M_PI * y[i]);

    return result;
  };

  // Exact Solution
  SpatialFunction c_ex = DirBC;

  // Gradients of the Exact Solution
  SpatialFunction dc_dx_ex = DirBC_dx;
  SpatialFunction dc_dy_ex = DirBC_dy;

  // Space discretization
  size_t degree = 4;
  Real penalty_coeff = 10.0;

  // Visualization settings
  bool PlotExact = true;
  bool PlotGridSol = true;
  bool PlotIniCond = true;
  int VisualizationStep = 10;
  int NqnVisualization = 5;

  // Save solution settings
  double SaveSolutionStep = 0.05;
  std::string VTKpFileName = "Distribution_p_" + std::to_string(degree) + "_t_";
};

} // namespace pacs

#endif