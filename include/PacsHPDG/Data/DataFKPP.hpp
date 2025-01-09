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
#include "../Fem/GeneralFunctor.hpp"
#include "../Geometry.hpp"

namespace pacs {

struct DataFKPP {

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
  SpatialTemporalFunction D_ext = [](Vector<Real> x, Vector<Real> y, Real t) {
    return 1.0 + 0.0 * x;
  };
  SpatialTemporalFunction alpha = [](Vector<Real> x, Vector<Real> y, Real t) {
    return 1.0 + 0.0 * x;
  };

  // Forcing Term
  bool homog_source_f = false;
  GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Real, Vector<Real>,
          Vector<Real>>
      source_f = [](Vector<Real> x, Vector<Real> y, Real t, Vector<Real> D,
                    Vector<Real> alpha) {
        Vector<Real> result{x.length};

        for (size_t i = 0; i < x.length; i++)
          result[i] =
              -(std::cos(2.0 * M_PI * x[i]) * std::cos(2 * M_PI * y[i]) + 2) +
              8 * M_PI * M_PI * D[i] * std::cos(2.0 * M_PI * x[i]) * std::cos(2.0 * M_PI * y[i]) * (1 - t) - 
              alpha[i] * (1-t) * (std::cos(2.0*M_PI*x[i]) * std::cos(2*M_PI*y[i])+2) * (1 - (1-t) *(std::cos(2.0*M_PI*x[i]) * std::cos(2.0*M_PI*y[i])+2));

        return result;
      };

  // Boundary Conditions
  SpatialTemporalFunction DirBC =
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
  SpatialTemporalFunction DirBC_dx =
      [](Vector<Real> x, Vector<Real> y, Real t) {
        Vector<Real> result{x.length};

        for (size_t i = 0; i < x.length; i++)
          result[i] = -2.0L * M_PI * std::sin(2.0 * M_PI * x[i]) *
                      std::cos(2.0 * M_PI * y[i]) * (1 - t);

        return result;
      };

  SpatialTemporalFunction DirBC_dy =
      [](Vector<Real> x, Vector<Real> y, Real t) {
        Vector<Real> result{x.length};

        for (size_t i = 0; i < x.length; i++)
          result[i] = -2.0L * M_PI * std::cos(2.0 * M_PI * x[i]) *
                      std::sin(2.0 * M_PI * y[i]) * (1 - t);

        return result;
      };

  SpatialFunction DirBC_dt = [](Vector<Real> x, Vector<Real> y) {
    Vector<Real> result{x.length};

    for (size_t i = 0; i < x.length; i++)
      result[i] =
          -(std::cos(2.0 * M_PI * x[i]) * std::cos(2.0 * M_PI * y[i]) + 2.0);

    return result;
  };

  // Exact Solution
  SpatialTemporalFunction c_ex = DirBC;

  // Gradients of the Exact Solution
  SpatialTemporalFunction dc_dx_ex = DirBC_dx;
  SpatialTemporalFunction dc_dy_ex = DirBC_dy;
  SpatialFunction dc_dt_ex = DirBC_dt;

  // Time discretization
  Real t_0 = 0.0;
  Real t_f = 2.0;
  Real dt = 1e-2;
  Real theta = 0.5;

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