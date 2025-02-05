/**
 * @file test_fit.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief
 * @date 2024-07-08
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <PacsHPDG.hpp>

#include <iostream>

int main() {

  using namespace pacs;

  // Domain;
  Point a{0.0, 0.0};
  Point b{0.0, 1.0};
  Point c{1.0, 1.0};
  Point d{0.0, 1.0};

  Polygon domain{{a, b, c, d}};

  // Mesh.
  Mesh mesh{domain, mesh_diagram("meshes/square/square_125.poly")};

  // Constructing some data.
  Vector<Real> x{4};
  Vector<Real> y{4};

  x[0] = 1.0L;
  x[1] = 2.0L;
  x[2] = 3.0L;
  x[3] = 4.0L;

  y[0] = 1.0L;
  y[1] = 4.0L;
  y[2] = 9.0L;
  y[3] = 16.0L;

  // Polynomial fit.
  LaplaceEstimator estimator(mesh);
  std::cout << estimator.polyfit(x, y, 2) << std::endl;
}