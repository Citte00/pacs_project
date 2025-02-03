/**
 * @file Solution.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2025-01-20
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <PacsHPDG.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace pacs {

/**
 * @brief Construct a new Laplace solution object.
 * 
 * @param data Laplace equation data struct.
 * @param mesh Mesh struct.
 * @param numerical numerical solution.
 */
void LaplaceSolution::computeSolution(const DataLaplace &data, const Mesh &mesh, const Vector<Real> &numerical) {

  // Number of quadrature nodes.
  std::size_t nqn = data.NqnVisualization;

  // Quadrature nodes.
  auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn);

  // Starting indices.
  std::vector<std::size_t> starts(mesh.elements.size());
  starts[0] = 0;

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts[j] = starts[j - 1] + mesh.elements[j - 1].dofs();

  // Local vectors indices.
  std::vector<std::size_t> local_indices(nqn*nqn);

  for (std::size_t h = 0; h < nqn * nqn; ++h)
    local_indices[h] = h;

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // Global matrix indices.
    std::vector<std::size_t> indices(mesh.elements[j].dofs());

    for (std::size_t k = 0; k < mesh.elements[j].dofs(); ++k)
      indices[k] = starts[j] + k;

    // Polygon.
    Polygon polygon = mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Loop over the sub-triangulation.
    for (std::size_t k = 0; k < triangles.size(); ++k) {

      // Triangle.
      Polygon triangle = triangles[k];

      // Jacobian.
      Matrix<Real> jacobian{2, 2};

      jacobian(0, 0) = triangle.points[1][0] - triangle.points[0][0];
      jacobian(0, 1) = triangle.points[2][0] - triangle.points[0][0];
      jacobian(1, 0) = triangle.points[1][1] - triangle.points[0][1];
      jacobian(1, 1) = triangle.points[2][1] - triangle.points[0][1];

      // Translation.
      Vector<Real> translation{2};

      translation[0] = triangle.points[0][0];
      translation[1] = triangle.points[0][1];

      // Physical nodes.
      Vector<Real> physical_x{nodes_x_2d.length};
      Vector<Real> physical_y{nodes_y_2d.length};

      for (std::size_t l = 0; l < physical_x.length; ++l) {
        Vector<Real> node{2};

        node[0] = nodes_x_2d[l];
        node[1] = nodes_y_2d[l];

        Vector<Real> transformed = jacobian * node + translation;

        physical_x[l] = transformed[0];
        physical_y[l] = transformed[1];
      }

      // Basis functions.
      Matrix<Real> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];

      // Numerical solution.
      Vector<Real> local_numerical = phi * numerical(indices);

      // Exact solution.
      Vector<Real> local_exact = data.c_ex(physical_x, physical_y);

      // Writing.
      this->m_x(local_indices, physical_x);
      this->m_y(local_indices, physical_y);
      this->m_numerical(local_indices, local_numerical);
      this->m_exact(local_indices, local_exact);

      // Local indices update.
      for (auto &index : local_indices)
        index += nqn * nqn;
    }
  }
};

/**
 * @brief Construct a new Heat solution object.
 *
 * @param data Heat equation data struct.
 * @param mesh Mesh struct.
 * @param numerical numerical solution.
 */
void HeatSolution::computeSolution(const DataHeat &data, const Mesh &mesh,
                                   const Heat &heat, const Vector<Real> &ch) {
#ifndef NVERBOSE
  std::cout << "Evaluating solutions." << std::endl;
#endif

  // Number of quadrature nodes.
  std::size_t nqn = GAUSS_ORDER;

  // Quadrature nodes.
  auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn);

  // Starting indices.
  std::vector<std::size_t> starts(mesh.elements.size());
  starts[0] = 0;

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts[j] = starts[j - 1] + mesh.elements[j - 1].dofs();

  // Local vectors indices.
  std::vector<std::size_t> local_indices(nqn * nqn);

  for (std::size_t h = 0; h < nqn * nqn; ++h)
    local_indices[h] = h;

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // Global matrix indices.
    std::vector<std::size_t> indices(mesh.elements[j].dofs());

    for (std::size_t k = 0; k < mesh.elements[j].dofs(); ++k)
      indices[k] = starts[j] + k;

    // Polygon.
    Polygon polygon = mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Loop over the sub-triangulation.
    for (std::size_t k = 0; k < triangles.size(); ++k) {

      // Triangle.
      Polygon triangle = triangles[k];

      // Jacobian.
      Matrix<Real> jacobian{2, 2};

      jacobian(0, 0) = triangle.points[1][0] - triangle.points[0][0];
      jacobian(0, 1) = triangle.points[2][0] - triangle.points[0][0];
      jacobian(1, 0) = triangle.points[1][1] - triangle.points[0][1];
      jacobian(1, 1) = triangle.points[2][1] - triangle.points[0][1];

      // Translation.
      Vector<Real> translation{2};

      translation[0] = triangle.points[0][0];
      translation[1] = triangle.points[0][1];

      // Physical nodes.
      Vector<Real> physical_x{nodes_x_2d.length};
      Vector<Real> physical_y{nodes_y_2d.length};

      for (std::size_t l = 0; l < physical_x.length; ++l) {
        Vector<Real> node{2};

        node[0] = nodes_x_2d[l];
        node[1] = nodes_y_2d[l];

        Vector<Real> transformed = jacobian * node + translation;

        physical_x[l] = transformed[0];
        physical_y[l] = transformed[1];
      }

      // Basis functions.
      Matrix<Real> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];

      // Numerical solution.
      Vector<Real> local_numerical = phi * ch(indices);

      // Exact solution.
      Vector<Real> local_exact = data.c_ex(physical_x, physical_y, heat.t());

      // Writing.
      this->m_x(local_indices, physical_x);
      this->m_y(local_indices, physical_y);
      this->m_numerical(local_indices, local_numerical);
      this->m_exact(local_indices, local_exact);

      // Local indices update.
      for (auto &index : local_indices)
        index += nqn * nqn;
    }
  }
};

/**
 * @brief Construct a new Fisher solution object.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh struct.
 * @param fisher Fisher equation object.
 */
void FisherSolution::computeSolution(const DataFKPP &data, const Mesh &mesh,
                                     const Fisher &fisher,
                                     const Vector<Real> &ch) {

  // Number of quadrature nodes.
  std::size_t nqn = GAUSS_ORDER;

  // Quadrature nodes.
  auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn);

  // Starting indices.
  std::vector<std::size_t> starts(mesh.elements.size());
  starts[0] = 0;

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts[j] = starts[j - 1] + mesh.elements[j - 1].dofs();

  // Local vectors indices.
  std::vector<std::size_t> local_indices(nqn * nqn);

  for (std::size_t h = 0; h < nqn * nqn; ++h)
    local_indices[h] = h;

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // Global matrix indices.
    std::vector<std::size_t> indices(mesh.elements[j].dofs());

    for (std::size_t k = 0; k < mesh.elements[j].dofs(); ++k)
      indices[k] = starts[j] + k;

    // Polygon.
    Polygon polygon = mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Loop over the sub-triangulation.
    for (std::size_t k = 0; k < triangles.size(); ++k) {

      // Triangle.
      Polygon triangle = triangles[k];

      // Jacobian.
      Matrix<Real> jacobian{2, 2};

      jacobian(0, 0) = triangle.points[1][0] - triangle.points[0][0];
      jacobian(0, 1) = triangle.points[2][0] - triangle.points[0][0];
      jacobian(1, 0) = triangle.points[1][1] - triangle.points[0][1];
      jacobian(1, 1) = triangle.points[2][1] - triangle.points[0][1];

      // Translation.
      Vector<Real> translation{2};

      translation[0] = triangle.points[0][0];
      translation[1] = triangle.points[0][1];

      // Physical nodes.
      Vector<Real> physical_x{nodes_x_2d.length};
      Vector<Real> physical_y{nodes_y_2d.length};

      for (std::size_t l = 0; l < physical_x.length; ++l) {
        Vector<Real> node{2};

        node[0] = nodes_x_2d[l];
        node[1] = nodes_y_2d[l];

        Vector<Real> transformed = jacobian * node + translation;

        physical_x[l] = transformed[0];
        physical_y[l] = transformed[1];
      }

      // Basis functions.
      Matrix<Real> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];

      // Numerical solution.
      Vector<Real> local_numerical = phi * ch(indices);

      // Exact solution.
      Vector<Real> local_exact = data.c_ex(physical_x, physical_y, fisher.t());

      // Writing.
      this->m_x(local_indices, physical_x);
      this->m_y(local_indices, physical_y);
      this->m_numerical(local_indices, local_numerical);
      this->m_exact(local_indices, local_exact);

      // Local indices update.
      for (auto &index : local_indices)
        index += nqn * nqn;
    }
  }
};

/**
 * @brief Outputs the solution to a .sol file.
 *
 * @param filename Filename.
 */
void LaplaceSolution::write(const std::string &filename) {

  // File loading.
  std::ofstream file{filename};

  file << "@ " << filename << "\n";
  file << "@ contplot.py readable mesh\n";
  file << "@ Structure: [x, y, numerical(x, y), exact(x, y)].\n";

  for (std::size_t j = 0; j < this->m_x.length; ++j) {
    file << this->m_x[j] << ",";
    file << this->m_y[j] << ",";
    file << this->m_numerical[j] << ",";
    file << this->m_exact[j] << "\n";
  }

  file.close();
};

} // namespace pacs