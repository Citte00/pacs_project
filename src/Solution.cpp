/**
 * @file Solution.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief
 * @date 2024-05-11
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <PacsHPDG.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace pacs {

/**
 * @brief Constructs a new solution object.
 *
 * @param mesh Mesh.
 * @param numerical Numerical solution.
 * @param exact Exact solution.
 */
Solution::Solution(const Mesh &mesh, const Vector<Real> &numerical,
                   const BiFunctor &exact)
    : x{mesh.entries}, y{mesh.entries}, numerical{mesh.entries},
      exact{mesh.entries} {

  // Number of quadrature nodes.
  std::size_t nqn = GAUSS_ORDER;

  // Quadrature nodes.
  auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn);

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.emplace_back(0);

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

  // Local vectors indices.
  std::vector<std::size_t> local_indices;

  for (std::size_t h = 0; h < nqn * nqn; ++h)
    local_indices.emplace_back(h);

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // Local dofs.
    std::size_t element_dofs = mesh.elements[j].dofs();

    // Global matrix indices.
    std::vector<std::size_t> indices;

    for (std::size_t k = 0; k < element_dofs; ++k)
      indices.emplace_back(starts[j] + k);

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
      Vector<Real> local_exact = exact(physical_x, physical_y);

      // Writing.
      this->x(local_indices, physical_x);
      this->y(local_indices, physical_y);
      this->numerical(local_indices, local_numerical);
      this->exact(local_indices, local_exact);

      // Local indices update.
      for (auto &index : local_indices)
        index += nqn * nqn;
    }
  }
}

/**
 * @brief Constructs a new Solution object for the heat equation.
 *
 * @param data Data struct.
 * @param mesh Mesh.
 * @param numerical Numerical solution.
 * @param t Time step.
 */
Solution::Solution(const DataHeat &data, const Mesh &mesh,
                   const Vector<Real> &numerical, const Real &t)
    : x{(data.NqnVisualization + 0) * data.NqnVisualization * mesh.entries},
      y{(data.NqnVisualization + 0) * data.NqnVisualization * mesh.entries},
      numerical{(data.NqnVisualization + 0) * data.NqnVisualization *
                mesh.entries},
      exact{(data.NqnVisualization + 0) * data.NqnVisualization *
            mesh.entries} {

  // Number of quadrature nodes.
  std::size_t nqn = data.NqnVisualization;

  // Quadrature nodes.
  auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn);
  auto [nodes_1d, weights_1d] = quadrature_1d(nqn);

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.emplace_back(0);

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

  // Local vectors indices.
  std::vector<std::size_t> local_indices;

  for (std::size_t h = 0; h < nqn * nqn; ++h)
    local_indices.emplace_back(h);

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // Local dofs.
    std::size_t element_dofs = mesh.elements[j].dofs();

    // Global matrix indices.
    std::vector<std::size_t> indices;

    for (std::size_t k = 0; k < element_dofs; ++k)
      indices.emplace_back(starts[j] + k);

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
      Vector<Real> local_exact = data.c_ex(physical_x, physical_y, t);

      // Writing.
      this->x(local_indices, physical_x);
      this->y(local_indices, physical_y);
      this->numerical(local_indices, local_numerical);
      this->exact(local_indices, local_exact);

      // Local indices update.
      for (auto &index : local_indices)
        index += nqn * nqn;
    }

    // Local neighbour indices.
    /*std::vector<std::size_t> neigh_indices;

    for (std::size_t i = 0; i < nqn+1; i++)
      neigh_indices.emplace_back(*local_indices.begin() + i);

    // Face integrals.

    // Element's neighbours.
    std::vector<std::array<int, 3>> element_neighbours = mesh.neighbours[j];

    // Edges.
    std::vector<Segment> edges{polygon.edges()};

    // Loop over faces.
    for (std::size_t k = 0; k < element_neighbours.size(); ++k) {

      // Edge geometry.
      Segment segment{edges[k]};

      // Edge's normal.
      Vector<Real> edge_vector{2};

      edge_vector[0] = segment[1][0] - segment[0][0];
      edge_vector[1] = segment[1][1] - segment[0][1];

      Vector<Real> normal_vector{2};

      normal_vector[0] = edge_vector[1];
      normal_vector[1] = -edge_vector[0];

      normal_vector /= norm(normal_vector);

      // Jacobian.
      Matrix<Real> jacobian{2, 2};

      jacobian(0, 0) = segment[1][0] - segment[0][0];
      jacobian(0, 1) = 0.5 * (segment[1][0] - segment[0][0]);
      jacobian(1, 0) = segment[1][1] - segment[0][1];
      jacobian(1, 1) = 0.5 * (segment[1][1] - segment[0][1]);

      // Translation.
      Vector<Real> translation{2};

      translation[0] = segment[0][0];
      translation[1] = segment[0][1];

      // Physical nodes.
      Vector<Real> physical_x{nodes_1d.length};
      Vector<Real> physical_y{nodes_1d.length};

      for (std::size_t l = 0; l < nodes_1d.length; ++l) {
        Vector<Real> node{2};

        node[0] = nodes_1d[l];

        Vector<Real> transformed = jacobian * node + translation;

        physical_x[l] = transformed[0];
        physical_y[l] = transformed[1];
      }

      // Basis functions.
      Matrix<Real> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];

      // Numerical solution.
      Vector<Real> local_numerical = phi * numerical(indices);

      // Exact solution.
      Vector<Real> local_exact = data.c_ex(physical_x, physical_y, t);

      // Writing.
      this->x(neigh_indices, physical_x);
      this->y(neigh_indices, physical_y);
      this->numerical(neigh_indices, local_numerical);
      this->exact(neigh_indices, local_exact);

      // Local indices update.
      for (auto &index : neigh_indices)
        index += nqn;

    }

    for (auto &index : local_indices)
      index += element_neighbours.size() * nqn;*/
  }
}

/**
 * @brief Constructs a new Solution object for the Fisher-KPP equation.
 *
 * @param data Data struct.
 * @param mesh Mesh.
 * @param numerical Numerical solution.
 * @param t Time step.
 */
Solution::Solution(const DataFKPP &data, const Mesh &mesh,
                   const Vector<Real> &numerical, const Real &t)
    : x{data.NqnVisualization * data.NqnVisualization * mesh.entries},
      y{data.NqnVisualization * data.NqnVisualization * mesh.entries},
      numerical{data.NqnVisualization * data.NqnVisualization * mesh.entries},
      exact{data.NqnVisualization * data.NqnVisualization * mesh.entries} {

  // Number of quadrature nodes.
  std::size_t nqn = data.NqnVisualization;

  // Quadrature nodes.
  auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn);
  auto [nodes_1d, weights_1d] = quadrature_1d(nqn);

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.emplace_back(0);

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

  // Local vectors indices.
  std::vector<std::size_t> local_indices;

  for (std::size_t h = 0; h < nqn * nqn; ++h)
    local_indices.emplace_back(h);

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // Local dofs.
    std::size_t element_dofs = mesh.elements[j].dofs();

    // Global matrix indices.
    std::vector<std::size_t> indices;

    for (std::size_t k = 0; k < element_dofs; ++k)
      indices.emplace_back(starts[j] + k);

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
      Vector<Real> local_exact = data.c_ex(physical_x, physical_y, t);

      // Writing.
      this->x(local_indices, physical_x);
      this->y(local_indices, physical_y);
      this->numerical(local_indices, local_numerical);
      this->exact(local_indices, local_exact);

      // Local indices update.
      for (auto &index : local_indices)
        index += nqn * nqn;
    }
  }
}

// OUTPUT.

/**
 * @brief Outputs the solution to a .sol file.
 *
 * @param filename Filename.
 */
void Solution::write(const std::string &filename) {
  // File loading.
  std::ofstream file{filename};

  file << "@ " << filename << "\n";
  file << "@ contplot.py readable mesh\n";
  file << "@ Structure: [x, y, numerical(x, y), exact(x, y)].\n";

  for (std::size_t j = 0; j < this->x.length; ++j) {
    file << this->x[j] << ",";
    file << this->y[j] << ",";
    file << this->numerical[j] << ",";
    file << this->exact[j] << "\n";
  }

  file.close();
}

/**
 * @brief Outputs the solution to a .vtk file
 *
 * @param filename Output file name.
 */
void Solution::writeVTK(const Mesh &mesh, const std::string &filename) {

  // Open the file.
  std::ofstream vtkFile{filename};

  // Check if the vtkFile was opened.
  if (!vtkFile.is_open()) {
    std::cerr << "Error opening vtkFile " << filename << std::endl;
    return;
  }

  // Headers.
  vtkFile << "# vtk DataFile Version 1.0" << std::endl;
  vtkFile << "Comment" << std::endl;
  vtkFile << "ASCII" << std::endl;
  vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

  // Points.
  vtkFile << "POINTS " << this->x.length << " float" << std::endl;
  for (std::size_t j = 0; j < this->x.length; ++j) {
    vtkFile << this->x[j] << " ";
    vtkFile << this->y[j] << " ";
    vtkFile << "0.0" << std::endl;
  }

  // Cells.
  size_t numCells = mesh.elements.size();
  size_t conn = 0;

  for (int i = 0; i < numCells; i++)
    conn += mesh.elements[i].edges.size();

  conn += numCells;

  vtkFile << "POLYGONS " << numCells << " " << conn << std::endl;
  // Loop over the cells.
  for (int i = 0; i < numCells; i++) {
    std::vector<std::array<int, 3>> element_neighbours = mesh.neighbours[i];
    vtkFile << element_neighbours.size() << " ";

    for (int k = 0; k < element_neighbours.size(); k++)
      vtkFile << element_neighbours[k][1] << " ";

    vtkFile << std::endl;
  }

  // Cell types.
  vtkFile << "CELL_TYPES " << numCells << std::endl;
  for (std::size_t j = 0; j < numCells; ++j) {
    vtkFile << "7" << std::endl;
  }

  // Point data.
  vtkFile << "POINT_DATA " << this->x.length << std::endl;
  vtkFile << "SCALARS " << "c_h " << "float 1" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;

  for (std::size_t j = 0; j < this->x.length; ++j)
    vtkFile << this->numerical[j] << std::endl;

  // Closing the file.
  vtkFile.close();
}

} // namespace pacs