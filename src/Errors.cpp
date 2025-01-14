/**
 * @file Errors.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

namespace pacs {

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Error structure.
     * 
     * @param mesh Mesh.
     * @param matrices Needed matrices (A, M).
     * @param numerical Numerical solution.
     * @param exact Exact solution.
     */
void LaplaceError::computeErrors(const DataLaplace &data, const Mesh &mesh, const Laplace &laplacian, const Vector<Real> &numerical) {
#ifndef NVERBOSE
  std::cout << "Evaluating errors." << std::endl;
#endif

  // Matrices.
  Sparse<Real> mass = laplacian.M();
  Sparse<Real> dg_stiff = laplacian.DG();

  // Mass blocks.
  auto blocks = laplacian.block_mass(mesh);

  // Error vector.
  Vector<Real> u_modals = modal(mesh, data.c_ex);
  Vector<Real> u_coeff = solve(mass, u_modals, blocks, DB);

  Vector<Real> error = u_coeff - numerical;

  // DG Error.
  this->m_dg_error = std::sqrt(dot(error, dg_stiff * error));

  // L2 Error.
  this->m_l2_error = std::sqrt(dot(error, mass * error));

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.reserve(mesh.elements.size());
  starts.emplace_back(0);

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

  // Quadrature nodes.
  std::vector<std::size_t> nqn(mesh.elements.size(), 0);
  std::transform(mesh.elements.begin(), mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

  // Sizes.
  Vector<Real> sizes{mesh.elements.size()};

  for (std::size_t j = 0; j < sizes.length; ++j) {
    Element element{mesh.elements[j]};

    for (const auto &p : element.element.points)
      for (const auto &q : element.element.points)
        sizes[j] = (distance(p, q) > sizes[j]) ? distance(p, q) : sizes[j];
  }

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // 2D quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

    // Local dofs.
    std::size_t element_dofs = mesh.elements[j].dofs();

    // Global matrix indices.
    std::vector<std::size_t> indices;
    indices.reserve(element_dofs);

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

      // Jacobian's determinant.
      Real jacobian_det =
          jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);

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

      // Weights scaling.
      Vector<Real> scaled = jacobian_det * weights_2d;

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      // Solutions.
      Vector<Real> u = data.c_ex(physical_x, physical_y);
      Vector<Real> uh = phi * numerical(indices);

      Vector<Real> grad_u = data.dc_dx_ex(physical_x, physical_y) +
                            data.dc_dy_ex(physical_x, physical_y);
      Vector<Real> grad_uh = (gradx_phi + grady_phi) * numerical(indices);

      // Local L2 error.
      this->m_l2_errors[j] += dot(scaled, (u - uh) * (u - uh));

      // Local H1 error.
      this->m_h1_errors[j] +=
          dot(scaled, (grad_u - grad_uh) * (grad_u - grad_uh));
    }

    this->m_l2_errors[j] = std::sqrt(this->m_l2_errors[j]);
    this->m_h1_errors[j] = std::sqrt(this->m_h1_errors[j]);
  }
};

} // namespace pacs