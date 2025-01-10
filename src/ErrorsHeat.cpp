/**
 * @file ErrorsHeat.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief
 * @date 2024-12-05
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <PacsHPDG.hpp>

namespace pacs {

/**
 * @brief Construct a new Heat Error object.
 *
 * @param data Data struct.
 * @param mesh Mesh struct.
 * @param matrices Matrices for the l2- and dg-error computation [mass,
 * dg_stiff].
 * @param numerical Numerical solution.
 * @param t Time step.
 */
HeatError::HeatError(const DataHeat &data, const Mesh &mesh,
                     const std::array<Sparse<Real>, 2> &matrices,
                     const Vector<Real> &numerical, const Real &t)
    : Error{mesh.elements.size()} {
#ifndef NVERBOSE
  std::cout << "Evaluating errors." << std::endl;
#endif

  // Matrices.
  auto [mass, dg_stiff] = matrices;

  // Mass blocks.
  Laplace eq{mesh};
  auto blocks = eq.block_mass(mesh);

  // Error vector.
  Vector<Real> u_modals = evaluateCoeff(mesh, data.c_ex, t);
  Vector<Real> u_coeff = solve(mass, u_modals, blocks, DB);

  Vector<Real> error = u_coeff - numerical;

  // DG Error.
  this->dg_error = std::sqrt(dot(error, dg_stiff * error));

  // L2 Error.
  this->l2_error = std::sqrt(dot(error, mass * error));

  // Dofs.
  this->dofs = mesh.dofs();

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.emplace_back(0);

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

  // Quadrature nodes.
  auto [nodes_1d, weights_1d] = quadrature_1d(GAUSS_ORDER);
  auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(GAUSS_ORDER);

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
      Vector<Real> u = data.c_ex(physical_x, physical_y, t);
      Vector<Real> uh = phi * numerical(indices);

      // Gradients.
      Vector<Real> grad_x = data.dc_dx_ex(physical_x, physical_y, t);
      Vector<Real> grad_y = data.dc_dy_ex(physical_x, physical_y, t);
      Vector<Real> grad_u = grad_x + grad_y;
      Vector<Real> grad_uh = (gradx_phi + grady_phi) * numerical(indices);

      // Local L2 error.
      this->l2_errors[j] += dot(scaled, (u - uh) * (u - uh));

      // Local H1 error.
      this->h1_errors[j] +=
          dot(scaled, (grad_u - grad_uh) * (grad_u - grad_uh));
    }

    this->l2_errors[j] = std::sqrt(this->l2_errors[j]);
    this->h1_errors[j] = std::sqrt(this->h1_errors[j]);
  }

  // Data.
  this->degree = 0;

  for (const auto &element : mesh.elements)
    this->degree =
        (element.degree > this->degree) ? element.degree : this->degree;

  this->size = 0.0;
  this->elements = mesh.elements.size();

  for (const auto &element : mesh.elements)
    for (const auto &p : element.element.points)
      for (const auto &q : element.element.points)
        this->size =
            (distance(p, q) > this->size) ? distance(p, q) : this->size;
}
}