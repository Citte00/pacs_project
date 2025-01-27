/**
 * @file Errors.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2025-01-15
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <PacsHPDG.hpp>

namespace pacs {
    /**
     * @brief Compute Laplace equation errors.
     * 
     * @param data Laplace equation data struct.
     * @param mesh Mesh struct.
     * @param laplacian Laplace equation object.
     * @param numerical Numerical solution.
     */
void LaplaceError::computeErrors(const DataLaplace &data, const Mesh &mesh, const Laplace &laplacian, const Vector<Real> &numerical) {
#ifndef NVERBOSE
  std::cout << "Evaluating errors." << std::endl;
#endif

  // Matrices.
  Sparse<Real> mass = laplacian.M();
  Sparse<Real> dg_stiff = laplacian.DG();

  // Error vector.
  Vector<Real> u_modals = laplacian.modal(mesh, data.c_ex);
  Vector<Real> error = u_modals - numerical;

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

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangles[k], {nodes_x_2d, nodes_y_2d});

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

/**
 * @brief Compute Heat equation errors.
 *
 * @param data Heat equation data struct.
 * @param mesh Mesh struct.
 * @param heat Heat equation object.
 */
void HeatError::computeErrors(const DataHeat &data, const Mesh &mesh,
                                 const Heat &heat) {
#ifndef NVERBOSE
  std::cout << "Evaluating errors." << std::endl;
#endif

  // Matrices.
  Sparse<Real> mass = heat.M();
  Sparse<Real> dg_stiff = heat.DG();

  // Error vector.
  Vector<Real> u_modals = heat.modal(mesh, data.c_ex);

  Vector<Real> error = u_modals - heat.ch();

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

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangles[k], {nodes_x_2d, nodes_y_2d});

      // Weights scaling.
      Vector<Real> scaled = jacobian_det * weights_2d;

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      // Solutions.
      Vector<Real> u = data.c_ex(physical_x, physical_y, heat.t());
      Vector<Real> uh = phi * heat.ch()(indices);

      Vector<Real> grad_u = data.dc_dx_ex(physical_x, physical_y, heat.t()) +
                            data.dc_dy_ex(physical_x, physical_y, heat.t());
      Vector<Real> grad_uh = (gradx_phi + grady_phi) * heat.ch()(indices);

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

/**
 * @brief Compute Fisher equation errors.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh struct.
 * @param fisher Fisher equation object.
 */
void FisherError::computeErrors(const DataFKPP &data, const Mesh &mesh,
                              const Fisher &fisher) {
#ifndef NVERBOSE
  std::cout << "Evaluating errors." << std::endl;
#endif

  // Matrices.
  Sparse<Real> mass = fisher.M();
  Sparse<Real> dg_stiff = fisher.DG();

  // Error vector.
  Vector<Real> u_modals = fisher.modal(mesh, data.c_ex);

  Vector<Real> error = u_modals - fisher.ch();

  // DG Error.
  this->m_dg_error = std::sqrt(dot(error, dg_stiff * error));

  // L2 Error.
  this->m_l2_error = std::sqrt(dot(error, mass * error));

  // Energy error.
  this->m_energy += data.dt * this->m_dg_error * this->m_dg_error;

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

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangles[k], {nodes_x_2d, nodes_y_2d});

      // Weights scaling.
      Vector<Real> scaled = jacobian_det * weights_2d;

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      // Solutions.
      Vector<Real> u = data.c_ex(physical_x, physical_y, fisher.t());
      Vector<Real> uh = phi * fisher.ch()(indices);

      Vector<Real> grad_u = data.dc_dx_ex(physical_x, physical_y, fisher.t()) +
                            data.dc_dy_ex(physical_x, physical_y, fisher.t());
      Vector<Real> grad_uh = (gradx_phi + grady_phi) * fisher.ch()(indices);

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