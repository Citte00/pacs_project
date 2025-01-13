/**
 * @file fisher.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief
 * @date 2025-01-10
 *
 * @copyright Copyright (c) 2025
 *
 */
#include <PacsHPDG.hpp>

namespace pacs {

/**
 * @brief System matrices assembler.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh Struct.
 */
void pacs::Fisher::assembly(const DataFKPP &data, const Mesh &mesh) {

  // Resize matrices to new mesh.
  this->m_mass.reshape(mesh.dofs(), mesh.dofs());
  this->m_nl_mass.reshape(mesh.dofs(), mesh.dofs());
  this->m_stiff.reshape(mesh.dofs(), mesh.dofs());
  this->m_DG_stiff.reshape(mesh.dofs(), mesh.dofs());

  // Number of quadrature nodes.
  std::vector<std::size_t> nqn(mesh.elements.size(), 0);
  std::transform(mesh.elements.begin(), mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Degrees of freedom.
  std::size_t dofs = mesh.dofs();

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

  // Matrices.
  Sparse<Real> M{dofs, dofs};
  Sparse<Real> A{dofs, dofs};
  Sparse<Real> IA{dofs, dofs};
  Sparse<Real> SA{dofs, dofs};

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.emplace_back(0);

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

  // Volume integrals.

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // 2D Quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

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

    // Local matrices.
    Matrix<Real> local_M_prj{element_dofs, element_dofs};
    Matrix<Real> local_M{element_dofs, element_dofs};
    Matrix<Real> local_A{element_dofs, element_dofs};

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

      // Param initialization.
      Vector<Real> Dext = data.D_ext(physical_x, physical_y, 0.0);
      Vector<Real> alpha = data.alpha(physical_x, physical_y, 0.0);

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      // Some products.
      Matrix<Real> scaled_phi{phi};
      Matrix<Real> scaled_phi_alpha{phi};
      Matrix<Real> scaled_gradx{gradx_phi};
      Matrix<Real> scaled_grady{grady_phi};

      for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
        scaled_phi.column(l, scaled_phi.column(l) * scaled);
        scaled_phi_alpha.column(l,
                                (alpha * scaled_phi_alpha.column(l)) * scaled);
        scaled_gradx.column(l, (Dext * scaled_gradx.column(l)) * scaled);
        scaled_grady.column(l, (Dext * scaled_grady.column(l)) * scaled);
      }

      // Local matrix assembly.
      local_M_prj += scaled_phi.transpose() * phi;
      local_M += scaled_phi_alpha.transpose() * phi;
      local_A += scaled_gradx.transpose() * gradx_phi +
                 scaled_grady.transpose() * grady_phi;
    }

    // Global matrix assembly.
    this->m_mass.insert(indices, indices, local_M_prj);
    M.insert(indices, indices, local_M);
    A.insert(indices, indices, local_A);

    // Face integrals.

    // Local matrices.
    Matrix<Real> local_IA{element_dofs, element_dofs};
    Matrix<Real> local_SA{element_dofs, element_dofs};

    // Element's neighbours.
    std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

    // Local matrices for neighbours.
    std::vector<Matrix<Real>> local_IAN;
    std::vector<Matrix<Real>> local_SAN;

    // Penalties.
    Vector<Real> penalties = penalty(mesh, j, data.penalty_coeff);

    // Edges.
    std::vector<Segment> edges{polygon.edges()};

    // Loop over faces.
    for (std::size_t k = 0; k < element_neighbours.size(); ++k) {

      // Neighbour information.
      auto [edge, neighbour, n_edge] = element_neighbours[k];

      // 1D Quadrature nodes and weights.
      auto [nodes_1d, weights_1d] =
          (neighbour > 0) ? quadrature_1d(std::max(nqn[j], nqn[neighbour]))
                          : quadrature_1d(nqn[j]);

      // Edge geometry.
      Segment segment{edges[k]};

      // Edge's normal.
      Vector<Real> normal_vector{2};

      normal_vector[0] = segment[1][1] - segment[0][1];
      normal_vector[1] = segment[0][0] - segment[1][0];

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

      // Weights scaling.
      Vector<Real> scaled = std::abs(segment) * weights_1d;

      // Param initialization.
      Vector<Real> Dext = data.D_ext(physical_x, physical_y, 0.0);

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      // Local matrix assembly.
      Matrix<Real> scaled_gradx{gradx_phi};
      Matrix<Real> scaled_grady{grady_phi};
      Matrix<Real> scaled_phi{phi};

      for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
        scaled_gradx.column(l, (Dext * scaled_gradx.column(l)) * scaled);
        scaled_grady.column(l, (Dext * scaled_grady.column(l)) * scaled);
        scaled_phi.column(l, scaled_phi.column(l) * scaled);
      }

      Matrix<Real> scaled_grad =
          normal_vector[0] * scaled_gradx + normal_vector[1] * scaled_grady;

      if (neighbour == -1) { // Boundary edge.

        local_IA += scaled_grad.transpose() * phi;
        local_SA += (penalties[k] * scaled_phi).transpose() * phi;

        // Empty small matrices.
        local_IAN.emplace_back(Matrix<Real>{1, 1});
        local_SAN.emplace_back(Matrix<Real>{1, 1});

      } else {

        local_IA += 0.5 * scaled_grad.transpose() * phi;
        local_SA += (penalties[k] * scaled_phi).transpose() * phi;

        // Neighbour's basis function.
        Matrix<Real> n_phi =
            basis_2d(mesh, neighbour, {physical_x, physical_y})[0];

        // Neighbour's local matrix.
        local_IAN.emplace_back(-0.5 * scaled_grad.transpose() * n_phi);
        local_SAN.emplace_back(-(penalties[k] * scaled_phi).transpose() *
                               n_phi);
      }
    }

    IA.insert(indices, indices, local_IA);
    SA.insert(indices, indices, local_SA);

    // Neighbouring DG matrices assembly.
    for (std::size_t k = 0; k < element_neighbours.size(); ++k) {
      if (element_neighbours[k][1] == -1)
        continue;

      std::vector<std::size_t> n_indices;
      std::size_t n_index = element_neighbours[k][1];
      std::size_t n_dofs = mesh.elements[n_index].dofs(); // Neighbour's dofs.

      for (std::size_t h = 0; h < n_dofs; ++h)
        n_indices.emplace_back(starts[n_index] + h);

      IA.add(indices, n_indices, local_IAN[k]);
      SA.add(indices, n_indices, local_SAN[k]);
    }
  }

  // Matrices.
  this->m_nl_mass = M;
  this->m_DG_stiff = A + SA;
  this->m_stiff = this->m_DG_stiff - IA - IA.transpose();

  // Compression.
  this->m_mass.compress();
  this->m_nl_mass.compress();
  this->m_DG_stiff.compress();
  this->m_stiff.compress();
};

/**
 * @brief Forcing term assembler.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh struct.
 */
void Fisher::assemblyforce(const DataFKPP &data, const Mesh &mesh) {
#ifndef NVERBOSE
  std::cout << "Computing the forcing term." << std::endl;
#endif

  // Resize forcing vector to new mesh.
  m_forcing.resize(mesh.dofs());

  // Number of quadrature nodes.
  std::vector<std::size_t> nqn(mesh.elements.size(), 0);
  std::transform(mesh.elements.begin(), mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.emplace_back(0);

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

  // Volume integrals.

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {
    // 2D Local quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

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

    // Local forcing term.
    Vector<Real> local_f{element_dofs};

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

      // Local source evaluation.
      Vector<Real> Dext = data.D_ext(physical_x, physical_y, this->m_t);
      Vector<Real> alpha = data.alpha(physical_x, physical_y, this->m_t);
      Vector<Real> local_source =
          data.source_f(physical_x, physical_y, this->m_t, Dext, alpha);

      // Basis functions.
      auto phi = basis_2d(mesh, j, {physical_x, physical_y})[0];
      Matrix<Real> scaled_phi{phi};

      for (std::size_t l = 0; l < scaled_phi.columns; ++l)
        scaled_phi.column(l, scaled_phi.column(l) * scaled);

      // Local forcing term.
      local_f += scaled_phi.transpose() * local_source;
    }

    // Face integrals.

    // Element's neighbours.
    std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

    // Penalties.
    Vector<Real> penalties = penalty(mesh, j, data.penalty_coeff);

    // Edges.
    std::vector<Segment> edges{polygon.edges()};

    // Loop over faces.
    for (std::size_t k = 0; k < element_neighbours.size(); ++k) {
      // 1D Quadrature nodes and weights.
      auto [nodes_1d, weights_1d] = quadrature_1d(nqn[j]);

      // Neighbour information.
      auto [edge, neighbour, n_edge] = element_neighbours[k];

      // Only domain's border,
      if (neighbour != -1)
        continue;

      // Edge geometry.
      Segment segment{edges[k]}; // Mesh's edges to be fixed. [!]

      // Edge's normal. Check the order. [!]
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

      // Weights scaling.
      Vector<Real> scaled = std::abs(segment) * weights_1d;

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      // Param initialisation.
      Vector<Real> Dext = data.D_ext(physical_x, physical_y, this->m_t);

      // Local matrix assembly.
      Matrix<Real> scaled_gradx{gradx_phi};
      Matrix<Real> scaled_grady{grady_phi};
      Matrix<Real> scaled_phi{phi};

      for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
        scaled_gradx.column(l, (Dext * scaled_gradx.column(l)) * scaled);
        scaled_grady.column(l, (Dext * scaled_grady.column(l)) * scaled);
        scaled_phi.column(l, scaled_phi.column(l) * scaled);
      }

      Matrix<Real> scaled_grad =
          normal_vector[0] * scaled_gradx + normal_vector[1] * scaled_grady;

      // Boundary conditions.
      Vector<Real> boundary = data.DirBC(physical_x, physical_y, this->m_t);

      // Local forcing term.
      local_f -= scaled_grad.transpose() * boundary;
      local_f += penalties[k] * scaled_phi.transpose() * boundary;
    }

    // Global forcing term.
    this->m_forcing(indices, local_f);
  }
}

/**
 * @brief Non-linear matrix assmbler.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh Struct.
 */
Sparse<Real> Fisher::assemblyNL(const DataFKPP &data, const Mesh &mesh,
                                const Vector<Real> &c_star) {

#ifndef NVERBOSE
  std::cout << "Computing the Fisher-KPP matrices." << std::endl;
#endif

  // Number of quadrature nodes.
  std::vector<std::size_t> nqn(mesh.elements.size(), 0);
  std::transform(mesh.elements.begin(), mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Degrees of freedom.
  std::size_t dofs = mesh.dofs();

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

  // Matrices.
  Sparse<Real> M_star{dofs, dofs};

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.emplace_back(0);

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

  // Volume integrals.

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {
    // 2D Quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

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

    // Local matrices.
    Matrix<Real> local_M{element_dofs, element_dofs};

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

      // Param initialization.
      Vector<Real> alpha = data.alpha(physical_x, physical_y, 0.0);
      Vector<Real> c_star = c_star(indices);

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      // Some products.
      Matrix<Real> scaled_phi{phi};
      Vector<Real> c_star_loc = phi * c_star;

      for (std::size_t l = 0; l < scaled_phi.columns; ++l) {
        scaled_phi.column(l,
                          (alpha * scaled_phi.column(l)) * c_star_loc * scaled);
      }

      // Local matrix assembly.
      local_M += scaled_phi.transpose() * phi;
    }

    // Global matrix assembly.
    M_star.insert(indices, indices, local_M);
  }

  // Matrices.
  Sparse<Real> non_lin = M_star;

  // Compression.
  non_lin.compress();

  return non_lin;
}

/**
 * @brief Solver method for the Fisher-KPP equation.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh Struct.
 * @param TOL Tolerance.
 */
void Fisher::solver(const DataFKPP &data, const Mesh &mesh, const Real &TOL) {

  // Mass blocks.
  auto blocks = block_mass(mesh);

  // Assembling the constant component of the matrices.
  Sparse<Real> LHS = this->m_mass + data.dt * data.theta * this->m_stiff;
  Sparse<Real> RHS =
      this->m_mass - data.dt * (1.0 - data.theta) * this->m_stiff;

  // Assembling the dynamic componenet of the marices.
  if (data.theta == 0) {

    Sparse<Real> MN = assemblyNL(data, mesh, this->m_ch);
    RHS += data.dt * (this->m_nl_mass - MN);

  } else if (data.theta == 1) {

    Sparse<Real> MN = assemblyNL(data, mesh, this->m_ch);
    LHS -= data.dt * (this->m_nl_mass - MN);

  } else if (data.theta == 0.5) {

    Vector<Real> c_star = 1.5 * this->m_ch - 0.5 * this->m_ch_old;
    Sparse<Real> MN = assemblyNL(data, mesh, c_star);
    LHS -= 0.5 * data.dt * (this->m_nl_mass - MN);
    RHS += 0.5 * data.dt * (this->m_nl_mass - MN);
  }

  LHS.compress();
  RHS.compress();

  // Construction of the complete RHS for the theta method.
  Vector<Real> forcing_old = this->m_forcing;
  assemblyforce(data, mesh);
  Vector<Real> F = RHS * this->m_ch + data.dt * data.theta * this->m_forcing +
                   data.dt * (1 - data.theta) * forcing_old;

  // Solves using GMRES.
  m_ch = solve(LHS, F, blocks, GMRES, DBI, TOL);
};

/**
 * @brief Get coefficients of the source function.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh Struct.
 * @param t Time step.
 * @return Vector<Real>
 */
Vector<Real> Fisher::evaluateSource(const DataFKPP &data, const Mesh &mesh,
                                    const Real &t) const {
  // Number of quadrature nodes.
  std::vector<std::size_t> nqn(mesh.elements.size(), 0);
  std::transform(mesh.elements.begin(), mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Coefficients.
  Vector<Real> coefficients{mesh.dofs()};

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.emplace_back(0);

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // Quadrature nodes.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

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

    // Local coefficients.
    Vector<Real> local_coefficients{element_dofs};

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
      Matrix<Real> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];
      Matrix<Real> scaled_phi = phi;

      for (std::size_t l = 0; l < scaled_phi.columns; ++l)
        scaled_phi.column(l, scaled_phi.column(l) * scaled);

      // function solution.
      Vector<Real> Dext = data.D_ext(physical_x, physical_y, t);
      Vector<Real> alpha = data.alpha(physical_x, physical_y, t);
      Vector<Real> local_function =
          data.source_f(physical_x, physical_y, t, alpha, Dext);

      // Local coefficients.
      local_coefficients += scaled_phi.transpose() * local_function;
    }

    // Update.
    coefficients(indices, local_coefficients);
  }

  return coefficients;
}

/**
 * @brief Evaluate initial conditions.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh Struct.
 */
void Fisher::evaluateIC(const DataFKPP &data, const Mesh &mesh) {
  // Mass blocks.
  auto blocks = block_mass(mesh);

  // Initial condition.
  Vector<Real> c_h = evaluateCoeff(mesh, data.c_ex, 0.0);
  Vector<Real> c_hh = evaluateCoeff(mesh, data.c_ex, 0.0 - data.dt);

  // Projection for modal coordinates.
  m_ch = (norm(c_h) > TOLERANCE) ? solve(this->m_mass, c_h, blocks, DB)
                                 : Vector<Real>{mesh.dofs()};
  m_ch_old = (norm(c_hh) > TOLERANCE) ? solve(this->m_mass, c_hh, blocks, DB)
                                      : Vector<Real>{mesh.dofs()};
};

/*
prolong solution;

 */

/**
 * @brief Compute Fisher equation errors.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh struct.
 * @param fisher Fisher equation object.
 */
void FisherError::computeErrors(const DataFKPP &data, const Mesh &mesh,
                                const Fisher &fisher) {

  // Compute the L2 and DG errors.
  this->computeError(mesh, fisher, data.c_ex);

  // Compute the int_error and the energy error.
  this->m_int_error += data.dt * std::pow(this->m_DG_error, 2);
  this->m_energy = std::sqrt(std::pow(this->m_L2_error, 2) + this->m_int_error);

  // Resize forcing vector to new mesh.
  this->m_L2_errors.resize(mesh.dofs());
  this->m_H1_errors.resize(mesh.dofs());

  // Number of quadrature nodes.
  std::vector<std::size_t> nqn{mesh.elements.size(), 0};
  std::transform(mesh.elements.begin(), mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.emplace_back(0);

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

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

    // 2D Local quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

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
      Vector<Real> u = data.c_ex(physical_x, physical_y, fisher.t());
      Vector<Real> uh = phi * fisher.ch()(indices);

      Vector<Real> grad_u = data.dc_dx_ex(physical_x, physical_y, fisher.t()) +
                            data.dc_dy_ex(physical_x, physical_y, fisher.t());
      Vector<Real> grad_uh = (gradx_phi + grady_phi) * fisher.ch()(indices);

      // Local L2 error.
      this->m_L2_errors[j] += dot(scaled, (u - uh) * (u - uh));

      // Local H1 error.
      this->m_H1_errors[j] +=
          dot(scaled, (grad_u - grad_uh) * (grad_u - grad_uh));
    }

    this->m_L2_errors[j] = std::sqrt(this->m_L2_errors[j]);
    this->m_H1_errors[j] = std::sqrt(this->m_H1_errors[j]);
  }
};

/**
 * @brief Compute Fisher-KPP equation error estimator.
 *
 * @param data Fisher-KPP equation data struct.
 * @param mesh Mesh struct.
 * @param fisher Fisher-KPP equation object.
 */
void FisherEstimator::computeEstimate(const DataFKPP &data, const Mesh &mesh,
                                      const Fisher &fisher) {
#ifndef NVERBOSE
  std::cout << "Evaluating estimates." << std::endl;
#endif

  // Number of quadrature nodes.
  std::vector<std::size_t> nqn(mesh.elements.size(), 0);
  std::transform(mesh.elements.begin(), mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.emplace_back(0);

  for (std::size_t j = 1; j < mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

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

  // Mass blocks.
  auto blocks = fisher.block_mass(mesh);

  // Coefficients.
  Vector<Real> f_coeff = fisher.evaluateSource(data, mesh, fisher.t());
  Vector<Real> g_coeff = fisher.evaluateCoeff(mesh, data.DirBC, fisher.t());

  Vector<Real> f_modals = (norm(f_coeff) > TOLERANCE)
                              ? solve(fisher.matrices()[0], f_coeff, blocks, DB)
                              : Vector<Real>{mesh.dofs()};
  Vector<Real> g_modals = (norm(g_coeff) > TOLERANCE)
                              ? solve(fisher.matrices()[0], g_coeff, blocks, DB)
                              : Vector<Real>{mesh.dofs()};

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // 2D Quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

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
      Matrix<Real> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];
      Matrix<Real> lap_phi = lap_basis_2d(mesh, j, {physical_x, physical_y});

      // Local numerical time derivative.
      Vector<Real> local_uh = phi * fisher.ch()(indices);
      Vector<Real> local_uh_old = phi * fisher.ch_old()(indices);
      Vector<Real> partial_uh = (local_uh - local_uh_old) / data.dt;

      // Local numerical laplacian.
      Vector<Real> lap_uh = lap_phi * fisher.ch()(indices);

      // Local exact source.
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y, fisher.t());
      Vector<Real> alpha = data.alpha(physical_x, physical_y, fisher.t());
      Vector<Real> f =
          data.source_f(physical_x, physical_y, fisher.t(), D_ext, alpha);

      // Local source approximation.
      Vector<Real> f_bar = phi * f_modals(indices);

      // Local non-linear term approx.
      Vector<Real> non_lin_uh = alpha * local_uh * (1.0 - local_uh);

      // Local estimator, R_{K, E}^2.
      this->m_estimates[j] +=
          sizes[j] * sizes[j] *
          dot(scaled, (f_bar + lap_uh - partial_uh + non_lin_uh) *
                          (f_bar + lap_uh - partial_uh + non_lin_uh));

      // Local data oscillation, O_{K, E}^2.
      this->m_estimates[j] +=
          sizes[j] * sizes[j] * dot(scaled, (f - f_bar) * (f - f_bar));
    }

    // Element's neighbours.
    std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

    // Penalties.
    Vector<Real> penalties = penalty(mesh, j, data.penalty_coeff);

    // Edges.
    std::vector<Segment> edges{polygon.edges()};

    // Loop over faces.
    for (std::size_t k = 0; k < element_neighbours.size(); ++k) {
      // Neighbour information.
      auto [edge, neighbour, n_edge] = element_neighbours[k];

      // 1D Quadrature nodes and weights.
      auto [nodes_1d, weights_1d] =
          (neighbour > 0) ? quadrature_1d(std::max(nqn[j], nqn[neighbour]))
                          : quadrature_1d(nqn[j]);

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
      edge_vector /= norm(edge_vector);

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

      // Weights scaling.
      Vector<Real> scaled = std::abs(segment) * weights_1d;

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      Matrix<Real> scaled_phi{phi};

      for (std::size_t l = 0; l < scaled_phi.columns; ++l)
        scaled_phi.column(l, scaled_phi.column(l) * scaled);

      // Local numerical solution and gradients.
      Vector<Real> uh = phi * fisher.ch()(indices);

      Matrix<Real> grad =
          normal_vector[0] * gradx_phi + normal_vector[1] * grady_phi;
      Vector<Real> grad_uh = grad * fisher.ch()(indices);

      Matrix<Real> grad_t =
          edge_vector[0] * gradx_phi + edge_vector[1] * grady_phi;
      Vector<Real> grad_uh_t = grad_t * fisher.ch()(indices);

      if (neighbour == -1) { // Boundary edge.

        // Local exact Dirichlet and gradient.
        Vector<Real> g = data.DirBC(physical_x, physical_y, fisher.t());
        Vector<Real> grad_g_t =
            edge_vector[0] * data.dc_dx_ex(physical_x, physical_y, fisher.t()) +
            edge_vector[1] * data.dc_dy_ex(physical_x, physical_y, fisher.t());

        // Approximate Dirichlet and gradient.
        Vector<Real> g_bar = phi * g_modals(indices);
        Vector<Real> grad_g_t_bar = grad_t * g_modals(indices);

        // Local estimator, R_{K, J}^2.
        this->m_estimates[j] +=
            penalties[k] * dot(scaled, (uh - g_bar) * (uh - g_bar));

        // Local estimator, R_{K, T}^2.
        this->m_estimates[j] +=
            sizes[j] * dot(scaled, (grad_uh_t - grad_g_t_bar) *
                                       (grad_uh_t - grad_g_t_bar));

        // Local data oscillation, O_{K, J}^2.
        this->m_estimates[j] +=
            penalties[k] * dot(scaled, (g - g_bar) * (g - g_bar));

        // Local data oscillation, O_{K, T}^2.
        this->m_estimates[j] +=
            sizes[j] *
            dot(scaled, (grad_g_t - grad_g_t_bar) * (grad_g_t - grad_g_t_bar));

      } else {
        // Neighbour's basis function.
        auto [n_phi, n_gradx_phi, n_grady_phi] =
            basis_2d(mesh, neighbour, {physical_x, physical_y});

        std::vector<std::size_t> n_indices;
        std::size_t n_index = element_neighbours[k][1];
        std::size_t n_dofs = mesh.elements[n_index].dofs(); // Neighbour's dofs.

        for (std::size_t h = 0; h < n_dofs; ++h)
          n_indices.emplace_back(starts[n_index] + h);

        // Neighbour's numerical solution and gradients.
        Vector<Real> n_uh = n_phi * fisher.ch()(n_indices);

        Matrix<Real> n_grad =
            normal_vector[0] * n_gradx_phi + normal_vector[1] * n_grady_phi;
        Vector<Real> n_grad_uh = n_grad * fisher.ch()(n_indices);

        Matrix<Real> n_grad_t =
            edge_vector[0] * n_gradx_phi + edge_vector[1] * n_grady_phi;
        Vector<Real> n_grad_uh_t = n_grad_t * fisher.ch()(n_indices);

        // Local estimator, R_{K, J}^2.
        this->m_estimates[j] +=
            penalties[k] * dot(scaled, (uh - n_uh) * (uh - n_uh));

        // Local estimator, R_{K, N}^2.
        this->m_estimates[j] +=
            sizes[j] *
            dot(scaled, (grad_uh - n_grad_uh) * (grad_uh - n_grad_uh));

        // Local estimator, R_{K, T}^2.
        this->m_estimates[j] +=
            sizes[j] *
            dot(scaled, (grad_uh_t - n_grad_uh_t) * (grad_uh_t - n_grad_uh_t));
      }
    }

    this->m_estimate += this->m_estimates[j];
    this->m_estimates[j] = std::sqrt(this->m_estimates[j]);

    // Degrees.
    Vector<Real> degrees{element_dofs};
    std::size_t counter = 0;

    for (std::size_t i = 0; i < mesh.elements[j].degree + 1; ++i)
      for (std::size_t k = 0; k < mesh.elements[j].degree + 1 - i; ++k) {
        degrees[counter] = static_cast<Real>(i + k);
        ++counter;
      }

    // Coefficients.
    Vector<Real> coefficients = fisher.ch()(indices);

    for (auto &coefficient : coefficients.elements)
      coefficient = std::log(std::abs(coefficient));

    // Fit.
    Vector<Real> fit = polyfit(degrees, coefficients, 1);
    this->m_fits[j] = -fit[1];
  }

  this->m_estimate = std::sqrt(this->m_estimate);
};
} // namespace pacs
