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
 * @brief Matrices assembly for Fisher-KPP equation.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh struct.
 */
void Fisher::assembly(const DataFKPP &data, const Mesh &mesh) {
#ifndef NVERBOSE
  std::cout << "Computing the Fisher-KPP matrices." << std::endl;
#endif

  // Number of elements.
  std::size_t num_elem = mesh.elements.size();

  // Starting indices.
  std::vector<std::size_t> starts(num_elem);
  starts[0] = 0;

  // Quadrature nodes.
  std::vector<std::size_t> nqn(num_elem);
  nqn[0] = 2 * mesh.elements[0].degree + 1;

  for (std::size_t j = 1; j < num_elem; ++j) {
    starts[j] = starts[j - 1] + mesh.elements[j - 1].dofs();
    nqn[j] = 2 * mesh.elements[j].degree + 1;
  }

  // Degrees of freedom.
  std::size_t dofs = mesh.dofs();

  // Matrices.
  Sparse<Real> M_prj{dofs, dofs};
  Sparse<Real> M{dofs, dofs};
  Sparse<Real> A{dofs, dofs};
  Sparse<Real> IA{dofs, dofs};
  Sparse<Real> SA{dofs, dofs};

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

  // Volume integrals.
  // Loop over the elements.
#pragma omp parallel for schedule(dynamic)
  for (std::size_t j = 0; j < num_elem; ++j) {

    // 2D Quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

    // Global matrix indices.
    std::vector<std::size_t> indices(mesh.elements[j].dofs());

    for (std::size_t k = 0; k < mesh.elements[j].dofs(); ++k)
      indices[k] = starts[j] + k;

    // Polygon.
    Polygon polygon = mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Local matrices.
    Matrix<Real> local_M_prj{indices.size(), indices.size()};
    Matrix<Real> local_M{indices.size(), indices.size()};
    Matrix<Real> local_A{indices.size(), indices.size()};

    // Loop over the sub-triangulation.
    for (const auto &triangle : triangles) {

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

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
#pragma omp critical
    {
      M_prj.insert(indices, indices, local_M_prj);
      M.insert(indices, indices, local_M);
      A.insert(indices, indices, local_A);
    }
    // Face integrals.

    // Local matrices.
    Matrix<Real> local_IA{indices.size(), indices.size()};
    Matrix<Real> local_SA{indices.size(), indices.size()};

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

      auto [normal_vector, edge_vector, physical_x, physical_y] =
          faces_physical_points(edges[k], nodes_1d);

      // Weights scaling.
      Vector<Real> scaled = std::abs(edges[k]) * weights_1d;

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

#pragma omp critical
    {
      IA.insert(indices, indices, local_IA);
      SA.insert(indices, indices, local_SA);
    }

    // Neighbouring DG matrices assembly.
    for (std::size_t k = 0; k < element_neighbours.size(); ++k) {
      if (element_neighbours[k][1] == -1)
        continue;

      std::vector<std::size_t> n_indices;
      std::size_t n_index = element_neighbours[k][1];
      std::size_t n_dofs = mesh.elements[n_index].dofs(); // Neighbour's dofs.

      for (std::size_t h = 0; h < n_dofs; ++h)
        n_indices.emplace_back(starts[n_index] + h);

#pragma omp critical
      {
        IA.add(indices, n_indices, local_IAN[k]);
        SA.add(indices, n_indices, local_SAN[k]);
      }
    }
  }

  // Matrices.
  this->m_mass = std::move(M_prj);
  this->m_nl_mass = std::move(M);
  this->m_dg_stiff = std::move(A) + SA;
  this->m_stiff = std::move(this->m_dg_stiff) - IA - IA.transpose();

  // Compression.
  this->m_mass.compress();
  this->m_nl_mass.compress();
  this->m_dg_stiff.compress();
  this->m_stiff.compress();
};

/**
 * @brief Assembly of the non-linear matrix for Fisher-KPP equation.
 *
 * @param data Data struct.
 * @param mesh Mesh struct.
 * @param uh Numerical solution.
 * @return Sparse<Real>
 */
Sparse<Real> Fisher::assembly_nl(const DataFKPP &data, const Mesh &mesh,
                                 const Vector<Real> &uh) {
#ifndef NVERBOSE
  std::cout << "Computing the Fisher-KPP matrices." << std::endl;
#endif

  // Number of elements.
  std::size_t num_elem = mesh.elements.size();

  // Starting indices.
  std::vector<std::size_t> starts(num_elem);
  starts[0] = 0;

  // Quadrature nodes.
  std::vector<std::size_t> nqn(num_elem);
  nqn[0] = 2 * mesh.elements[0].degree + 1;

  for (std::size_t j = 1; j < num_elem; ++j) {
    starts[j] = starts[j - 1] + mesh.elements[j - 1].dofs();
    nqn[j] = 2 * mesh.elements[j].degree + 1;
  }

  // Degrees of freedom.
  std::size_t dofs = mesh.dofs();

  // Matrices.
  Sparse<Real> M_star{dofs, dofs};

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

// Loop over the elements.
#pragma omp parallel for schedule(dynamic)
  for (std::size_t j = 0; j < num_elem; ++j) {

    // 2D Quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

    // Global matrix indices.
    std::vector<std::size_t> indices(mesh.elements[j].dofs());

    for (std::size_t k = 0; k < mesh.elements[j].dofs(); ++k)
      indices[k] = starts[j] + k;

    // Polygon.
    Polygon polygon = mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Local matrices.
    Matrix<Real> local_M{indices.size(), indices.size()};

    // Loop over the sub-triangulation.
    for (const auto &triangle : triangles) {

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

      // Weights scaling.
      Vector<Real> scaled = jacobian_det * weights_2d;

      // Param initialization.
      Vector<Real> alpha = data.alpha(physical_x, physical_y, 0.0);
      Vector<Real> c_star = uh(indices);

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
#pragma omp critical
    M_star.insert(indices, indices, local_M);
  }

  // Matrices.
  Sparse<Real> non_lin = std::move(M_star);

  // Compression.
  non_lin.compress();

  return non_lin;
}

/**
 * @brief Forcing term assembly.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh struct.
 */
void Fisher::assembly_force(const DataFKPP &data, const Mesh &mesh) {
#ifndef NVERBOSE
  std::cout << "Computing the forcing term." << std::endl;
#endif

  // Number of elements.
  std::size_t num_elem = mesh.elements.size();

  // Starting indices.
  std::vector<std::size_t> starts(num_elem);
  starts[0] = 0;

  // Quadrature nodes.
  std::vector<std::size_t> nqn(num_elem);
  nqn[0] = 2 * mesh.elements[0].degree + 1;

  for (std::size_t j = 1; j < num_elem; ++j) {
    starts[j] = starts[j - 1] + mesh.elements[j - 1].dofs();
    nqn[j] = 2 * mesh.elements[j].degree + 1;
  }

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

  // Loop over the elements.
#pragma omp parallel for schedule(dynamic)
  for (std::size_t j = 0; j < num_elem; ++j) {

    // 2D Local quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

    // Global matrix indices.
    std::vector<std::size_t> indices(mesh.elements[j].dofs());

    for (std::size_t k = 0; k < mesh.elements[j].dofs(); ++k)
      indices[k] = starts[j] + k;

    // Polygon.
    Polygon polygon = mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Local forcing term.
    Vector<Real> local_f{indices.size()};

    // Loop over the sub-triangulation.
    for (const auto &triangle : triangles) {

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

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

      auto [normal_vector, edge_vector, physical_x, physical_y] =
          faces_physical_points(edges[k], nodes_1d);

      // Weights scaling.
      Vector<Real> scaled = std::abs(edges[k]) * weights_1d;

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
#pragma omp critical
    this->m_forcing(indices, local_f);
  }
};

/**
 * @brief Fisher-KPP equation solver.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh struct.
 * @param ch_oold Solution vector at previous time step.
 * @param forcing_old Forcing term.
 * @param TOL Tolerance.
 */
Vector<Real> Fisher::solver(const DataFKPP &data, const Mesh &mesh,
                            const Vector<Real> &ch_oold,
                            const Vector<Real> &forcing_old, const Real &TOL) {
#ifndef NVERBOSE
  std::cout << "Solving the algebraic system." << std::endl;
#endif

  // Mass blocks.
  auto blocks = block_mass(mesh);

  // Assembling the constant component of the matrices.
  Sparse<Real> LHS = this->m_mass + data.dt * data.theta * this->m_stiff;
  Sparse<Real> RHS =
      this->m_mass - data.dt * (1.0 - data.theta) * this->m_stiff;

  // Assembling the dynamic componenet of the marices.
  if (data.theta == 0) {

    Sparse<Real> MN = assembly_nl(data, mesh, this->m_ch_old);
    RHS += data.dt * (this->m_mass - MN);

  } else if (data.theta == 1) {

    Sparse<Real> MN = assembly_nl(data, mesh, this->m_ch_old);
    LHS -= data.dt * (this->m_mass - MN);

  } else if (data.theta == 0.5) {

    Vector<Real> c_star = 1.5 * this->m_ch_old - 0.5 * ch_oold;
    Sparse<Real> MN = assembly_nl(data, mesh, c_star);
    LHS -= 0.5 * data.dt * (this->m_mass - MN);
    RHS += 0.5 * data.dt * (this->m_mass - MN);
  }

  LHS.compress();
  RHS.compress();

  // Construction of the complete RHS for the theta method.
  Vector<Real> F = RHS * this->m_ch_old +
                   data.dt * data.theta * this->m_forcing +
                   data.dt * (1.0 - data.theta) * forcing_old;

  // Solves using GMRES.
  return solve(LHS, F, blocks, GMRES, DBI, TOL);
}

/**
 * @brief Evaluate the modal coefficient of the Fisher source term.
 *
 * @param data Fisher equation data struct.
 * @param mesh Mesh struct.
 * @return Vector<Real>
 */
Vector<Real> Fisher::modal_source(const DataFKPP &data,
                                  const Mesh &mesh) const {
#ifndef NVERBOSE
  std::cout << "Retrieving modal coeffficient of the source term." << std::endl;
#endif

  // Number of elements.
  std::size_t num_elem = mesh.elements.size();

  // Starting indices.
  std::vector<std::size_t> starts(num_elem);
  starts[0] = 0;

  // Quadrature nodes.
  std::vector<std::size_t> nqn(num_elem);
  nqn[0] = 2 * mesh.elements[0].degree + 1;

  for (std::size_t j = 1; j < num_elem; ++j) {
    starts[j] = starts[j - 1] + mesh.elements[j - 1].dofs();
    nqn[j] = 2 * mesh.elements[j].degree + 1;
  }

  // Coefficients.
  Vector<Real> coefficients{mesh.dofs()};

// Loop over the elements.
#pragma omp parallel for schedule(dynamic)
  for (std::size_t j = 0; j < num_elem; ++j) {

    // Quadrature nodes.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

    // Global matrix indices.
    std::vector<std::size_t> indices(mesh.elements[j].dofs());

    for (std::size_t k = 0; k < mesh.elements[j].dofs(); ++k)
      indices[k] = starts[j] + k;

    // Polygon.
    Polygon polygon = mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Local coefficients.
    Vector<Real> local_coefficients{indices.size()};

    // Loop over the sub-triangulation.
    for (const auto &triangle : triangles) {

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

      // Weights scaling.
      Vector<Real> scaled = jacobian_det * weights_2d;

      // Basis functions.
      Matrix<Real> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];
      Matrix<Real> scaled_phi = phi;

      for (std::size_t l = 0; l < scaled_phi.columns; ++l)
        scaled_phi.column(l, scaled_phi.column(l) * scaled);

      // function solution.
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y, m_t);
      Vector<Real> alpha = data.alpha(physical_x, physical_y, m_t);
      Vector<Real> local_function =
          data.source_f(physical_x, physical_y, m_t, D_ext, alpha);

      // Local coefficients.
      local_coefficients += scaled_phi.transpose() * local_function;
    }

// Update.
#pragma omp critical
    coefficients(indices, local_coefficients);
  }

  auto blocks = block_mass(mesh);

  return (norm(coefficients) > TOLERANCE)
             ? solve(this->m_mass, coefficients, blocks, DB)
             : Vector<Real>{mesh.dofs()};
};

} // namespace pacs