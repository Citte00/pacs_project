/**
 * @file heat.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Contains the implementation of all Heat object methods.
 * @date 2024-12-16
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <PacsHPDG.hpp>

namespace pacs {

/**
 * @brief Update Heat matrices and forcing term in adaptive framework.
 *
 * @param data Heat equation data structure.
 * @param mesh Mesh structure.
 */
void Heat::update(const DataHeat &data, const Mesh &mesh) {
  // Dofs.
  std::size_t dofs = mesh.dofs();

  // Update matrices.
  m_mass.reshape(dofs, dofs);
  m_stiff.reshape(dofs, dofs);
  m_dg_stiff.reshape(dofs, dofs);

  assembly(data, mesh);

  // Update forcing term.
  m_forcing.resize(dofs);

  assembly_force(data, mesh);
};

/**
 * @brief Assembly of the heat equation matrices.
 *
 * @param data Heat equation data structure.
 * @param mesh Mesh structure.
 * @return std::array<Sparse<Real>, 3>
 */
void Heat::assembly(const DataHeat &data, const Mesh &mesh) {
#ifndef NVERBOSE
  std::cout << "Computing the Heat equation matrices." << std::endl;
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
  Sparse<Real> M{dofs, dofs};
  Sparse<Real> A{dofs, dofs};
  Sparse<Real> IA{dofs, dofs};
  Sparse<Real> SA{dofs, dofs};

  // Neighbours.
  const std::vector<std::vector<std::array<int, 3>>> &neighbours =
      mesh.neighbours;

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
    const Polygon &polygon = mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Local matrices.
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
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y, 0.0);

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      // Some products.
      Matrix<Real> scaled_phi{phi};
      Matrix<Real> scaled_gradx{gradx_phi};
      Matrix<Real> scaled_grady{grady_phi};

      for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
        scaled_phi.column(l, scaled_phi.column(l) * scaled);
        scaled_gradx.column(l, (D_ext * scaled_gradx.column(l)) * scaled);
        scaled_grady.column(l, (D_ext * scaled_grady.column(l)) * scaled);
      }

      // Local matrix assembly.
      local_M += scaled_phi.transpose() * phi;
      local_A += scaled_gradx.transpose() * gradx_phi +
                 scaled_grady.transpose() * grady_phi;
    }

// Global matrix assembly.
#pragma omp critical
    {
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
    local_IAN.reserve(element_neighbours.size());
    std::vector<Matrix<Real>> local_SAN;
    local_SAN.reserve(element_neighbours.size());

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
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y, 0.0);

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      // Local matrix assembly.
      Matrix<Real> scaled_gradx{gradx_phi};
      Matrix<Real> scaled_grady{grady_phi};
      Matrix<Real> scaled_phi{phi};

      for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
        scaled_gradx.column(l, (D_ext * scaled_gradx.column(l)) * scaled);
        scaled_grady.column(l, (D_ext * scaled_grady.column(l)) * scaled);
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

      n_indices.reserve(n_dofs);
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
  this->m_mass = std::move(M);
  this->m_dg_stiff = std::move(A) + SA;
  this->m_stiff = this->m_dg_stiff - IA - IA.transpose();

  // Compression.
  this->m_mass.compress();
  this->m_dg_stiff.compress();
  this->m_stiff.compress();
};

/**
 * @brief Assembly of the forcing term for the Heat eqaution.
 *
 * @param data Heat equation data structure.
 * @param mesh Mesh structure.
 * @return Vector<Real>
 */
void Heat::assembly_force(const DataHeat &data, const Mesh &mesh) {
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

  // Volume integrals.
  Vector<Real> forcing{mesh.dofs()};

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
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y, m_t);
      Vector<Real> local_source =
          data.source_f(physical_x, physical_y, m_t, D_ext);

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
      // Neighbour information.
      auto [edge, neighbour, n_edge] = element_neighbours[k];

      // 1D Quadrature nodes and weights.
      auto [nodes_1d, weights_1d] = quadrature_1d(nqn[j]);

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

      // Local matrix assembly.
      Matrix<Real> scaled_gradx{gradx_phi};
      Matrix<Real> scaled_grady{grady_phi};
      Matrix<Real> scaled_phi{phi};

      // Boundary conditions.
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y, m_t);
      Vector<Real> boundary = data.DirBC(physical_x, physical_y, m_t);

      for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
        scaled_gradx.column(l, (D_ext * scaled_gradx.column(l)) * scaled);
        scaled_grady.column(l, (D_ext * scaled_grady.column(l)) * scaled);
        scaled_phi.column(l, scaled_phi.column(l) * scaled);
      }

      Matrix<Real> scaled_grad =
          normal_vector[0] * scaled_gradx + normal_vector[1] * scaled_grady;

      // Local forcing term.
      local_f -= scaled_grad.transpose() * boundary;
      local_f += penalties[k] * scaled_phi.transpose() * boundary;
    }

// Global forcing term.
#pragma omp critical
    forcing(indices, local_f);
  }
  this->m_forcing = std::move(forcing);
};

/**
 * @brief Solver of the heat equation matricial system.
 *
 * @param data Data struct.
 * @param mesh Mesh struct.
 * @param Matrices System matrices, [mass, stiff, dg_stiff].
 * @param c_old Soltuion at previous time step.
 * @param forcing Array of forcing terms.
 * @param TOL Tolerance.
 * @return Vector<Real>
 */
Vector<Real> Heat::solver(const DataHeat &data, const Mesh &mesh,
                          const Vector<Real> &c_old,
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

  LHS.compress();
  RHS.compress();

  // Construction of the complete RHS for the theta method.
  Vector<Real> F = RHS * c_old + data.dt * data.theta * this->m_forcing +
                   data.dt * (1 - data.theta) * forcing_old;

  // Solves using BICGSTAB.
  return solve(LHS, F, blocks, BICGSTAB, DBI, TOL);
};

/**
 * @brief Evaluate the modal coefficients of a function.
 *
 * @param mesh Mesh.
 * @param function Function to evaluate.
 * @return Vector<Real>
 */
Vector<Real> Heat::modal(const Mesh &mesh, const TriFunctor &function) const {
#ifndef NVERBOSE
  std::cout << "Retrieving modal coefficients." << std::endl;
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
      Vector<Real> local_function = function(physical_x, physical_y, m_t);

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

/**
 * @brief Evaluate the modal coefficient of the Heat source term.
 *
 * @param data Heat equation data struct.
 * @param mesh Mesh struct.
 * @return Vector<Real>
 */
Vector<Real> Heat::modal_source(const DataHeat &data, const Mesh &mesh) const {
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
      Vector<Real> local_function =
          data.source_f(physical_x, physical_y, m_t, D_ext);

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

/**
 * @brief Compute the indices related to bases for each degree.
 *
 * @param mesh Mesh struct.
 * @return Matrix<int>
 */
Matrix<int> Heat::transition(const std::size_t &degree) const {
#ifndef NVERBOSE
  std::cout << "Computing transition matrix." << std::endl;
#endif

  // Computation of the number of required basis
  std::size_t N_max = degree + 1;
  std::size_t length = (N_max * (N_max + 1)) / 2;

  // First indexes vector.
  std::vector<int> idx1;
  idx1.reserve(length);
  for (size_t i = 0; i < N_max; i++)
    idx1.insert(idx1.end(), N_max - i, i);

  // Second indexes vector.
  std::vector<int> idx2;
  idx2.reserve(length);
  for (size_t i = 0; i < N_max; i++)
    for (size_t j = i; j < N_max; j++)
      idx2.emplace_back(j);

  // Concatenate the two vector in a matrix.
  Matrix<int> idxglob{2, idx1.size()};
  idxglob.row(0, Vector<int>(length, idx1));
  idxglob.row(1, Vector<int>(length, idx2));

  Matrix<int> transition{N_max - 1, length};

  for (size_t i = 1; i < N_max; i++) {
    std::size_t N = i + 1;
    std::size_t len = N * (N + 1) / 2;

    // First indexes vector.
    std::vector<int> idx1_loc;
    idx1_loc.reserve(len);
    for (size_t i = 0; i < N; i++)
      idx1_loc.insert(idx1_loc.end(), N - i, i);

    // Second indexes vector.
    std::vector<int> idx2_loc;
    idx2_loc.reserve(len);
    for (size_t i = 0; i < N; i++)
      for (size_t j = i; j < N; j++)
        idx2_loc.emplace_back(j);

    // Concatenate the two vector in a matrix.
    Matrix<int> idxgloc{2, idx1_loc.size()};
    idxgloc.row(0, Vector<int>(len, idx1_loc));
    idxgloc.row(1, Vector<int>(len, idx2_loc));

    std::vector<int> idx(length, 0);
    for (std::size_t j = 0; j < len; j++) {

      std::vector<int> col;
      col.reserve(length);

      for (std::size_t k = 0; k < length; k++)
        col.emplace_back((idxglob(0, k) == idxgloc(0, j)) &&
                         (idxglob(1, k) == idxgloc(1, j)));

      std::transform(idx.begin(), idx.end(), col.begin(), idx.begin(),
                     std::plus<>());
    }

    transition.row(i - 1, Vector<int>(length, idx));
  }
  return transition;
};

/**
 * @brief Adapt solution vector to refined polynomial degrees.
 *
 * @param new_mesh New mesh struct.
 * @param old_mesh Old mesh struct.
 * @param mask_p Vector identifying element to refine through degree.
 */
Vector<Real> Heat::prolong_solution_p(const Mesh &new_mesh,
                                      const Mesh &old_mesh,
                                      const Vector<Real> &ch,
                                      const Mask &mask_p) const {
#ifndef NVERBOSE
  std::cout << "Prolong solution to new dofs." << std::endl;
#endif

  // Number of elements.
  std::size_t old_elem = old_mesh.elements.size();
  std::size_t new_elem = new_mesh.elements.size();

  // Number of quadrature nodes.
  std::vector<std::size_t> nqn(new_elem);
  nqn[0] = 2 * new_mesh.elements[0].degree + 1;

  // Resize new solution vector.
  Vector<Real> new_ch(new_mesh.dofs());

  // Starting indices.
  std::vector<std::size_t> old_starts(old_elem);
  old_starts[0] = 0;

  for (std::size_t j = 1; j < old_elem; ++j)
    old_starts[j] = old_starts[j - 1] + old_mesh.elements[j - 1].dofs();

  std::vector<std::size_t> new_starts(new_elem);
  new_starts[0] = 0;

  for (std::size_t j = 1; j < new_elem; ++j) {
    new_starts[j] = new_starts[j - 1] + new_mesh.elements[j - 1].dofs();
    nqn[j] = 2 * new_mesh.elements[j].degree + 1;
  }

// Loop over elements.
#pragma omp parallel for schedule(dynamic)
  for (std::size_t j = 0; j < new_elem; j++) {

    // 2D quadrature weights and nodes.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

    // Local dofs.
    std::size_t old_element_dofs = old_mesh.elements[j].dofs();
    std::size_t new_element_dofs = new_mesh.elements[j].dofs();

    // Global matrix indices.
    std::vector<std::size_t> old_indices(old_element_dofs);
    for (std::size_t k = 0; k < old_element_dofs; ++k)
      old_indices[k] = old_starts[j] + k;

    std::vector<std::size_t> new_indices(new_element_dofs);
    for (std::size_t k = 0; k < new_element_dofs; ++k)
      new_indices[k] = new_starts[j] + k;

    // Check to refine.
    if (!mask_p[j]) {
      new_ch(new_indices, ch(old_indices));
      continue;
    }

    // Polygon.
    Polygon polygon = new_mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Local mass matrix.
    Matrix<Real> local_mass{new_indices.size(), new_indices.size()};

    // Local coefficients.
    Vector<Real> local_coefficients{new_indices.size()};

    // Loop over the sub-triangulation.
    for (const auto &triangle : triangles) {

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

      // Weights scaling.
      Vector<Real> scaled = jacobian_det * weights_2d;

      // Basis functions.
      Matrix<Real> phi = basis_2d(new_mesh, j, {physical_x, physical_y})[0];
      Matrix<Real> phi_old = basis_2d(old_mesh, j, {physical_x, physical_y})[0];
      Matrix<Real> scaled_phi{phi};

      for (std::size_t l = 0; l < scaled_phi.columns; ++l)
        scaled_phi.column(l, scaled_phi.column(l) * scaled);

      // Computing the transition to identify basis indices.
      Matrix<int> transitions = transition(new_mesh.elements[j].degree);

      // I have to identify the indices that corresponds to the same basis.
      std::vector<std::size_t> indexes;
      indexes.reserve(transitions.columns);
      for (size_t i = 0; i < transitions.columns; i++)
        if (transitions(old_mesh.elements[j].degree - 1, i) == 1)
          indexes.emplace_back(i);

      // Local coefficients.
      local_mass += scaled_phi.transpose() * phi;
      local_coefficients += scaled_phi.transpose() * (phi_old * ch(old_indices));
    }

    // Update the solution.
    local_coefficients = solve(local_mass, local_coefficients);
#pragma omp critical
    new_ch(new_indices, local_coefficients);
  }
  return new_ch;
};

/**
 * @brief Adapt solution vector to refined mesh.
 *
 * @param new_mesh New mesh struct.
 * @param old_mesh Old mesh struct.
 * @param mask_h Vector identifying element to refine through size.
 */
Vector<Real> Heat::prolong_solution_h(const Mesh &new_mesh,
                                      const Mesh &old_mesh,
                                      const Vector<Real> &ch,
                                      const Mask &mask_h) const {
#ifndef NVERBOSE
  std::cout << "Prolong solution to new elements." << std::endl;
#endif

  // Number of elements.
  std::size_t old_elem = old_mesh.elements.size();
  std::size_t new_elem = new_mesh.elements.size();

  // Number of quadrature nodes.
  std::vector<std::size_t> nqn_new(new_elem);
  nqn_new[0] = 2 * new_mesh.elements[0].degree + 1;

  // Starting indices for old mesh.
  std::vector<std::size_t> old_starts(old_elem);
  old_starts[0] = 0;

  // Starting indices for new mesh.
  std::vector<std::size_t> new_starts(new_elem);
  new_starts[0] = 0;

  for (std::size_t j = 1; j < old_elem; ++j) {
    old_starts[j] = old_starts[j - 1] + old_mesh.elements[j - 1].dofs();
  }
  
  for (std::size_t j = 1; j < new_elem; ++j) {
    new_starts[j] = new_starts[j - 1] + new_mesh.elements[j - 1].dofs();
    nqn_new[j] = 2 * new_mesh.elements[j].degree + 1;
  }

  // Creating new solution vector.
  Vector<Real> new_ch{new_mesh.dofs()};
  std::size_t new_elem_idx = 0;

  // Copy non-refined element.
  for (std::size_t j = 0; j < old_elem; ++j) {
    if (!mask_h[j]) {
      std::size_t start = old_starts[j];
      std::size_t count = old_mesh.elements[j].dofs();
      std::copy(ch.elements.begin() + start,
                ch.elements.begin() + start + count,
                new_ch.elements.begin() + new_starts[new_elem_idx]);
      new_elem_idx++;
    }
  }

  // Loop over old elements.
  for (std::size_t j = 0; j < old_elem; j++) {

    // Refinement check.
    if (!mask_h[j])
      continue;

    // Global matrix indices.
    std::vector<std::size_t> old_indices(old_mesh.elements[j].dofs());
    for (std::size_t k = 0; k < old_mesh.elements[j].dofs(); ++k)
      old_indices[k] = old_starts[j] + k;

    // New elements.
    std::size_t new_elements = (old_mesh.elements[j].edges.size() <= 4)
                                   ? old_mesh.elements[j].edges.size()
                                   : old_mesh.elements[j].edges.size() + 1;

    // Loop over new elements.
    for (std::size_t i = 0; i < new_elements; ++i, ++new_elem_idx) {

      // 2D quadrature weights and nodes.
      auto [nodes_x_2d, nodes_y_2d, weights_2d] =
          quadrature_2d(nqn_new[new_elem_idx]);

      // Local dofs.
      std::vector<std::size_t> new_indices(
          new_mesh.elements[new_elem_idx].dofs());

      for (std::size_t k = 0; k < new_mesh.elements[new_elem_idx].dofs(); ++k)
        new_indices[k] = new_starts[new_elem_idx] + k;

      // Polygon.
      Polygon polygon = new_mesh.element(new_elem_idx);

      // Element sub-triangulation.
      std::vector<Polygon> triangles = triangulate(polygon);

      // Local mass matrix.
      Matrix<Real> local_mass{new_indices.size(), new_indices.size()};

      // Local coefficients.
      Vector<Real> local_coefficients{new_indices.size()};

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Weights scaling.
        Vector<Real> scaled = jacobian_det * weights_2d;

        // Basis functions.
        Matrix<Real> phi =
            basis_2d(new_mesh, new_elem_idx, {physical_x, physical_y})[0];
        Matrix<Real> phi_old =
            basis_2d(old_mesh, j, {physical_x, physical_y})[0];
        Matrix<Real> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_phi.columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // Local coefficients.
        local_mass += scaled_phi.transpose() * phi;
        local_coefficients +=
            scaled_phi.transpose() * (phi_old * ch(old_indices));
      }
      // Update solution vector.
      local_coefficients = solve(local_mass, local_coefficients);
      new_ch(new_indices, local_coefficients);
    }
  }
  // Return update solution.
  return new_ch;
};

} // namespace pacs
