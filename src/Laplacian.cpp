/**
 * @file Laplacian.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Contains the implementation of all Laplace objet methods.
 * @date 2025-01-24
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <PacsHPDG.hpp>

namespace pacs {
/**
 * @brief Returns the matrix for the Laplacian operator.
 *
 * @param mesh Mesh.
 * @param penalty_coefficient Penalty coefficient.
 * @return std::array<Sparse<Real>, 3>
 */
void Laplace::assembly(const DataLaplace &data, const Mesh &mesh) {
#ifndef NVERBOSE
  std::cout << "Computing the laplacian matrix." << std::endl;
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
  std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

  // Volume integrals.

  // Loop over the elements.
#pragma omp parallel for schedule(dynamic)
  for (std::size_t j = 0; j < num_elem; ++j) {

    // 2D quadrature nodes and weights.
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
    Matrix<Real> local_A{indices.size(), indices.size()};

    // Loop over the sub-triangulation.
    for (const auto &triangle : triangles) {

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

      // Weights scaling.
      Vector<Real> scaled = jacobian_det * weights_2d;

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      // Some products.
      Matrix<Real> scaled_gradx{gradx_phi};
      Matrix<Real> scaled_grady{grady_phi};
      Matrix<Real> scaled_phi{phi};

      for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
        scaled_gradx.column(l, scaled_gradx.column(l) * scaled);
        scaled_grady.column(l, scaled_grady.column(l) * scaled);
        scaled_phi.column(l, scaled_phi.column(l) * scaled);
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
    std::vector<Matrix<Real>> local_SAN;

    // Penalties.
    Vector<Real> penalties = penalty(mesh, j, data.penalty_coeff);

    // Edges.
    std::vector<Segment> edges{polygon.edges()};

    // Loop over faces.
    for (std::size_t k = 0; k < element_neighbours.size(); ++k) {

      // Neighbour information.
      auto [edge, neighbour, n_edge] = element_neighbours[k];

      // 1D quadrature nodes and weights.
      auto [nodes_1d, weights_1d] =
          (neighbour > 0) ? quadrature_1d(std::max(nqn[j], nqn[neighbour]))
                          : quadrature_1d(nqn[j]);

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

      for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
        scaled_gradx.column(l, scaled_gradx.column(l) * scaled);
        scaled_grady.column(l, scaled_grady.column(l) * scaled);
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
  this->m_mass = std::move(M);
  this->m_dg_stiff = std::move(A) + SA;
  this->m_stiff = std::move(this->m_dg_stiff) - IA - IA.transpose();

  // Compression.
  this->m_mass.compress();
  this->m_dg_stiff.compress();
  this->m_stiff.compress();
};

/**
 * @brief Extrapolates blocks (indices) based on mass structure.
 *
 * @param mesh Mesh.
 * @return std::vector<std::array<std::vector<std::size_t>, 2>>
 */
std::vector<std::array<std::vector<std::size_t>, 2>>
Laplace::block_mass(const Mesh &mesh) const {
#ifndef NVERBOSE
  std::cout << "Evaluating the mass blocks." << std::endl;
#endif

  // Blocks.
  std::vector<std::array<std::vector<std::size_t>, 2>> blocks(
      mesh.elements.size());
  std::size_t start = 0;

  // Precomputing dofs.
  std::vector<std::size_t> dofs(mesh.elements.size());

  for (std::size_t j = 0; j < mesh.elements.size(); ++j)
    dofs[j] = mesh.elements[j].dofs();

  // Evaluating blocks.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {
    
    std::vector<std::size_t> indices(dofs[j]);
    for (std::size_t k = 0; k < dofs[j]; ++k)
      indices[k] = start + k;

    blocks[j] = std::array<std::vector<std::size_t>, 2>{indices, indices};
    start += dofs[j];
  }

  return blocks;
};

/**
 * @brief Assemblies the RHS.
 *
 * @param mesh Mesh.
 * @param penalty_coefficient Penalty coefficient.
 * @return Vector<Real>
 */
Vector<Real> Laplace::assembly_force(const DataLaplace &data,
                                     const Mesh &mesh) const {
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

  // Degrees of freedom.
  std::size_t dofs = mesh.dofs();

  // Forcing term.
  Vector<Real> forcing{dofs};

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

  // Volume integrals.

  // Loop over the elements.
#pragma omp parallel for schedule(dynamic)
  for (std::size_t j = 0; j < num_elem; ++j) {

    // 2D quadrature nodes and weights.
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
      Vector<Real> local_source = data.source_f(physical_x, physical_y);

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

      // 1D quadrature nodes and weights.
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

      for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
        scaled_gradx.column(l, scaled_gradx.column(l) * scaled);
        scaled_grady.column(l, scaled_grady.column(l) * scaled);
        scaled_phi.column(l, scaled_phi.column(l) * scaled);
      }

      Matrix<Real> scaled_grad =
          normal_vector[0] * scaled_gradx + normal_vector[1] * scaled_grady;

      // Boundary conditions.
      Vector<Real> boundary = data.DirBC(physical_x, physical_y);

      // Local forcing term.
      local_f -= scaled_grad.transpose() * boundary;
      local_f += penalties[k] * scaled_phi.transpose() * boundary;
    }

// Global forcing term.
#pragma omp critical
    forcing(indices, local_f);
  }
  return forcing;
};

/**
 * @brief Custom Laplacian solver.
 *
 * @param mesh Mesh.
 * @param A Matrix.
 * @param b Vector.
 * @param TOL Tolerance.
 * @return Vector<Real>
 */
Vector<Real> Laplace::solver(const Mesh &mesh, const Vector<Real> &b,
                             const Real &TOL) const {
#ifndef NVERBOSE
  std::cout << "Solving the algebraic system." << std::endl;
#endif

  // Mass blocks.
  auto blocks = block_mass(mesh);

  // Solves using BICGSTAB and DBI preconditioner.
  return solve(this->m_stiff, b, blocks, GMRES, DBI, TOL);
};

/**
 * @brief Returns the modal coefficients of a function.
 *
 * @param mesh Mesh.
 * @param function Function.
 * @return Vector<Real>
 */
Vector<Real> Laplace::modal(const Mesh &mesh, const BiFunctor &function) const {
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
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // 2D quadrature nodes and weights.
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
      Vector<Real> local_function = function(physical_x, physical_y);

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