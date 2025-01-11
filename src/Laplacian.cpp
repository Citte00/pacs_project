/**
 * @file laplacian.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Src file with all methods related to the laplace equation problem.
 * @date 2025-01-10
 *
 * @copyright Copyright (c) 2025
 *
 */
#include <PacsHPDG.hpp>

namespace pacs {

/**
 * @brief System matrices assembly.
 *
 * @param data Laplace equation data struct.
 * @param mesh Mesh struct.
 */
void Laplace::assembly(const DataLap &data, const Mesh &mesh) {
#ifndef NVERBOSE
  std::cout << "Computing the laplacian matrix." << std::endl;
#endif

  // Reshape matrices to the new mesh.
  this->m_mass.reshape(mesh.dofs(), mesh.dofs());
  this->m_stiff.reshape(mesh.dofs(), mesh.dofs());
  this->m_DG_stiff.reshape(mesh.dofs(), mesh.dofs());

  // Number of quadrature nodes.
  std::vector<std::size_t> nqn{mesh.elements.size(), 0};
  std::transform(mesh.elements.begin(), mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Degrees of freedom.
  std::size_t dofs = mesh.dofs();

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

  // Matrices.
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
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y);

      // Basis functions.
      auto [phi, gradx_phi, grady_phi] =
          basis_2d(mesh, j, {physical_x, physical_y});

      // Some products.
      Matrix<Real> scaled_gradx{gradx_phi};
      Matrix<Real> scaled_grady{grady_phi};
      Matrix<Real> scaled_phi{phi};

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
    this->m_mass.insert(indices, indices, local_M);
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

      // Param initialization.
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y);

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
  this->m_DG_stiff = A + SA;
  this->m_stiff = this->m_DG_stiff - IA - IA.transpose();

  // Compression.
  this->m_mass.compress();
  this->m_DG_stiff.compress();
  this->m_stiff.compress();
}

/**
 * @brief Extrapolates blocks (indices) based on mass structure.
 *
 * @param mesh Mesh struct.
 * @return std::vector<std::array<std::vector<std::size_t>, 2>>
 */
std::vector<std::array<std::vector<std::size_t>, 2>>
Laplace::block_mass(const Mesh &mesh) const {

#ifndef NVERBOSE
  std::cout << "Evaluating the mass blocks." << std::endl;
#endif

  // Elements.
  std::size_t elements = mesh.elements.size();

  // Blocks.
  std::vector<std::array<std::vector<std::size_t>, 2>> blocks;
  blocks.resize(elements);
  std::size_t start = 0;

  // Precomputing dofs.
  std::vector<std::size_t> dofs;
  dofs.resize(elements);

  for (std::size_t j = 0; j < elements; ++j)
    dofs[j] = mesh.elements[j].dofs();

  // Evaluating blocks.
  for (std::size_t j = 0; j < elements; ++j) {
    std::vector<std::size_t> indices;
    indices.resize(dofs[j]);

    for (std::size_t k = 0; k < dofs[j]; ++k)
      indices[k] = start + k;

    blocks[j] = std::array<std::vector<std::size_t>, 2>{indices, indices};
    start += dofs[j];
  }

  return blocks;
}

/**
 * @brief Assemby the forcing term.
 *
 * @param data Laplace equation data struct.
 * @param mesh Mesh struct.
 */
void Laplace::assemblyforce(const DataLap &data, const Mesh &mesh) {
#ifndef NVERBOSE
  std::cout << "Computing the forcing term." << std::endl;
#endif

  // Resize forcing vector to new mesh.
  this->m_forcing.resize(mesh.dofs());

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
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y);
      Vector<Real> local_source = data.source_f(physical_x, physical_y, D_ext);

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

      // Local matrix assembly.
      Matrix<Real> scaled_gradx{gradx_phi};
      Matrix<Real> scaled_grady{grady_phi};
      Matrix<Real> scaled_phi{phi};

      // Boundary conditions.
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y);
      Vector<Real> boundary = data.DirBC(physical_x, physical_y);

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
    this->m_forcing(indices, local_f);
  }
}

/**
 * @brief Laplace equation solver.
 *
 * @param mesh Mesh Struct.
 * @param TOL Tolerance.
 * @return Vector<Real>
 */
void Laplace::solver(const Mesh &mesh, const Real &TOL) {
  // Mass blocks.
  auto blocks = block_mass(mesh);

  // Solves using BICGSTAB and DBI preconditioner.
  this->m_ch = solve(this->m_stiff, this->m_forcing, blocks, GMRES, DBI, TOL);
}

/**
 * @brief Retrieve coefficients of a generic function.
 *
 * @param mesh Mesh Struct.
 * @param function Generic function.
 * @return Vector<Real>
 */
Vector<Real> Laplace::evaluateCoeff(const Mesh &mesh,
                                    const BiFunctor &function) const {

  // Number of quadrature nodes.
  std::vector<std::size_t> nqn{mesh.elements.size(), 0};
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
      Vector<Real> local_function = function(physical_x, physical_y);

      // Local coefficients.
      local_coefficients += scaled_phi.transpose() * local_function;
    }

    // Update.
    coefficients(indices, local_coefficients);
  }

  return coefficients;
}

/**
 * @brief Retrieve coefficients of the source term.
 *
 * @param data Laplace equation data struct.
 * @param mesh Mesh Struct.
 * @return Vector<Real>
 */
Vector<Real> Laplace::evaluateSource(const DataLap &data,
                                     const Mesh &mesh) const {
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
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y);
      Vector<Real> local_function =
          data.source_f(physical_x, physical_y, D_ext);

      // Local coefficients.
      local_coefficients += scaled_phi.transpose() * local_function;
    }

    // Update.
    coefficients(indices, local_coefficients);
  }

  return coefficients;
}

/**
 * @brief Compute DG and L2 scalar errors.
 *
 * @tparam Functor
 * @param mesh Mesh struct.
 * @param laplace Equation object.
 * @param exact Exact solution.
 */
template <typename Functor>
void LaplaceError::computeError(const Mesh &mesh, const Laplace &laplace,
                                const Functor &exact) {
#ifndef NVERBOSE
  std::cout << "Evaluating errors." << std::endl;
#endif

  // Matrices.
  Sparse<Real> mass = laplace.matrices()[0];
  Sparse<Real> dg_stiff = laplace.matrices()[2];

  // Mass blocks.
  auto blocks = laplace.block_mass(mesh);

  // Error vector.
  Vector<Real> u_modals = laplace.evaluateCoeff(mesh, exact);
  Vector<Real> u_coeff = solve(mass, u_modals, blocks, DB);

  Vector<Real> error = u_coeff - laplace.ch();

  // DG Error.
  this->m_DG_error = std::sqrt(dot(error, dg_stiff * error));

  // L2 Error.
  this->m_L2_error = std::sqrt(dot(error, mass * error));
}

/**
 * @brief Compute laplace equation errors.
 *
 * @param data Laplace equation data struct.
 * @param mesh Mesh struct.
 * @param laplace Laplace equation object.
 */
void LaplaceError::computeErrors(const DataLap &data, const Mesh &mesh,
                                 const Laplace &laplace) {

  // Compute L2 and DG error.
  this->computeError(mesh, laplace, data.c_ex);

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
      Vector<Real> u = data.c_ex(physical_x, physical_y);
      Vector<Real> uh = phi * laplace.ch()(indices);

      Vector<Real> grad_u = data.dc_dx_ex(physical_x, physical_y) +
                            data.dc_dy_ex(physical_x, physical_y);
      Vector<Real> grad_uh = (gradx_phi + grady_phi) * laplace.ch()(indices);

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

} // namespace pacs