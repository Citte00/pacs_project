/**
 * @file Laplacian.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Laplace equation class.
 * @date 2025-01-14
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_SOLVERS_LAPLACIAN_HPP
#define INCLUDE_PACSHPDG_SOLVERS_LAPLACIAN_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"

namespace pacs {

/**
 * @brief Laplace equation class
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class Laplace {
protected:
  Sparse<T> m_mass, m_stiff, m_dg_stiff;

public:
  // CONSTRUCTOR.
  Laplace(const Mesh &mesh_)
      : m_mass{mesh_.dofs(), mesh_.dofs()}, m_stiff{mesh_.dofs(), mesh_.dofs()},
        m_dg_stiff{mesh_.dofs(), mesh_.dofs()} {

    this->m_mass.inner.reserve(DOFS_MAX);
    this->m_mass.outer.reserve(DOFS_MAX);
    this->m_mass.values.reserve(DOFS_MAX);

    this->m_stiff.inner.reserve(DOFS_MAX);
    this->m_stiff.outer.reserve(DOFS_MAX);
    this->m_stiff.values.reserve(DOFS_MAX);

    this->m_dg_stiff.inner.reserve(DOFS_MAX);
    this->m_dg_stiff.outer.reserve(DOFS_MAX);
    this->m_dg_stiff.values.reserve(DOFS_MAX);
  };

  // GETTERS.
  Sparse<T> M() const { return this->m_mass; };
  Sparse<T> A() const { return this->m_stiff; };
  Sparse<T> DG() const { return this->m_dg_stiff; };

  // METHODS.

  /**
   * @brief Extrapolates blocks (indices) based on mass structure.
   *
   * @param mesh_ Mesh structure.
   * @return std::vector<std::array<std::vector<std::size_t>, 2>>
   */
  std::vector<std::array<std::vector<std::size_t>, 2>>
  block_mass(const Mesh &mesh_) const {
#ifndef NVERBOSE
    std::cout << "Evaluating the mass blocks." << std::endl;
#endif

    // Blocks.
    std::vector<std::array<std::vector<std::size_t>, 2>> blocks(
        mesh_.elements.size());
    std::size_t start = 0;

    // Precomputing dofs.
    std::vector<std::size_t> dofs(mesh_.elements.size());

    for (std::size_t j = 0; j < mesh_.elements.size(); ++j)
      dofs[j] = mesh_.elements[j].dofs();

    // Evaluating blocks.
    for (std::size_t j = 0; j < mesh_.elements.size(); ++j) {

      std::vector<std::size_t> indices(dofs[j]);
      for (std::size_t k = 0; k < dofs[j]; ++k)
        indices[k] = start + k;

      blocks[j] = std::array<std::vector<std::size_t>, 2>{indices, indices};
      start += dofs[j];
    }

    return blocks;
  };

  /**
   * @brief Assembly the laplace equation system matrices.
   *
   * @param data_ Laplace equation data structure.
   * @param mesh_ Mesh structure.
   */
  void assembly(const DataLaplace &data_, const Mesh &mesh_) {
#ifndef NVERBOSE
    std::cout << "Computing the laplacian matrix." << std::endl;
#endif

    // Number of elements.
    std::size_t num_elem = mesh_.elements.size();

    // Starting indices.
    std::vector<std::size_t> starts(num_elem);
    starts[0] = 0;

    // Quadrature nodes.
    std::vector<std::size_t> nqn(num_elem);
    nqn[0] = 2 * mesh_.elements[0].degree + 1;

    for (std::size_t j = 1; j < num_elem; ++j) {
      starts[j] = starts[j - 1] + mesh_.elements[j - 1].dofs();
      nqn[j] = 2 * mesh_.elements[j].degree + 1;
    }

    // Degrees of freedom.
    std::size_t dofs = mesh_.dofs();

    // Matrices.
    Sparse<T> M{dofs, dofs};
    Sparse<T> A{dofs, dofs};
    Sparse<T> IA{dofs, dofs};
    Sparse<T> SA{dofs, dofs};

    // Neighbours.
    std::vector<std::vector<std::array<int, 3>>> neighbours = mesh_.neighbours;

    // Volume integrals.

    // Loop over the elements.
#pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < num_elem; ++j) {

      // 2D quadrature nodes and weights.
      auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

      // Global matrix indices.
      std::vector<std::size_t> indices(mesh_.elements[j].dofs());
      for (std::size_t k = 0; k < mesh_.elements[j].dofs(); ++k)
        indices[k] = starts[j] + k;

      // Polygon.
      Polygon polygon = mesh_.element(j);

      // Element sub-triangulation.
      std::vector<Polygon> triangles = triangulate(polygon);

      // Local matrices.
      Matrix<T> local_M{indices.size(), indices.size()};
      Matrix<T> local_A{indices.size(), indices.size()};

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Some products.
        Matrix<T> scaled_gradx{gradx_phi};
        Matrix<T> scaled_grady{grady_phi};
        Matrix<T> scaled_phi{phi};

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
      Matrix<T> local_IA{indices.size(), indices.size()};
      Matrix<T> local_SA{indices.size(), indices.size()};

      // Element's neighbours.
      std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

      // Local matrices for neighbours.
      std::vector<Matrix<T>> local_IAN;
      std::vector<Matrix<T>> local_SAN;

      // Penalties.
      Vector<Real> penalties = penalty(mesh_, j, data_.penalty_coeff);

      // Edges.
      std::vector<Segment> edges{polygon.edges()};

      // Loop over faces.
      for (std::size_t k = 0; k < element_neighbours.size(); ++k) {

        // Neighbour information.
        auto [edge, neighbour, n_edge] = element_neighbours[k];

        // 1D quadrature nodes and weights.
        auto [nodes_1d, weights_1d] = quadrature_1d(nqn[j]);

        auto [normal_vector, edge_vector, physical_x, physical_y] =
            faces_physical_points(edges[k], nodes_1d);

        // Weights scaling.
        Vector<T> scaled = std::abs(edges[k]) * weights_1d;

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Local matrix assembly.
        Matrix<T> scaled_gradx{gradx_phi};
        Matrix<T> scaled_grady{grady_phi};
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
          scaled_gradx.column(l, scaled_gradx.column(l) * scaled);
          scaled_grady.column(l, scaled_grady.column(l) * scaled);
          scaled_phi.column(l, scaled_phi.column(l) * scaled);
        }

        Matrix<T> scaled_grad =
            normal_vector[0] * scaled_gradx + normal_vector[1] * scaled_grady;

        if (neighbour == -1) { // Boundary edge.

          local_IA += scaled_grad.transpose() * phi;
          local_SA += (penalties[k] * scaled_phi).transpose() * phi;

          // Empty small matrices.
          local_IAN.emplace_back(Matrix<T>{1, 1});
          local_SAN.emplace_back(Matrix<T>{1, 1});

        } else {

          local_IA += 0.5 * scaled_grad.transpose() * phi;
          local_SA += (penalties[k] * scaled_phi).transpose() * phi;

          // Neighbour's basis function.
          Matrix<T> n_phi =
              basis_2d(mesh_, neighbour, {physical_x, physical_y})[0];

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
        std::size_t n_dofs =
            mesh_.elements[n_index].dofs(); // Neighbour's dofs.

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
   * @brief Assembly the forcing term.
   *
   * @param data_ Laplace equation data structure.
   * @param mesh_ Mesh structure.
   * @return Vector<T>
   */
  Vector<T> assembly_force(const DataLaplace &data_, const Mesh &mesh_) const {
#ifndef NVERBOSE
    std::cout << "Computing the forcing term." << std::endl;
#endif

    // Number of elements.
    std::size_t num_elem = mesh_.elements.size();

    // Starting indices.
    std::vector<std::size_t> starts(num_elem);
    starts[0] = 0;

    // Quadrature nodes.
    std::vector<std::size_t> nqn(num_elem);
    nqn[0] = 2 * mesh_.elements[0].degree + 1;

    for (std::size_t j = 1; j < num_elem; ++j) {
      starts[j] = starts[j - 1] + mesh_.elements[j - 1].dofs();
      nqn[j] = 2 * mesh_.elements[j].degree + 1;
    }

    // Degrees of freedom.
    std::size_t dofs = mesh_.dofs();

    // Forcing term.
    Vector<T> forcing{dofs};

    // Neighbours.
    std::vector<std::vector<std::array<int, 3>>> neighbours = mesh_.neighbours;

    // Volume integrals.

    // Loop over the elements.
#pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < num_elem; ++j) {

      // 2D quadrature nodes and weights.
      auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

      // Global matrix indices.
      std::vector<std::size_t> indices(mesh_.elements[j].dofs());
      for (std::size_t k = 0; k < mesh_.elements[j].dofs(); ++k)
        indices[k] = starts[j] + k;

      // Polygon.
      Polygon polygon = mesh_.element(j);

      // Element sub-triangulation.
      std::vector<Polygon> triangles = triangulate(polygon);

      // Local forcing term.
      Vector<T> local_f{indices.size()};

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Local source evaluation.
        Vector<T> local_source = data_.source_f(physical_x, physical_y);

        // Basis functions.
        auto phi = basis_2d(mesh_, j, {physical_x, physical_y})[0];
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_phi.columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // Local forcing term.
        local_f += scaled_phi.transpose() * local_source;
      }

      // Face integrals.

      // Element's neighbours.
      std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

      // Penalties.
      Vector<Real> penalties = penalty(mesh_, j, data_.penalty_coeff);

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
        Vector<T> scaled = std::abs(edges[k]) * weights_1d;

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Local matrix assembly.
        Matrix<T> scaled_gradx{gradx_phi};
        Matrix<T> scaled_grady{grady_phi};
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
          scaled_gradx.column(l, scaled_gradx.column(l) * scaled);
          scaled_grady.column(l, scaled_grady.column(l) * scaled);
          scaled_phi.column(l, scaled_phi.column(l) * scaled);
        }

        Matrix<T> scaled_grad =
            normal_vector[0] * scaled_gradx + normal_vector[1] * scaled_grady;

        // Boundary conditions.
        Vector<T> boundary = data_.DirBC(physical_x, physical_y);

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
   * @brief Laplace equation solver.
   *
   * @param mesh_ Mesh structure.
   * @param RHS_ System rhs term vector.
   * @param tol_ Tolerance.
   * @return Vector<T>
   */
  Vector<T> solver(const Mesh &mesh_, const Vector<T> &RHS_,
                   const T &tol_ = 1E-15) const {
#ifndef NVERBOSE
    std::cout << "Solving the algebraic system." << std::endl;
#endif

    // Mass blocks.
    auto blocks = block_mass(mesh_);

    // Solves using BICGSTAB and DBI preconditioner.
    return solve(this->m_stiff, RHS_, blocks, GMRES, DBI, tol_);
  };

  /**
   * @brief Get functions modal coefficients.
   *
   * @param mesh_ Mesh structure.
   * @param function_ Function to evaluate.
   * @return Vector<T>
   */
  Vector<T>
  modal(const Mesh &mesh_,
        const Functor<Vector<T>, Vector<T>, Vector<T>> &function_) const {
#ifndef NVERBOSE
    std::cout << "Retrieving modal coefficients." << std::endl;
#endif

    // Number of elements.
    std::size_t num_elem = mesh_.elements.size();

    // Starting indices.
    std::vector<std::size_t> starts(num_elem);
    starts[0] = 0;

    // Quadrature nodes.
    std::vector<std::size_t> nqn(num_elem);
    nqn[0] = 2 * mesh_.elements[0].degree + 1;

    for (std::size_t j = 1; j < num_elem; ++j) {
      starts[j] = starts[j - 1] + mesh_.elements[j - 1].dofs();
      nqn[j] = 2 * mesh_.elements[j].degree + 1;
    }

    // Coefficients.
    Vector<T> coefficients{mesh_.dofs()};

// Loop over the elements.
#pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < mesh_.elements.size(); ++j) {

      // 2D quadrature nodes and weights.
      auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

      // Global matrix indices.
      std::vector<std::size_t> indices(mesh_.elements[j].dofs());
      for (std::size_t k = 0; k < mesh_.elements[j].dofs(); ++k)
        indices[k] = starts[j] + k;

      // Polygon.
      Polygon polygon = mesh_.element(j);

      // Element sub-triangulation.
      std::vector<Polygon> triangles = triangulate(polygon);

      // Local coefficients.
      Vector<T> local_coefficients{indices.size()};

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Basis functions.
        Matrix<T> phi = basis_2d(mesh_, j, {physical_x, physical_y})[0];
        Matrix<T> scaled_phi = phi;

        for (std::size_t l = 0; l < scaled_phi.columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // function solution.
        Vector<T> local_function = function_(physical_x, physical_y);

        // Local coefficients.
        local_coefficients += scaled_phi.transpose() * local_function;
      }

// Update.
#pragma omp critical
      coefficients(indices, local_coefficients);
    }

    auto blocks = block_mass(mesh_);

    return (norm(coefficients) > TOLERANCE)
               ? solve(this->m_mass, coefficients, blocks, DB)
               : Vector<T>{mesh_.dofs()};
  };
};
} // namespace pacs

#endif