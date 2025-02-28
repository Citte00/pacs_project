/**
 * @file Heat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief
 * @date 2024-12-16
 *
 * @copyright Copyright (c) 2024
 *
 */
#ifndef INCLUDE_PACSHPDG_SOLVERS_HEAT_HPP
#define INCLUDE_PACSHPDG_SOLVERS_HEAT_HPP

#include "./Laplacian.hpp"

namespace pacs {

/**
 * @brief Heat equation class
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class Heat : public Laplace<T> {
protected:
  // Forcing term and numerical solution.
  Vector<T> m_forcing;

  // Time step.
  T m_t;

public:
  // CONSTRUCTOR.
  Heat(const Mesh &mesh_)
      : Laplace<T>(mesh_), m_forcing{mesh_.dofs()}, m_t{0.0} {
    this->m_forcing.elements.reserve(DOFS_MAX);
  };

  // GETTERS.
  Vector<T> forcing() const { return this->m_forcing; };
  Vector<T> &forcing() { return this->m_forcing; };
  T t() const { return this->m_t; };
  T &t() { return this->m_t; };

  // METHODS.

  /**
   * @brief Update matrices and forcing term in adaptive framework.
   *
   * @param data_ Heat equation data structure.
   * @param mesh_ Mesh structure.
   */
  void update(const DataHeat &data_, const Mesh &mesh_) {
#ifndef NVERBOSE
    std::cout << "Updating the Heat equation matrices." << std::endl;
#endif

    // Dofs.
    std::size_t dofs = mesh_.dofs();

    // Update matrices.
    this->m_mass.reshape(dofs, dofs);
    this->m_stiff.reshape(dofs, dofs);
    this->m_dg_stiff.reshape(dofs, dofs);

    assembly(data_, mesh_);

    // Update forcing term.
    this->m_forcing.resize(dofs);

    assembly_force(data_, mesh_);
  };

  /**
   * @brief Assembly the heat equation system matrices.
   *
   * @param data_ Heat equation data structure.
   * @param mesh_ Mesh structure.
   */
  void assembly(const DataHeat &data_, const Mesh &mesh_) {
#ifndef NVERBOSE
    std::cout << "Computing the Heat equation matrices." << std::endl;
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
    const std::vector<std::vector<std::array<int, 3>>> &neighbours =
        mesh_.neighbours;

    // Volume integrals.

    // Loop over the elements.
#pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < num_elem; ++j) {

      // Quadrature nodes and weights.
      auto [nodes_1d, weights_1d, nodes_2d, weights_2d] = Quadrature(nqn[j]);

      // Global matrix indices.
      std::vector<std::size_t> indices(mesh_.elements[j].dofs());
      for (std::size_t k = 0; k < mesh_.elements[j].dofs(); ++k)
        indices[k] = starts[j] + k;

      // Polygon.
      const Polygon &polygon = mesh_.element(j);

      // Element sub-triangulation.
      std::vector<Polygon> triangles = triangulate(polygon);

      // Local matrices.
      Matrix<T> local_M{indices.size(), indices.size()};
      Matrix<T> local_A{indices.size(), indices.size()};

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, nodes_2d);

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Param initialization.
        Vector<T> D_ext = data_.D_ext(physical_x, physical_y, 0.0);

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Some products.
        Matrix<T> scaled_phi{phi};
        Matrix<T> scaled_gradx{gradx_phi};
        Matrix<T> scaled_grady{grady_phi};

        for (std::size_t l = 0; l < scaled_gradx.m_columns; ++l) {
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
      Matrix<T> local_IA{indices.size(), indices.size()};
      Matrix<T> local_SA{indices.size(), indices.size()};

      // Element's neighbours.
      std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

      // Local matrices for neighbours.
      std::vector<Matrix<T>> local_IAN;
      local_IAN.reserve(element_neighbours.size());
      std::vector<Matrix<T>> local_SAN;
      local_SAN.reserve(element_neighbours.size());

      // Penalties.
      Vector<Real> penalties = penalty(mesh_, j, data_.penalty_coeff);

      // Edges.
      std::vector<Segment> edges{polygon.edges()};

      // Loop over faces.
      for (std::size_t k = 0; k < element_neighbours.size(); ++k) {
        // Neighbour information.
        auto [edge, neighbour, n_edge] = element_neighbours[k];

        // Physical points.
        auto [normal_vector, edge_vector, physical_x, physical_y] =
            faces_physical_points(edges[k], nodes_1d);

        // Weights scaling.
        Vector<T> scaled = std::abs(edges[k]) * weights_1d;

        // Param initialization.
        Vector<T> D_ext = data_.D_ext(physical_x, physical_y, 0.0);

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Local matrix assembly.
        Matrix<T> scaled_gradx{gradx_phi};
        Matrix<T> scaled_grady{grady_phi};
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_gradx.m_columns; ++l) {
          scaled_gradx.column(l, (D_ext * scaled_gradx.column(l)) * scaled);
          scaled_grady.column(l, (D_ext * scaled_grady.column(l)) * scaled);
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

        std::size_t n_index = element_neighbours[k][1];
        std::size_t n_dofs =
            mesh_.elements[n_index].dofs(); // Neighbour's dofs.
        std::vector<std::size_t> n_indices(n_dofs);

        for (std::size_t h = 0; h < n_dofs; ++h)
          n_indices[h] = starts[n_index] + h;

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
   * @brief Assembly the forcing term.
   *
   * @param data_ Heat equation data structure.
   * @param mesh_ Mesh structure.
   */
  void assembly_force(const DataHeat &data_, const Mesh &mesh_) {
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

    // Neighbours.
    std::vector<std::vector<std::array<int, 3>>> neighbours = mesh_.neighbours;

    // Volume integrals.
    Vector<T> forcing{mesh_.dofs()};

// Loop over the elements.
#pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < num_elem; ++j) {

      // Local quadrature nodes and weights.
      auto [nodes_1d, weights_1d, nodes_2d, weights_2d] = Quadrature(nqn[j]);

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
            get_Jacobian_physical_points(triangle, nodes_2d);

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Local source evaluation.
        Vector<T> D_ext = data_.D_ext(physical_x, physical_y, m_t);
        Vector<T> local_source =
            data_.source_f(physical_x, physical_y, m_t, D_ext);

        // Basis functions.
        auto phi = basis_2d(mesh_, j, {physical_x, physical_y})[0];
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_phi.m_columns; ++l)
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

        // Boundary conditions.
        Vector<T> D_ext = data_.D_ext(physical_x, physical_y, m_t);
        Vector<T> boundary = data_.DirBC(physical_x, physical_y, m_t);

        for (std::size_t l = 0; l < scaled_gradx.m_columns; ++l) {
          scaled_gradx.column(l, (D_ext * scaled_gradx.column(l)) * scaled);
          scaled_grady.column(l, (D_ext * scaled_grady.column(l)) * scaled);
          scaled_phi.column(l, scaled_phi.column(l) * scaled);
        }

        Matrix<T> scaled_grad =
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
   * @brief Heat equation solver.
   *
   * @param data_ Heat equation data structure.
   * @param mesh_ Mesh structure.
   * @param ch_old_ Numerical solution at previuos time step.
   * @param forcing_old_ Forcing term at previous time step.
   * @param tol_ Tolerance.
   * @return Vector<T>
   */
  Vector<T> solver(const DataHeat &data_, const Mesh &mesh_,
                   const Vector<T> &ch_old_, const Vector<T> &forcing_old_,
                   const T &tol_ = 1E-15) const {
#ifndef NVERBOSE
    std::cout << "Solving the algebraic system." << std::endl;
#endif

    // Mass blocks.
    auto blocks = this->block_mass(mesh_);

    // Assembling the constant component of the matrices.
    Sparse<T> LHS = this->m_mass + data_.dt * data_.theta * this->m_stiff;
    Sparse<T> RHS =
        this->m_mass - data_.dt * (1.0 - data_.theta) * this->m_stiff;

    LHS.compress();
    RHS.compress();

    // Construction of the complete RHS for the theta method.
    Vector<T> F = RHS * ch_old_ + data_.dt * data_.theta * this->m_forcing +
                  data_.dt * (1 - data_.theta) * forcing_old_;

    // Solves using BICGSTAB.
    return solve(LHS, F, blocks, GMRES, DBI, tol_);
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
        const Functor<Vector<T>, Vector<T>, Vector<T>, T> &function_) const {
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
    for (std::size_t j = 0; j < num_elem; ++j) {

      // Quadrature nodes.
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

        for (std::size_t l = 0; l < scaled_phi.m_columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // function solution.
        Vector<T> local_function = function_(physical_x, physical_y, m_t);

        // Local coefficients.
        local_coefficients += scaled_phi.transpose() * local_function;
      }

      // Update.
#pragma omp critical
      coefficients(indices, local_coefficients);
    }

    auto blocks = this->block_mass(mesh_);

    return (norm(coefficients) > TOLERANCE)
               ? solve(this->m_mass, coefficients, blocks, DB)
               : Vector<T>{mesh_.dofs()};
  };

  /**
   * @brief Get source function modal coefficient.
   *
   * @param data_ Heat equation data structure.
   * @param mesh_ Mesh structure.
   * @return Vector<T>
   */
  Vector<T> modal_source(const DataHeat &data_, const Mesh &mesh_) const {
#ifndef NVERBOSE
    std::cout << "Retrieving modal coeffficient of the source term."
              << std::endl;
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
    for (std::size_t j = 0; j < num_elem; ++j) {
      // Quadrature nodes.
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

        for (std::size_t l = 0; l < scaled_phi.m_columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // function solution.
        Vector<T> D_ext = data_.D_ext(physical_x, physical_y, m_t);
        Vector<T> local_function =
            data_.source_f(physical_x, physical_y, m_t, D_ext);

        // Local coefficients.
        local_coefficients += scaled_phi.transpose() * local_function;
      }

      // Update.
#pragma omp critical
      coefficients(indices, local_coefficients);
    }

    auto blocks = this->block_mass(mesh_);

    return (norm(coefficients) > TOLERANCE)
               ? solve(this->m_mass, coefficients, blocks, DB)
               : Vector<T>{mesh_.dofs()};
  };

  // hp-adaptive methods.
  /**
   * @brief Prolong solution to p-refined mesh.
   *
   * @param new_mesh_ Refined mesh structure.
   * @param old_mesh_ Unrefined mesh structure.
   * @param ch_ Numerical solution to prolong.
   * @param p_mask_ Vector identifying element to refine through degree.
   * @return Vector<T>
   */
  Vector<T> prolong_solution_p(const Mesh &new_mesh_, const Mesh &old_mesh_,
                               const Vector<T> &ch_,
                               const Mask &p_mask_) const {
#ifndef NVERBOSE
    std::cout << "Prolong solution to new dofs." << std::endl;
#endif

    // Number of elements.
    std::size_t old_elem = old_mesh_.elements.size();
    std::size_t new_elem = new_mesh_.elements.size();

    // Number of quadrature nodes.
    std::vector<std::size_t> nqn(new_elem);
    nqn[0] = 2 * new_mesh_.elements[0].degree + 1;

    // Resize new solution vector.
    Vector<T> new_ch(new_mesh_.dofs());

    // Starting indices.
    std::vector<std::size_t> old_starts(old_elem);
    old_starts[0] = 0;

    for (std::size_t j = 1; j < old_elem; ++j)
      old_starts[j] = old_starts[j - 1] + old_mesh_.elements[j - 1].dofs();

    std::vector<std::size_t> new_starts(new_elem);
    new_starts[0] = 0;

    for (std::size_t j = 1; j < new_elem; ++j) {
      new_starts[j] = new_starts[j - 1] + new_mesh_.elements[j - 1].dofs();
      nqn[j] = 2 * new_mesh_.elements[j].degree + 1;
    }

// Loop over elements.
#pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < new_elem; j++) {

      // 2D quadrature weights and nodes.
      auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

      // Local dofs.
      std::size_t old_element_dofs = old_mesh_.elements[j].dofs();
      std::size_t new_element_dofs = new_mesh_.elements[j].dofs();

      // Global matrix indices.
      std::vector<std::size_t> old_indices(old_element_dofs);
      for (std::size_t k = 0; k < old_element_dofs; ++k)
        old_indices[k] = old_starts[j] + k;

      std::vector<std::size_t> new_indices(new_element_dofs);
      for (std::size_t k = 0; k < new_element_dofs; ++k)
        new_indices[k] = new_starts[j] + k;

      // Check to refine.
      if (!p_mask_[j]) {
        new_ch(new_indices, ch_(old_indices));
        continue;
      }

      // Polygon.
      Polygon polygon = new_mesh_.element(j);

      // Element sub-triangulation.
      std::vector<Polygon> triangles = triangulate(polygon);

      // Local mass matrix.
      Matrix<T> local_mass{new_indices.size(), new_indices.size()};

      // Local coefficients.
      Vector<T> local_coefficients{new_indices.size()};

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Basis functions.
        Matrix<T> phi = basis_2d(new_mesh_, j, {physical_x, physical_y})[0];
        Matrix<T> phi_old = basis_2d(old_mesh_, j, {physical_x, physical_y})[0];
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_phi.m_columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // Local coefficients.
        local_mass += scaled_phi.transpose() * phi;
        local_coefficients +=
            scaled_phi.transpose() * (phi_old * ch_(old_indices));
      }

      // Update the solution.
      local_coefficients = solve(local_mass, local_coefficients);
#pragma omp critical
      new_ch(new_indices, local_coefficients);
    }
    return new_ch;
  };

  /**
   * @brief Prolong solution to h-refined mesh.
   *
   * @param new_mesh_ Refined mesh structure.
   * @param old_mesh_ Unrefined mesh structure.
   * @param ch_ Numerical solution to prolong.
   * @param h_mask_ Vector identifying element to refine through size.
   * @return Vector<T>
   */
  Vector<T> prolong_solution_h(const Mesh &new_mesh_, const Mesh &old_mesh_,
                               const Vector<T> &ch_,
                               const Mask &h_mask_) const {
#ifndef NVERBOSE
    std::cout << "Prolong solution to new elements." << std::endl;
#endif

    // Number of elements.
    std::size_t old_elems = old_mesh_.elements.size();
    std::size_t new_elems = new_mesh_.elements.size();

    // Number of quadrature nodes.
    std::vector<std::size_t> nqn_new(new_elems);
    nqn_new[0] = 2 * new_mesh_.elements[0].degree + 1;

    // Starting indices for old mesh.
    std::vector<std::size_t> old_starts(old_elems);
    old_starts[0] = 0;

    // Starting indices for new mesh.
    std::vector<std::size_t> new_starts(new_elems);
    new_starts[0] = 0;

    for (std::size_t j = 1; j < old_elems; ++j) {
      old_starts[j] = old_starts[j - 1] + old_mesh_.elements[j - 1].dofs();
    }

    for (std::size_t j = 1; j < new_elems; ++j) {
      new_starts[j] = new_starts[j - 1] + new_mesh_.elements[j - 1].dofs();
      nqn_new[j] = 2 * new_mesh_.elements[j].degree + 1;
    }

    // Creating new solution vector.
    Vector<T> new_ch{new_mesh_.dofs()};
    std::size_t new_elem_idx = 0;

    // Copy non-refined element.
    for (std::size_t j = 0; j < old_elems; ++j) {
      if (!h_mask_[j]) {
        std::size_t start = old_starts[j];
        std::size_t count = old_mesh_.elements[j].dofs();
        std::copy(ch_.elements.begin() + start,
                  ch_.elements.begin() + start + count,
                  new_ch.elements.begin() + new_starts[new_elem_idx]);
        new_elem_idx++;
      }
    }

    // Store where each old element's new elements begin andhow many of them are
    // created.
    std::vector<std::size_t> new_elem_offset(old_elems, 0);
    std::vector<std::size_t> new_elements(old_elems, 0);

    for (std::size_t j = 0; j < old_elems; ++j) {
      if (h_mask_[j]) {
        std::size_t elements = (old_mesh_.elements[j].edges.size() <= 4)
                                   ? old_mesh_.elements[j].edges.size()
                                   : old_mesh_.elements[j].edges.size() + 1;
        new_elements[j] = elements;
        new_elem_offset[j] = new_elem_idx;
        new_elem_idx += elements;
      } else {
        new_elem_offset[j] = j; // Preserved elements keep their index.
      }
    }

// Loop over old elements.
#pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < old_elems; j++) {

      // Refinement check.
      if (!h_mask_[j])
        continue;

      // Global matrix indices.
      std::vector<std::size_t> old_indices(old_mesh_.elements[j].dofs());
      for (std::size_t k = 0; k < old_mesh_.elements[j].dofs(); ++k)
        old_indices[k] = old_starts[j] + k;

      // Loop over new elements.
      for (std::size_t i = 0; i < new_elements[j]; ++i) {

        // Indices.
        std::size_t index = new_elem_offset[j] + i;

        // 2D quadrature weights and nodes.
        auto [nodes_x_2d, nodes_y_2d, weights_2d] =
            quadrature_2d(nqn_new[index]);

        // Local dofs.
        std::vector<std::size_t> new_indices(new_mesh_.elements[index].dofs());

        for (std::size_t k = 0; k < new_mesh_.elements[index].dofs(); ++k)
          new_indices[k] = new_starts[index] + k;

        // Polygon.
        Polygon polygon = new_mesh_.element(index);

        // Element sub-triangulation.
        std::vector<Polygon> triangles = triangulate(polygon);

        // Local mass matrix.
        Matrix<T> local_mass{new_indices.size(), new_indices.size()};

        // Local coefficients.
        Vector<T> local_coefficients{new_indices.size()};

        // Loop over the sub-triangulation.
        for (const auto &triangle : triangles) {

          // Jacobian's determinant and physical nodes.
          auto [jacobian_det, physical_x, physical_y] =
              get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

          // Weights scaling.
          Vector<T> scaled = jacobian_det * weights_2d;

          // Basis functions.
          Matrix<T> phi =
              basis_2d(new_mesh_, index, {physical_x, physical_y})[0];
          Matrix<T> phi_old =
              basis_2d(old_mesh_, j, {physical_x, physical_y})[0];
          Matrix<T> scaled_phi{phi};

          for (std::size_t l = 0; l < scaled_phi.m_columns; ++l)
            scaled_phi.column(l, scaled_phi.column(l) * scaled);

          // Local coefficients.
          local_mass += scaled_phi.transpose() * phi;
          local_coefficients +=
              scaled_phi.transpose() * (phi_old * ch_(old_indices));
        }
        // Update solution vector.
        local_coefficients = solve(local_mass, local_coefficients);
#pragma omp critical
        new_ch(new_indices, local_coefficients);
      }
    }
    // Return update solution.
    return new_ch;
  };
};
} // namespace pacs

#endif