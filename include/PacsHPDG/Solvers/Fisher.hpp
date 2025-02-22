/**
 * @file Fisher.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief
 * @date 2024-11-28
 *
 * @copyright Copyright (c) 2024
 *
 */
#ifndef INCLUDE_PACSHPDG_SOLVERS_FISHER_HPP
#define INCLUDE_PACSHPDG_SOLVERS_FISHER_HPP

#include "./Heat.hpp"

namespace pacs {

/**
 * @brief Fisher-KPP equation class
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class Fisher : public Heat<T> {
protected:
  // Non-linear matrix.
  Sparse<T> m_nl_mass;
  Vector<T> m_ch_old;

public:
  // CONSTRUCTOR.
  Fisher(const DataFKPP &data_, const Mesh &mesh_)
      : Heat<T>(mesh_), m_nl_mass{mesh_.dofs(), mesh_.dofs()},
        m_ch_old{mesh_.dofs()} {
    this->m_t = -data_.dt;

    // Reserve space for nl_mass matrix.
    this->m_nl_mass.inner.reserve(DOFS_MAX);
    this->m_nl_mass.outer.reserve(DOFS_MAX);
    this->m_nl_mass.values.reserve(DOFS_MAX);

    // Reserve space for ch_old vector.
    this->m_ch_old.elements.reserve(DOFS_MAX);
  };

  // GETTER.
  Sparse<T> M_alpha() const { return this->m_nl_mass; };
  Vector<T> ch_old() const { return this->m_ch_old; };
  Vector<T> &ch_old() { return this->m_ch_old; };

  // METHODS.
  /**
   * @brief Update matrices and forcing term in adaptive framework.
   *
   * @param data_ Fisher-KPP equation data structure.
   * @param mesh_ Mesh structure.
   */
  void update(const DataFKPP &data_, const Mesh &mesh_) {
#ifndef NVERBOSE
    std::cout << "Updating the Fisher-KPP equation matrices." << std::endl;
#endif

    // Dofs.
    std::size_t dofs = mesh_.dofs();

    // Update matrices.
    this->m_mass.reshape(dofs, dofs);
    this->m_stiff.reshape(dofs, dofs);
    this->m_dg_stiff.reshape(dofs, dofs);
    this->m_nl_mass.reshape(dofs, dofs);

    assembly(data_, mesh_);

    // Update forcing term.
    this->m_forcing.resize(dofs);
    assembly_force(data_, mesh_);
  };

  /**
   * @brief Assembly the fisher equation system matrices.
   *
   * @param data_ Fisher-KPP equation data structure.
   * @param mesh_ Mesh structure.
   */
  void assembly(const DataFKPP &data_, const Mesh &mesh_) {
#ifndef NVERBOSE
    std::cout << "Computing the Fisher-KPP matrices." << std::endl;
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
    Sparse<T> M_prj{dofs, dofs};
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

      // 2D Quadrature nodes and weights.
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
      Matrix<T> local_M_prj{indices.size(), indices.size()};
      Matrix<T> local_M{indices.size(), indices.size()};
      Matrix<T> local_A{indices.size(), indices.size()};

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Param initialization.
        Vector<T> Dext = data_.D_ext(physical_x, physical_y, 0.0);
        Vector<T> alpha = data_.alpha(physical_x, physical_y, 0.0);

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Some products.
        Matrix<T> scaled_phi{phi};
        Matrix<T> scaled_phi_alpha{phi};
        Matrix<T> scaled_gradx{gradx_phi};
        Matrix<T> scaled_grady{grady_phi};

        for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
          scaled_phi.column(l, scaled_phi.column(l) * scaled);
          scaled_phi_alpha.column(l, (alpha * scaled_phi_alpha.column(l)) *
                                         scaled);
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

        // 1D Quadrature nodes and weights.
        auto [nodes_1d, weights_1d] =
            (neighbour > 0) ? quadrature_1d(std::max(nqn[j], nqn[neighbour]))
                            : quadrature_1d(nqn[j]);

        auto [normal_vector, edge_vector, physical_x, physical_y] =
            faces_physical_points(edges[k], nodes_1d);

        // Weights scaling.
        Vector<T> scaled = std::abs(edges[k]) * weights_1d;

        // Param initialization.
        Vector<T> Dext = data_.D_ext(physical_x, physical_y, 0.0);

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Local matrix assembly.
        Matrix<T> scaled_gradx{gradx_phi};
        Matrix<T> scaled_grady{grady_phi};
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
          scaled_gradx.column(l, (Dext * scaled_gradx.column(l)) * scaled);
          scaled_grady.column(l, (Dext * scaled_grady.column(l)) * scaled);
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
   * @brief Assembly the non-linear term.
   *
   * @param data_ Fisher-KPP equation data structure.
   * @param mesh_ Mesh structure.
   * @param ch_ Numerical solution.
   * @return Sparse<T>
   */
  Sparse<T> assembly_nl(const DataFKPP &data_, const Mesh &mesh_,
                        const Vector<T> &ch_) const {
#ifndef NVERBOSE
    std::cout << "Computing the Fisher-KPP matrices." << std::endl;
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
    Sparse<T> M_star{dofs, dofs};

    // Neighbours.
    std::vector<std::vector<std::array<int, 3>>> neighbours = mesh_.neighbours;

// Loop over the elements.
#pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < num_elem; ++j) {

      // 2D Quadrature nodes and weights.
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

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Param initialization.
        Vector<T> alpha = data_.alpha(physical_x, physical_y, 0.0);
        Vector<T> c_star = ch_(indices);

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Some products.
        Matrix<T> scaled_phi{phi};
        Vector<T> c_star_loc = phi * c_star;

        for (std::size_t l = 0; l < scaled_phi.columns; ++l) {
          scaled_phi.column(l, (alpha * scaled_phi.column(l)) * c_star_loc *
                                   scaled);
        }

        // Local matrix assembly.
        local_M += scaled_phi.transpose() * phi;
      }

// Global matrix assembly.
#pragma omp critical
      M_star.insert(indices, indices, local_M);
    }

    // Matrices.
    Sparse<T> non_lin = std::move(M_star);

    // Compression.
    non_lin.compress();

    return non_lin;
  };

  /**
   * @brief Assembly the forcing term.
   *
   * @param data_ Fisher-KPP equation data structure.
   * @param mesh_ Mesh structure.
   */
  void assembly_force(const DataFKPP &data_, const Mesh &mesh_) {
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

    // Loop over the elements.
#pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < num_elem; ++j) {

      // 2D Local quadrature nodes and weights.
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
        Vector<T> Dext = data_.D_ext(physical_x, physical_y, this->m_t);
        Vector<T> alpha = data_.alpha(physical_x, physical_y, this->m_t);
        Vector<T> local_source =
            data_.source_f(physical_x, physical_y, this->m_t, Dext, alpha);

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
        Vector<T> scaled = std::abs(edges[k]) * weights_1d;

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Param initialisation.
        Vector<T> Dext = data_.D_ext(physical_x, physical_y, this->m_t);

        // Local matrix assembly.
        Matrix<T> scaled_gradx{gradx_phi};
        Matrix<T> scaled_grady{grady_phi};
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_gradx.columns; ++l) {
          scaled_gradx.column(l, (Dext * scaled_gradx.column(l)) * scaled);
          scaled_grady.column(l, (Dext * scaled_grady.column(l)) * scaled);
          scaled_phi.column(l, scaled_phi.column(l) * scaled);
        }

        Matrix<T> scaled_grad =
            normal_vector[0] * scaled_gradx + normal_vector[1] * scaled_grady;

        // Boundary conditions.
        Vector<T> boundary = data_.DirBC(physical_x, physical_y, this->m_t);

        // Local forcing term.
        local_f -= scaled_grad.transpose() * boundary;
        local_f += penalties[k] * scaled_phi.transpose() * boundary;
      }

// Global forcing term.
#pragma omp critical
      this->m_forcing(indices, local_f);
    }
  };

  //
  /**
   * @brief Fisher-KPP equation solver.
   *
   * @param data_ Fisher-KPP equation data structure.
   * @param mesh_ Mesh structure.
   * @param ch_oold_ Numerical solution at previous time step.
   * @param forcing_old_ Forcing term at previous time step.
   * @param tol_ Tolerance.
   * @return Vector<T>
   */
  Vector<T> solver(const DataFKPP &data_, const Mesh &mesh_,
                   const Vector<T> &ch_oold_, const Vector<T> &forcing_old_,
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

    // Assembling the dynamic componenet of the marices.
    if (data_.theta == 0) {

      Sparse<T> MN = assembly_nl(data_, mesh_, this->m_ch_old);
      RHS += data_.dt * (this->m_mass - MN);

    } else if (data_.theta == 1) {

      Sparse<T> MN = assembly_nl(data_, mesh_, this->m_ch_old);
      LHS -= data_.dt * (this->m_mass - MN);

    } else if (data_.theta == 0.5) {

      Vector<T> c_star = 1.5 * this->m_ch_old - 0.5 * ch_oold_;
      Sparse<T> MN = assembly_nl(data_, mesh_, c_star);
      LHS -= 0.5 * data_.dt * (this->m_mass - MN);
      RHS += 0.5 * data_.dt * (this->m_mass - MN);
    }

    LHS.compress();
    RHS.compress();

    // Construction of the complete RHS for the theta method.
    Vector<T> F = RHS * this->m_ch_old +
                  data_.dt * data_.theta * this->m_forcing +
                  data_.dt * (1.0 - data_.theta) * forcing_old_;

    // Solves using BICGSTAB.
    return solve(LHS, F, blocks, BICGSTAB, DBI, tol_);
  };

  /**
   * @brief Get source function modal coefficient.
   *
   * @param data_ Fisher-KPP equation data structure.
   * @param mesh_ Mesh structure.
   * @return Vector<T>
   */
  Vector<T> modal_source(const DataFKPP &data_, const Mesh &mesh_) const {
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

        for (std::size_t l = 0; l < scaled_phi.columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // function solution.
        Vector<T> D_ext = data_.D_ext(physical_x, physical_y, this->m_t);
        Vector<T> alpha = data_.alpha(physical_x, physical_y, this->m_t);
        Vector<T> local_function =
            data_.source_f(physical_x, physical_y, this->m_t, D_ext, alpha);

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
};
} // namespace pacs

#endif