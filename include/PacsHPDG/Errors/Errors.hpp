/**
 * @file Errors.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Errors computation classes.
 * @date 2024-12-22
 *
 * @copyright Copyright (c) 2024
 *
 */
#ifndef INCLUDE_PACSHPDG_ERRORS_ERRORS_HPP
#define INCLUDE_PACSHPDG_ERRORS_ERRORS_HPP

#include <array>
#include <iostream>

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"
#include "../Solvers.hpp"

namespace pacs {

/**
 * @brief Laplace equation errors computation class.
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class LaplaceError {
protected:
  std::size_t m_degree, m_dofs;
  T m_size;

  T m_dg_error;
  T m_l2_error;

  Vector<T> m_l2_errors;
  Vector<T> m_h1_errors;
  Vector<T> m_dg_errors;

public:
  // CONSTRUCTOR.
  LaplaceError(const Mesh &mesh_)
      : m_l2_errors{mesh_.elements.size()}, m_h1_errors{mesh_.elements.size()},
        m_dg_errors{mesh_.elements.size()} {

    this->m_dofs = mesh_.dofs();

    this->m_degree = 0;
    for (const auto &element : mesh_.elements)
      this->m_degree =
          (element.degree > this->m_degree) ? element.degree : this->m_degree;

    this->m_size = 0.0;
    for (const auto &element : mesh_.elements)
      for (const auto &p : element.element.points)
        for (const auto &q : element.element.points)
          this->m_size =
              (distance(p, q) > this->m_size) ? distance(p, q) : this->m_size;
  };

  // GETTERS.
  std::size_t dofs() const { return this->m_dofs; };
  std::size_t p() const { return this->m_degree; };
  T h() const { return this->m_size; };
  T L2error() const { return this->m_l2_error; };
  T DGerror() const { return this->m_dg_error; };
  Vector<T> L2errors() const { return this->m_l2_errors; };
  Vector<T> H1errors() const { return this->m_h1_errors; };
  Vector<T> DGerrors() const { return this->m_dg_errors; };

  // METHODS.
  /**
   * @brief Compute Laplace equation L2 and DG errors.
   *
   * @param data_ Laplace equation data structure.
   * @param mesh_ Mesh structure.
   * @param laplacian_ Laplace equation object.
   * @param ch_ Numerical solution.
   */
  void error(const DataLaplace &data_, const Mesh &mesh_,
             const Laplace<T> &laplacian_, const Vector<T> &ch_) {
#ifndef NVERBOSE
    std::cout << "Evaluating L2 and DG errors." << std::endl;
#endif

    // Matrices.
    Sparse<T> mass = laplacian_.M();
    Sparse<T> dg_stiff = laplacian_.DG();

    // Error vector.
    Vector<T> u_modals = laplacian_.modal(mesh_, data_.c_ex);
    Vector<T> error = u_modals - ch_;

    // DG Error.
    this->m_dg_error = std::sqrt(dot(error, dg_stiff * error));

    // L2 Error.
    this->m_l2_error = std::sqrt(dot(error, mass * error));
  };

  /**
   * @brief Compute Laplace equation L2 and H1 errors.
   *
   * @param data_ Laplace equation data structure.
   * @param mesh_ Mesh structure.
   * @param laplacian_ Laplace equation object.
   * @param ch_ Numerical solution.
   */
  void errors(const DataLaplace &data_, const Mesh &mesh_,
              const Laplace<T> &laplacian_, const Vector<T> &ch_) {
#ifndef NVERBOSE
    std::cout << "Evaluating L2 and H1 errors." << std::endl;
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

      // 2D quadrature nodes and weights.
      auto [nodes_1d, weights_1d, nodes_2d, weights_2d] = Quadrature(nqn[j]);

      // Global matrix indices.
      std::vector<std::size_t> indices(mesh_.elements[j].dofs());

      for (std::size_t k = 0; k < mesh_.elements[j].dofs(); ++k)
        indices[k] = starts[j] + k;

      // Polygon.
      Polygon polygon = mesh_.element(j);

      // Element sub-triangulation.
      std::vector<Polygon> triangles = triangulate(polygon);

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, nodes_2d);

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Solutions.
        Vector<T> u = data_.c_ex(physical_x, physical_y);
        Vector<T> uh = phi * ch_(indices);

        Vector<T> grad_u = data_.dc_dx_ex(physical_x, physical_y) +
                           data_.dc_dy_ex(physical_x, physical_y);
        Vector<T> grad_uh = (gradx_phi + grady_phi) * ch_(indices);

        // Local L2 error.
        this->m_l2_errors[j] += dot(scaled, (u - uh) * (u - uh));

        // Local H1 error.
        this->m_h1_errors[j] +=
            dot(scaled, (grad_u - grad_uh) * (grad_u - grad_uh));
      }

      // Element's neighbours.
      std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

      // Penalties.
      Vector<Real> penalties = penalty(mesh_, j, data_.penalty_coeff);

      // Edges.
      std::vector<Segment> edges{polygon.edges()};

      // Jump error.
      T jump = 0.0;

      // Loop over faces.
      for (std::size_t k = 0; k < element_neighbours.size(); ++k) {

        // Physical points.
        auto [normal_vector, edge_vector, physical_x, physical_y] =
            faces_physical_points(edges[k], nodes_1d);

        // Weights scaling.
        Vector<T> scaled = std::abs(edges[k]) * weights_1d;

        // Basis functions.
        auto phi = basis_2d(mesh_, j, {physical_x, physical_y})[0];
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_phi.m_columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // Local solutions.
        Vector<T> u = data_.c_ex(physical_x, physical_y);
        Vector<T> uh = phi * ch_(indices);

        // Jump error.
        jump += penalties[k] * dot(scaled, (uh - u) * (uh - u));
      }

      this->m_l2_errors[j] = std::sqrt(this->m_l2_errors[j]);
      this->m_h1_errors[j] = std::sqrt(this->m_h1_errors[j]);
      this->m_dg_errors[j] = this->m_h1_errors[j] + std::sqrt(jump);
    }
  };

  /**
   * @brief Friend operator<< for polymorphic printing.
   *
   * @param ost File on which to print the error.
   * @param error Error object.
   * @return std::ostream&
   */
  friend std::ostream &operator<<(std::ostream &ost,
                                  const LaplaceError<T> &error) {
    ost << "Elements: " << error.L2errors().size() << "\n";
    ost << "Dofs: " << error.dofs() << "\n";
    ost << "Degree (p): " << error.p() << "\n";
    ost << "Size (h): " << error.h() << "\n";
    ost << "L2 Error: " << error.L2error() << "\n";
    return ost << "DG Error: " << error.DGerror() << std::endl;
  };
};

/**
 * @brief Heat equation errors computation class.
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class HeatError : public LaplaceError<T> {
public:
  // CONSTRUCTOR.
  HeatError(const Mesh &mesh_) : LaplaceError<T>(mesh_) {};

  // METHODS.
  /**
   * @brief Compute heat equation L2 and DG errors.
   *
   * @param data_ Heat equation data struct.
   * @param mesh_ Mesh struct.
   * @param heat_ Heat equation object.
   * @param ch_ Numerical solution.
   */
  void error(const DataHeat &data_, const Mesh &mesh_, const Heat<T> &heat_,
             const Vector<T> &ch_) {
#ifndef NVERBOSE
    std::cout << "Evaluating L2 and DG errors." << std::endl;
#endif

    // Matrices.
    Sparse<T> mass = heat_.M();
    Sparse<T> dg_stiff = heat_.DG();

    // Error vector.
    Vector<T> u_modals = heat_.modal(mesh_, data_.c_ex);
    Vector<T> error = u_modals - ch_;

    // DG Error.
    this->m_dg_error = std::sqrt(dot(error, dg_stiff * error));

    // L2 Error.
    this->m_l2_error = std::sqrt(dot(error, mass * error));

    // Number of elements.
    std::size_t num_elem = mesh_.elements.size();

    // Starting indices.
    std::vector<std::size_t> starts(num_elem);
    starts[0] = 0;

    for (std::size_t j = 1; j < num_elem; ++j)
      starts[j] = starts[j - 1] + mesh_.elements[j - 1].dofs();

      // Loop over the elements.
#pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < num_elem; ++j) {

      // Global matrix indices.
      std::vector<std::size_t> indices(mesh_.elements[j].dofs());

      for (std::size_t k = 0; k < mesh_.elements[j].dofs(); ++k)
        indices[k] = starts[j] + k;

      Matrix<T> stiff = dg_stiff(indices, indices);
      Vector<T> e = error(indices);
      this->m_dg_errors[j] = std::sqrt(dot(e, stiff * e));
    }
  };

  /**
   * @brief Compute heat equation L2 and H1 errors.
   *
   * @param data_ Heat equation data struct.
   * @param mesh_ Mesh struct.
   * @param heat_ Heat equation object.
   * @param ch_ Numerical solution.
   */
  void errors(const DataHeat &data_, const Mesh &mesh_, const Heat<T> &heat_,
              const Vector<T> &ch_) {
#ifndef NVERBOSE
    std::cout << "Evaluating L2 and H1 errors." << std::endl;
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

      // 2D quadrature nodes and weights.
      auto [nodes_1d, weights_1d, nodes_2d, weights_2d] = Quadrature(nqn[j]);

      // Global matrix indices.
      std::vector<std::size_t> indices(mesh_.elements[j].dofs());

      for (std::size_t k = 0; k < mesh_.elements[j].dofs(); ++k)
        indices[k] = starts[j] + k;

      // Polygon.
      Polygon polygon = mesh_.element(j);

      // Element sub-triangulation.
      std::vector<Polygon> triangles = triangulate(polygon);

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, nodes_2d);

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Solutions.
        Vector<T> u = data_.c_ex(physical_x, physical_y, heat_.t());
        Vector<T> uh = phi * ch_(indices);

        Vector<T> grad_u = data_.dc_dx_ex(physical_x, physical_y, heat_.t()) +
                           data_.dc_dy_ex(physical_x, physical_y, heat_.t());
        Vector<T> grad_uh = (gradx_phi + grady_phi) * ch_(indices);

        // Local L2 error.
        this->m_l2_errors[j] += dot(scaled, (u - uh) * (u - uh));

        // Local H1 error.
        this->m_h1_errors[j] +=
            dot(scaled, (grad_u - grad_uh) * (grad_u - grad_uh));
      }

      // Element's neighbours.
      std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

      // Penalties.
      Vector<Real> penalties = penalty(mesh_, j, data_.penalty_coeff);

      // Edges.
      std::vector<Segment> edges{polygon.edges()};

      // Jump error.
      T jump = 0.0;

      // Loop over faces.
      for (std::size_t k = 0; k < element_neighbours.size(); ++k) {

        // Physical points.
        auto [normal_vector, edge_vector, physical_x, physical_y] =
            faces_physical_points(edges[k], nodes_1d);

        // Weights scaling.
        Vector<T> scaled = std::abs(edges[k]) * weights_1d;

        // Basis functions.
        auto phi = basis_2d(mesh_, j, {physical_x, physical_y})[0];

        Vector<T> uh = phi * ch_(indices);

        if( element_neighbours[k][1] == -1){

          // Boundary term.
          Vector<T> g = data_.DirBC(physical_x, physical_y, heat_.t());

          jump += penalties[k] * dot(scaled, (uh - g) * (uh - g));
        
        } else {

          // Neighbour's basis function.
          auto [n_phi, n_gradx_phi, n_grady_phi] =
              basis_2d(mesh, neighbour, {physical_x, physical_y});

          std::size_t n_index = element_neighbours[k][1];
          std::size_t n_dofs =
              mesh.elements[n_index].dofs(); // Neighbour's dofs.

          std::vector<std::size_t> n_indices(n_dofs);
          for (std::size_t h = 0; h < n_dofs; ++h)
            n_indices[h] = starts[n_index] + h;

          // Neighbour's numerical solution and gradients.
          Vector<T> n_uh = n_phi * ch_(n_indices);

          // Jump error.
          jump += penalties[k] * dot(scaled, (uh - n_uh) * (uh - n_uh));
        }
      }

      this->m_l2_errors[j] = std::sqrt(this->m_l2_errors[j]);
      this->m_h1_errors[j] = std::sqrt(this->m_h1_errors[j]);
      this->m_dg_errors[j] = this->m_h1_errors[j] + std::sqrt(jump);
    }
  };
};

/**
 * @brief Fisher equation errors computation class.
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class FisherError : public HeatError<T> {
protected:
  // Energy error.
  T m_energy;

public:
  // CONSTRUCTOR.
  FisherError(const Mesh &mesh_) : HeatError<T>(mesh_), m_energy{0.0} {};

  // GETTER.
  T energy() const { return this->m_energy; };
  T &energy() { return this->m_energy; };

  // METHODS.

  /**
   * @brief Compute Fisher-KPP equation L2, DG and energy errors.
   *
   * @param data_ Fisher equation data struct.
   * @param mesh_ Mesh struct.
   * @param fisher_ Fisher equation object.
   * @param ch_ Numerical Solution.
   */
  void error(const DataFKPP &data_, const Mesh &mesh_, const Fisher<T> &fisher_,
             const Vector<T> &ch_) {
#ifndef NVERBOSE
    std::cout << "Evaluating L2, DG and energy errors." << std::endl;
#endif

    // Matrices.
    Sparse<T> mass = fisher_.M();
    Sparse<T> dg_stiff = fisher_.DG();

    // Error vector.
    Vector<T> u_modals = fisher_.modal(mesh_, data_.c_ex);
    Vector<T> error = u_modals - ch_;

    // DG Error.
    this->m_dg_error = std::sqrt(dot(error, dg_stiff * error));

    // L2 Error.
    this->m_l2_error = std::sqrt(dot(error, mass * error));

    // Energy error.
    this->m_energy += data_.dt * this->m_dg_error * this->m_dg_error;
  };

  /**
   * @brief Compute Fisher-KPP equation L2 and H1 errors.
   *
   * @param data_ Fisher equation data struct.
   * @param mesh_ Mesh struct.
   * @param fisher_ Fisher equation object.
   * @param ch_ Numerical Solution.
   */
  void errors(const DataFKPP &data_, const Mesh &mesh_,
              const Fisher<T> &fisher_, const Vector<T> &ch_) {
#ifndef NVERBOSE
    std::cout << "Evaluating L2 and H1 errors." << std::endl;
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

      // 2D quadrature nodes and weights.
      auto [nodes_1d, weights_1d, nodes_2d, weights_2d] = Quadrature(nqn[j]);

      // Global matrix indices.
      std::vector<std::size_t> indices(mesh_.elements[j].dofs());

      for (std::size_t k = 0; k < mesh_.elements[j].dofs(); ++k)
        indices[k] = starts[j] + k;

      // Polygon.
      Polygon polygon = mesh_.element(j);

      // Element sub-triangulation.
      std::vector<Polygon> triangles = triangulate(polygon);

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, nodes_2d);

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh_, j, {physical_x, physical_y});

        // Solutions.
        Vector<T> u = data_.c_ex(physical_x, physical_y, fisher_.t());
        Vector<T> uh = phi * ch_(indices);

        Vector<T> grad_u = data_.dc_dx_ex(physical_x, physical_y, fisher_.t()) +
                           data_.dc_dy_ex(physical_x, physical_y, fisher_.t());
        Vector<T> grad_uh = (gradx_phi + grady_phi) * ch_(indices);

        // Local L2 error.
        this->m_l2_errors[j] += dot(scaled, (u - uh) * (u - uh));

        // Local H1 error.
        this->m_h1_errors[j] +=
            dot(scaled, (grad_u - grad_uh) * (grad_u - grad_uh));
      }

      // Element's neighbours.
      std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

      // Penalties.
      Vector<Real> penalties = penalty(mesh_, j, data_.penalty_coeff);

      // Edges.
      std::vector<Segment> edges{polygon.edges()};

      // Jump error.
      T jump = 0.0;

      // Loop over faces.
      for (std::size_t k = 0; k < element_neighbours.size(); ++k) {

        // Physical points.
        auto [normal_vector, edge_vector, physical_x, physical_y] =
            faces_physical_points(edges[k], nodes_1d);

        // Weights scaling.
        Vector<T> scaled = std::abs(edges[k]) * weights_1d;

        // Basis functions.
        auto phi = basis_2d(mesh_, j, {physical_x, physical_y})[0];
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_phi.m_columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // Local solutions.
        Vector<T> u = data_.c_ex(physical_x, physical_y, fisher_.t());
        Vector<T> uh = phi * ch_(indices);

        // Jump error.
        jump += penalties[k] * dot(scaled, (uh - u) * (uh - u));
      }

      this->m_l2_errors[j] = std::sqrt(this->m_l2_errors[j]);
      this->m_h1_errors[j] = std::sqrt(this->m_h1_errors[j]);
      this->m_dg_errors[j] = this->m_h1_errors[j] + std::sqrt(jump);
    }
  };

  /**
   * @brief Friend operator<< for polymorphic printing.
   *
   * @param ost File on which to print the error.
   * @param error Error object.
   * @return std::ostream&
   */
  friend std::ostream &operator<<(std::ostream &ost,
                                  const FisherError<T> &error) {
    ost << "Elements: " << error.L2errors().size() << std::endl;
    ost << "Dofs: " << error.dofs() << std::endl;
    ost << "Degree (p): " << error.p() << std::endl;
    ost << "Size (h): " << error.h() << std::endl;
    ost << "L2 Error: " << error.L2error() << std::endl;
    ost << "DG Error: " << error.DGerror() << std::endl;
    return ost << "Energy Error: "
               << std::sqrt(std::pow(error.L2error(), 2) + error.energy())
               << std::endl;
  };
};

} // namespace pacs

#endif