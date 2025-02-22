/**
 * @file Estimators.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Error estimates computation classes.
 * @date 2025-01-20
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_ERRORS_ESTIMATORS_HPP
#define INCLUDE_PACSHPDG_ERRORS_ESTIMATORS_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"
#include "../Solvers.hpp"

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace pacs {

/**
 * @brief Laplace equation error estimator class.
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class LaplaceEstimator {
protected:
  // DOFs.
  std::size_t m_dofs;
  // Estimates.
  T m_estimate;
  Vector<T> m_estimates, m_fits;
  // Refined mesh.
  Mesh m_mesh;

public:
  // CONSTRUCTOR.
  LaplaceEstimator(const Mesh &mesh_)
      : m_estimates{mesh_.elements.size()}, m_fits{mesh_.elements.size()},
        m_mesh{mesh_} {
    this->m_dofs = mesh_.dofs();
    this->m_estimate = 0.0;
  };

  // GETTERS.
  std::size_t dofs() const { return this->m_dofs; };
  T estimate() const { return this->m_estimate; };
  Vector<T> estimates() const { return this->m_estimates; };
  Vector<T> fits() const { return this->m_fits; };
  Mesh mesh() const { return this->m_mesh; };

  // METHODS.

  /**
   * @brief Polynomial fit.
   *
   * @param x_ X.
   * @param y_ Y.
   * @param p_ Algorithm order.
   * @return Vector<T>
   */
  Vector<T> polyfit(const Vector<T> &x_, const Vector<T> &y_,
                    const std::size_t &p_) const {
#ifndef NDEBUG // Integrity check.
    assert(p_ > 0);
    assert(x_.length == y_.length);
#endif

    // X.
    Matrix<T> X{x_.length, p_ + 1};

    // Building X.
    for (std::size_t j = 0; j < p_ + 1; ++j) {
      X.column(j, Vector<T>{x_.length, 1.0});

      for (std::size_t k = 0; k < j; ++k) {
        X.column(j, X.column(j) * x_);
      }
    }

    // Solution.
    return solve(X.transpose() * X, X.transpose() * y_, QRD);
  };

  /**
   * @brief Compute Laplace equation error estimates.
   *
   * @param data_ Laplace equation data structure.
   * @param laplacian Laplace equation object.
   * @param numerical Numerical solution.
   */
  void computeEstimates(const DataLaplace &data_, const Laplace<T> &laplacian_,
                        const Vector<T> &ch_) {
#ifndef NVERBOSE
    std::cout << "Evaluating estimates." << std::endl;
#endif
    // Number of elements.
    Mesh mesh = this->m_mesh;
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

    // Sizes.
    Vector<T> sizes{num_elem};

    for (std::size_t j = 0; j < sizes.length; ++j) {
      Element element{mesh.elements[j]};

      for (const auto &p : element.element.points)
        for (const auto &q : element.element.points)
          sizes[j] = (distance(p, q) > sizes[j]) ? distance(p, q) : sizes[j];
    }

    // Mass blocks.
    auto blocks = laplacian_.block_mass(mesh);

    // Coefficients.
    Vector<T> f_modals = laplacian_.modal(mesh, data_.source_f);
    Vector<T> g_modals = laplacian_.modal(mesh, data_.DirBC);

    // Estimate.
    T estimate = 0.0;

    // Loop over the elements.
#pragma omp parallel for reduction(+ : estimate) schedule(dynamic)
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

      // Thread-local variable to accumulate estimates.
      T local_estimate = 0.0;

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Basis functions.
        Matrix<T> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];
        Matrix<T> lap_phi = lap_basis_2d(mesh, j, {physical_x, physical_y});

        // Local numerical laplacian.
        Vector<T> lap_uh = lap_phi * ch_(indices);

        // Local exact source.
        Vector<T> f = data_.source_f(physical_x, physical_y);

        // Local source approximation.
        Vector<T> f_bar = phi * f_modals(indices);

        // Local estimator, R_{K, E}^2.
        local_estimate += sizes[j] * sizes[j] *
                          dot(scaled, (f_bar + lap_uh) * (f_bar + lap_uh));

        // Local data oscillation, O_{K, E}^2.
        local_estimate +=
            sizes[j] * sizes[j] * dot(scaled, (f - f_bar) * (f - f_bar));
      }

      // Element's neighbours.
      std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

      // Penalties.
      Vector<Real> penalties = penalty(mesh, j, data_.penalty_coeff);

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
        Vector<T> scaled = std::abs(edges[k]) * weights_1d;

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh, j, {physical_x, physical_y});
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_phi.columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // Local numerical solution and gradients.
        Vector<T> uh = phi * ch_(indices);

        Matrix<T> grad =
            normal_vector[0] * gradx_phi + normal_vector[1] * grady_phi;
        Vector<T> grad_uh = grad * ch_(indices);

        Matrix<T> grad_t =
            edge_vector[0] * gradx_phi + edge_vector[1] * grady_phi;
        Vector<T> grad_uh_t = grad_t * ch_(indices);

        if (neighbour == -1) { // Boundary edge.

          // Local exact Dirichlet and gradient.
          Vector<T> g = data_.DirBC(physical_x, physical_y);
          Vector<T> grad_g_t =
              edge_vector[0] * data_.dc_dx_ex(physical_x, physical_y) +
              edge_vector[1] * data_.dc_dy_ex(physical_x, physical_y);

          // Approximate Dirichlet and gradient.
          Vector<T> g_bar = phi * g_modals(indices);
          Vector<T> grad_g_t_bar = grad_t * g_modals(indices);

          // Local estimator, R_{K, J}^2.
          local_estimate +=
              penalties[k] * dot(scaled, (uh - g_bar) * (uh - g_bar));

          // // Local estimator, R_{K, N}^2.
          // this->m_estimates[j] += sizes[j] * dot(scaled, grad_uh * grad_uh);

          // Local estimator, R_{K, T}^2.
          local_estimate +=
              sizes[j] * dot(scaled, (grad_uh_t - grad_g_t_bar) *
                                         (grad_uh_t - grad_g_t_bar));

          // Local data oscillation, O_{K, J}^2.
          local_estimate +=
              penalties[k] * dot(scaled, (g - g_bar) * (g - g_bar));

          // Local data oscillation, O_{K, T}^2.
          local_estimate +=
              sizes[j] * dot(scaled, (grad_g_t - grad_g_t_bar) *
                                         (grad_g_t - grad_g_t_bar));
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

          Matrix<T> n_grad =
              normal_vector[0] * n_gradx_phi + normal_vector[1] * n_grady_phi;
          Vector<T> n_grad_uh = n_grad * ch_(n_indices);

          Matrix<T> n_grad_t =
              edge_vector[0] * n_gradx_phi + edge_vector[1] * n_grady_phi;
          Vector<T> n_grad_uh_t = n_grad_t * ch_(n_indices);

          // Local estimator, R_{K, J}^2.
          local_estimate +=
              penalties[k] * dot(scaled, (uh - n_uh) * (uh - n_uh));

          // Local estimator, R_{K, N}^2.
          local_estimate += sizes[j] * dot(scaled, (grad_uh - n_grad_uh) *
                                                       (grad_uh - n_grad_uh));

          // Local estimator, R_{K, T}^2.
          local_estimate +=
              sizes[j] * dot(scaled, (grad_uh_t - n_grad_uh_t) *
                                         (grad_uh_t - n_grad_uh_t));
        }
      }

      this->m_estimates[j] = local_estimate;
      estimate += local_estimate;

      // Degrees.
      Vector<T> degrees{indices.size()};
      std::size_t counter = 0;

      for (std::size_t i = 0; i < mesh.elements[j].degree + 1; ++i)
        for (std::size_t k = 0; k < mesh.elements[j].degree + 1 - i; ++k) {
          degrees[counter] = static_cast<T>(i + k);
          ++counter;
        }

      // Coefficients.
      Vector<T> coefficients = ch_(indices);

      for (auto &coefficient : coefficients.elements)
        coefficient = std::log(std::abs(coefficient));

      // Fit.
      Vector<T> fit = this->polyfit(degrees, coefficients, 1);
      this->m_fits[j] = -fit[1];
    }

    this->m_estimate = std::sqrt(estimate);
  };

  // Refinement.
  /**
   * @brief Refines specified elements' size for a mesh.
   *
   * @param mask_ Vector to identify which element to refine.
   */
  void mesh_refine_size(const Mask &mask_) {
#ifndef NDEBUG // Integrity check.
    assert(mask_.size() == this->m_mesh.elements.size());
#endif

#ifndef NVERBOSE
    std::cout << "Refining mesh size: " << std::flush;
#endif

    // Degrees.
    std::vector<std::size_t> degrees;

    for (std::size_t j = 0; j < mask_.size(); ++j)
      if (!(mask_[j]))
        degrees.emplace_back(this->m_mesh.elements[j].degree);

    // Polygons.
    std::vector<Polygon> diagram;
    std::vector<Polygon> refine;
    std::vector<Polygon> refined;

    for (std::size_t j = 0; j < mask_.size(); ++j)
      if (mask_[j]) {
        refine.emplace_back(this->m_mesh.element(j));

        for (std::size_t k = 0;
             k < this->m_mesh.elements[j].edges.size() +
                     static_cast<std::size_t>(
                         this->m_mesh.elements[j].edges.size() > 4);
             ++k)
          degrees.emplace_back(this->m_mesh.elements[j].degree);

      } else
        diagram.emplace_back(this->m_mesh.element(j));

    // Refine.
    for (const auto &polygon : refine) {
      Point centroid = polygon.centroid();
      std::vector<Point> points;

      for (const auto &edge : polygon.edges()) {

        // New point.
        Point point = (edge[0] + edge[1]) * 0.5;
        points.emplace_back(point);

        for (auto &element : diagram) {
          std::vector<Segment> edges = element.edges();

          for (std::size_t k = 0; k < edges.size(); ++k)
            if (edges[k] == edge) {

              // Diagram editing.
              element.points.insert(element.points.begin() + k + 1, point);
              break;
            }
        }
      }

      // New polygons.
      std::vector<Segment> edges = polygon.edges();
      std::vector<Point> centrals;

      for (std::size_t j = 0; j < points.size(); ++j) {
        Point first = points[j];
        Point second = (j < points.size() - 1) ? points[j + 1] : points[0];

        std::vector<Point> vertices{first, edges[j][1], second};

        if (edges.size() <= 4)
          vertices.emplace_back(centroid);
        else {
          vertices.emplace_back((second + centroid) * 0.5);
          vertices.emplace_back((first + centroid) * 0.5);
          centrals.emplace_back((second + centroid) * 0.5);
        }

        refined.emplace_back(Polygon{vertices});
      }

      // Central Polygon.
      if (edges.size() > 4)
        refined.emplace_back(Polygon{centrals});
    }

    // Update.
    for (const auto &polygon : refined)
      diagram.emplace_back(polygon);

#ifndef NVERBOSE
    std::cout << refine.size() << " elements." << std::endl;
#endif

    // Refinement.
    const Mesh new_mesh{this->m_mesh.domain, diagram, degrees};
    this->m_mesh = std::move(new_mesh);
  };

  /**
   * @brief Refines specified elements' degree for a mesh.
   *
   * @param mask_ Vector to identify which element to refine.
   */
  void mesh_refine_degree(const Mask &mask_) {
#ifndef NDEBUG // Integrity check.
    assert(mask_.size() == this->m_mesh.elements.size());
#endif

#ifndef NVERBOSE
    std::cout << "Refining mesh degree: " << std::flush;
    std::size_t counter = 0;
#endif

    for (std::size_t j = 0; j < mask_.size(); ++j)
      if (mask_[j]) {
        ++(this->m_mesh).elements[j].degree;

#ifndef NVERBOSE
        ++counter;
#endif
      }

#ifndef NVERBOSE
    std::cout << counter << " elements." << std::endl;
#endif
  };

  /**
   * @brief Choose elements to refine.
   *
   * @param refine_ Refinement percentage.
   * @param speed_ Solution's smoothness.
   * @return std::array<Mask, 2>
   */
  std::array<Mask, 2> find_elem_to_refine(const T &refine_ = 0.75,
                                          const T &speed_ = 1.0) const {
#ifndef NDEBUG // Integrity check.
    assert((refine_ > 0.0L) && (refine_ < 1.0L));
    assert(speed_ > 0.0L);
#endif

    // Masks.
    Mask p_mask = this->m_fits > speed_;
    Mask h_mask = this->m_estimates > refine_ * sum(this->m_estimates) /
                                          this->m_mesh.elements.size();

    // Strategy.
    for (std::size_t j = 0; j < this->m_mesh.elements.size(); ++j) {
      if (!h_mask[j] || this->m_mesh.elements[j].degree >=
                            10) // p-Refine only error-marked elements.
        p_mask[j] = false;

      if (p_mask[j] && h_mask[j]) // p > h.
        h_mask[j] = false;
    }

    return {h_mask, p_mask};
  };

  /**
   * @brief Friend operator<< for output printing.
   *
   * @param ost File on which to print the error estimates.
   * @param estimator Estimator object.
   * @return std::ostream&
   */
  friend std::ostream &operator<<(std::ostream &ost,
                                  const LaplaceEstimator<T> &estimator) {
    ost << "Dofs: " << estimator.dofs() << std::endl;
    ost << "Estimate: " << estimator.estimate() << std::endl;
    return ost;
  }

  /**
   * @brief Outputs the estimates distribution to a .poly file.
   *
   * @param filename_ File name.
   */
  void write(const std::string &filename_) {
    // File loading.
    std::ofstream file{filename_};

    file << "@ " << filename_ << "\n";
    file << "@ polyplot.py readable mesh\n";

    // Stats.
    file << "@ Domain: " << this->m_mesh.domain << "\n";
    file << "@ Elements number: " << this->m_mesh.elements.size() << "\n";

    // Polygons.
    file << "@ Elements: \n";

    int count = 0;

    for (const auto &element : this->m_mesh.elements) {
      Polygon polygon = element.element;

      for (const auto &vertex : polygon.points)
        file << std::setprecision(12) << vertex[0] << " "
             << std::setprecision(12) << vertex[1] << " ";

      file << this->m_estimates[count];

      file << "\n";
      count++;
    }
  };
};

/**
 * @brief Heat equation error estimator class.
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class HeatEstimator : public LaplaceEstimator<T> {
public:
  // CONSTRUCTOR.
  HeatEstimator(const Mesh &mesh_) : LaplaceEstimator<T>(mesh_) {};

  // METHODS.

  /**
   * @brief Compute heat equation error estimates.
   *
   * @param data_ Heat equation data structure.
   * @param heat_ Heat equation object.
   * @param ch_ Numerical solution at current time step.
   * @param ch_old_ Numerical solution at previous time step.
   */
  void computeEstimates(const DataHeat &data_, const Heat<T> &heat_,
                        const Vector<T> &ch_, const Vector<T> &ch_old_) {
#ifndef NVERBOSE
    std::cout << "Evaluating estimates." << std::endl;
#endif

    // Number of elements.
    Mesh mesh = this->m_mesh;
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

    // Sizes.
    Vector<T> sizes{num_elem};

    for (std::size_t j = 0; j < sizes.length; ++j) {
      Element element{mesh.elements[j]};

      for (const auto &p : element.element.points)
        for (const auto &q : element.element.points)
          sizes[j] = (distance(p, q) > sizes[j]) ? distance(p, q) : sizes[j];
    }

    // Coefficients.
    Vector<T> f_modals = heat_.modal_source(data_, mesh);
    Vector<T> g_modals = heat_.modal(mesh, data_.DirBC);

    // Estimate.
    T estimate = 0.0;

    // Loop over the elements.
#pragma omp parallel for reduction(+ : estimate) schedule(dynamic)
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

      // Local estimate.
      T local_estimate = 0.0;

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Basis functions.
        Matrix<T> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];
        Matrix<T> lap_phi = lap_basis_2d(mesh, j, {physical_x, physical_y});

        // Local exact source.
        Vector<T> D_ext = data_.D_ext(physical_x, physical_y, heat_.t());
        Vector<T> f = data_.source_f(physical_x, physical_y, heat_.t(), D_ext);

        // Local numerical heat.
        Vector<T> lap_uh = D_ext * (lap_phi * ch_(indices));

        // Local time derivative.
        Vector<T> local_uh = phi * ch_(indices);
        Vector<T> local_uh_old = phi * ch_old_(indices);
        Vector<T> partial_uh_t = (local_uh - local_uh_old) / data_.dt;

        // Local source approximation.
        Vector<T> f_bar = phi * f_modals(indices);

        // Local estimator, R_{K, E}^2.
        local_estimate += sizes[j] * sizes[j] *
                          dot(scaled, (f_bar - partial_uh_t + lap_uh) *
                                          (f_bar - partial_uh_t + lap_uh));

        // Local data oscillation, O_{K, E}^2.
        local_estimate +=
            sizes[j] * sizes[j] * dot(scaled, (f - f_bar) * (f - f_bar));
      }

      // Element's neighbours.
      std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

      // Penalties.
      Vector<Real> penalties = penalty(mesh, j, data_.penalty_coeff);

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
        Vector<T> scaled = std::abs(edges[k]) * weights_1d;

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh, j, {physical_x, physical_y});
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_phi.columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // Local numerical solution and gradients.
        Vector<T> uh = phi * ch_(indices);

        Matrix<T> grad =
            normal_vector[0] * gradx_phi + normal_vector[1] * grady_phi;
        Vector<T> grad_uh = grad * ch_(indices);

        Matrix<T> grad_t =
            edge_vector[0] * gradx_phi + edge_vector[1] * grady_phi;
        Vector<T> grad_uh_t = grad_t * ch_(indices);

        if (neighbour == -1) { // Boundary edge.

          // Local exact Dirichlet and gradient.
          Vector<T> g = data_.DirBC(physical_x, physical_y, heat_.t());
          Vector<T> grad_g_t =
              edge_vector[0] *
                  data_.dc_dx_ex(physical_x, physical_y, heat_.t()) +
              edge_vector[1] *
                  data_.dc_dy_ex(physical_x, physical_y, heat_.t());

          // Approximate Dirichlet and gradient.
          Vector<T> g_bar = phi * g_modals(indices);
          Vector<T> grad_g_t_bar = grad_t * g_modals(indices);

          // Local estimator, R_{K, J}^2.
          local_estimate +=
              penalties[k] * dot(scaled, (uh - g_bar) * (uh - g_bar));

          // // Local estimator, R_{K, N}^2.
          // this->m_estimates[j] += sizes[j] * dot(scaled, grad_uh * grad_uh);

          // Local estimator, R_{K, T}^2.
          local_estimate +=
              sizes[j] * dot(scaled, (grad_uh_t - grad_g_t_bar) *
                                         (grad_uh_t - grad_g_t_bar));

          // Local data oscillation, O_{K, J}^2.
          local_estimate +=
              penalties[k] * dot(scaled, (g - g_bar) * (g - g_bar));

          // Local data oscillation, O_{K, T}^2.
          local_estimate +=
              sizes[j] * dot(scaled, (grad_g_t - grad_g_t_bar) *
                                         (grad_g_t - grad_g_t_bar));
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

          Matrix<T> n_grad =
              normal_vector[0] * n_gradx_phi + normal_vector[1] * n_grady_phi;
          Vector<T> n_grad_uh = n_grad * ch_(n_indices);

          Matrix<T> n_grad_t =
              edge_vector[0] * n_gradx_phi + edge_vector[1] * n_grady_phi;
          Vector<T> n_grad_uh_t = n_grad_t * ch_(n_indices);

          // Local estimator, R_{K, J}^2.
          local_estimate +=
              penalties[k] * dot(scaled, (uh - n_uh) * (uh - n_uh));

          // Local estimator, R_{K, N}^2.
          local_estimate += sizes[j] * dot(scaled, (grad_uh - n_grad_uh) *
                                                       (grad_uh - n_grad_uh));

          // Local estimator, R_{K, T}^2.
          local_estimate +=
              sizes[j] * dot(scaled, (grad_uh_t - n_grad_uh_t) *
                                         (grad_uh_t - n_grad_uh_t));
        }
      }

      this->m_estimates[j] = local_estimate;
      estimate += local_estimate;

      // Degrees.
      Vector<T> degrees{indices.size()};
      std::size_t counter = 0;

      for (std::size_t i = 0; i < mesh.elements[j].degree + 1; ++i)
        for (std::size_t k = 0; k < mesh.elements[j].degree + 1 - i; ++k) {
          degrees[counter] = static_cast<T>(i + k);
          ++counter;
        }

      // Coefficients.
      Vector<T> coefficients = ch_(indices);

      for (auto &coefficient : coefficients.elements)
        coefficient = std::log(std::abs(coefficient));

      // Fit.
      Vector<T> fit = this->polyfit(degrees, coefficients, 1);
      this->m_fits[j] = -fit[1];
    }

    this->m_estimate = std::sqrt(estimate);
  };
};

/**
 * @brief Fisher equation error estimator class.
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class FisherEstimator : public HeatEstimator<T> {
public:
  // CONSTRUCTOR.
  FisherEstimator(const Mesh &mesh_) : HeatEstimator<T>(mesh_) {};

  // METHODS.

  /**
   * @brief Compute Fisher-KPP error estimates.
   *
   * @param data_ Fisher equation data structure.
   * @param fisher_ Fisher equation object.
   * @param ch_ Numerical solution at current time step.
   */
  void computeEstimates(const DataFKPP &data_, const Fisher<T> &fisher_,
                        const Vector<T> &ch_) {
#ifndef NVERBOSE
    std::cout << "Evaluating estimates." << std::endl;
#endif

    // Number of elements.
    Mesh mesh = this->m_mesh;
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

    // Sizes.
    Vector<T> sizes{mesh.elements.size()};

    for (std::size_t j = 0; j < sizes.length; ++j) {
      Element element{mesh.elements[j]};

      for (const auto &p : element.element.points)
        for (const auto &q : element.element.points)
          sizes[j] = (distance(p, q) > sizes[j]) ? distance(p, q) : sizes[j];
    }

    // Coefficients.
    Vector<T> f_modals = fisher_.modal_source(data_, mesh);
    Vector<T> g_modals = fisher_.modal(mesh, data_.DirBC);

    // Estimate.
    T estimate = 0.0;

    // Loop over the elements.
#pragma omp parallel for reduction(+ : estimate) schedule(dynamic)
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

      // Local estimate.
      T local_estimate = 0.0;

      // Loop over the sub-triangulation.
      for (const auto &triangle : triangles) {

        // Jacobian's determinant and physical nodes.
        auto [jacobian_det, physical_x, physical_y] =
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Weights scaling.
        Vector<T> scaled = jacobian_det * weights_2d;

        // Basis functions.
        Matrix<T> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];
        Matrix<T> lap_phi = lap_basis_2d(mesh, j, {physical_x, physical_y});

        // Coefficients.
        Vector<T> D_ext = data_.D_ext(physical_x, physical_y, fisher_.t());
        Vector<T> alpha = data_.alpha(physical_x, physical_y, fisher_.t());

        // Local time derivative.
        Vector<T> local_uh = phi * ch_(indices);
        Vector<T> local_uh_old = phi * fisher_.ch_old()(indices);
        Vector<T> partial_uh_t = (local_uh - local_uh_old) / data_.dt;

        // Local numerical laplacian.
        Vector<T> lap_uh = lap_phi * ch_(indices);

        // Local numerical non-linear term.
        Vector<T> nl_uh = alpha * local_uh * (1.0 - local_uh);

        // Local exact source.
        Vector<T> f =
            data_.source_f(physical_x, physical_y, fisher_.t(), D_ext, alpha);

        // Local source approximation.
        Vector<T> f_bar = phi * f_modals(indices);

        // Local estimator, R_{K, E}^2.
        local_estimate +=
            sizes[j] * sizes[j] *
            dot(scaled, (f_bar - partial_uh_t + lap_uh + nl_uh) *
                            (f_bar - partial_uh_t + lap_uh + nl_uh));

        // Local data oscillation, O_{K, E}^2.
        local_estimate +=
            sizes[j] * sizes[j] * dot(scaled, (f - f_bar) * (f - f_bar));
      }

      // Element's neighbours.
      std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

      // Penalties.
      Vector<Real> penalties = penalty(mesh, j, data_.penalty_coeff);

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
        Vector<T> scaled = std::abs(edges[k]) * weights_1d;

        // Basis functions.
        auto [phi, gradx_phi, grady_phi] =
            basis_2d(mesh, j, {physical_x, physical_y});
        Matrix<T> scaled_phi{phi};

        for (std::size_t l = 0; l < scaled_phi.columns; ++l)
          scaled_phi.column(l, scaled_phi.column(l) * scaled);

        // Local numerical solution and gradients.
        Vector<T> uh = phi * ch_(indices);

        Matrix<T> grad =
            normal_vector[0] * gradx_phi + normal_vector[1] * grady_phi;
        Vector<T> grad_uh = grad * ch_(indices);

        Matrix<T> grad_t =
            edge_vector[0] * gradx_phi + edge_vector[1] * grady_phi;
        Vector<T> grad_uh_t = grad_t * ch_(indices);

        if (neighbour == -1) { // Boundary edge.

          // Local exact Dirichlet and gradient.
          Vector<T> g = data_.DirBC(physical_x, physical_y, fisher_.t());
          Vector<T> grad_g_t =
              edge_vector[0] *
                  data_.dc_dx_ex(physical_x, physical_y, fisher_.t()) +
              edge_vector[1] *
                  data_.dc_dy_ex(physical_x, physical_y, fisher_.t());

          // Approximate Dirichlet and gradient.
          Vector<T> g_bar = phi * g_modals(indices);
          Vector<T> grad_g_t_bar = grad_t * g_modals(indices);

          // Local estimator, R_{K, J}^2.
          local_estimate +=
              penalties[k] * dot(scaled, (uh - g_bar) * (uh - g_bar));

          // // Local estimator, R_{K, N}^2.
          // this->m_estimates[j] += sizes[j] * dot(scaled, grad_uh * grad_uh);

          // Local estimator, R_{K, T}^2.
          local_estimate +=
              sizes[j] * dot(scaled, (grad_uh_t - grad_g_t_bar) *
                                         (grad_uh_t - grad_g_t_bar));

          // Local data oscillation, O_{K, J}^2.
          local_estimate +=
              penalties[k] * dot(scaled, (g - g_bar) * (g - g_bar));

          // Local data oscillation, O_{K, T}^2.
          local_estimate +=
              sizes[j] * dot(scaled, (grad_g_t - grad_g_t_bar) *
                                         (grad_g_t - grad_g_t_bar));
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

          Matrix<T> n_grad =
              normal_vector[0] * n_gradx_phi + normal_vector[1] * n_grady_phi;
          Vector<T> n_grad_uh = n_grad * ch_(n_indices);

          Matrix<T> n_grad_t =
              edge_vector[0] * n_gradx_phi + edge_vector[1] * n_grady_phi;
          Vector<T> n_grad_uh_t = n_grad_t * ch_(n_indices);

          // Local estimator, R_{K, J}^2.
          local_estimate +=
              penalties[k] * dot(scaled, (uh - n_uh) * (uh - n_uh));

          // Local estimator, R_{K, N}^2.
          local_estimate += sizes[j] * dot(scaled, (grad_uh - n_grad_uh) *
                                                       (grad_uh - n_grad_uh));

          // Local estimator, R_{K, T}^2.
          local_estimate +=
              sizes[j] * dot(scaled, (grad_uh_t - n_grad_uh_t) *
                                         (grad_uh_t - n_grad_uh_t));
        }
      }

      this->m_estimates[j] = local_estimate;
      estimate += local_estimate;

      // Degrees.
      Vector<T> degrees{indices.size()};
      std::size_t counter = 0;

      for (std::size_t i = 0; i < mesh.elements[j].degree + 1; ++i)
        for (std::size_t k = 0; k < mesh.elements[j].degree + 1 - i; ++k) {
          degrees[counter] = static_cast<T>(i + k);
          ++counter;
        }

      // Coefficients.
      Vector<T> coefficients = ch_(indices);

      for (auto &coefficient : coefficients.elements)
        coefficient = std::log(std::abs(coefficient));

      // Fit.
      Vector<T> fit = this->polyfit(degrees, coefficients, 1);
      this->m_fits[j] = -fit[1];
    }

    this->m_estimate = std::sqrt(estimate);
  };
};

} // namespace pacs

#endif