/**
 * @file Estimators.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief
 * @date 2025-01-14
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <PacsHPDG.hpp>
#include <fstream>
#include <iomanip>

namespace pacs {

/**
 * @brief Polynomial fit.
 *
 * @param x X.
 * @param y Y.
 * @param p Algorithm order.
 * @return Vector<Real>
 */
Vector<Real> LaplaceEstimator::polyfit(const Vector<Real> &x,
                                       const Vector<Real> &y,
                                       const std::size_t &p) const {
#ifndef NDEBUG // Integrity check.
  assert(p > 0);
  assert(x.length == y.length);
#endif

  // X.
  Matrix<Real> X{x.length, p + 1};

  // Building X.
  for (std::size_t j = 0; j < p + 1; ++j) {
    X.column(j, Vector<Real>{x.length, 1.0});

    for (std::size_t k = 0; k < j; ++k) {
      X.column(j, X.column(j) * x);
    }
  }

  // Solution.
  return solve(X.transpose() * X, X.transpose() * y, QRD);
};

/**
 * @brief Refines specified elements' size for a mesh.
 *
 * @param mask Vector to identify which element to refine.
 */
void LaplaceEstimator::mesh_refine_size(const Mask &mask) {
#ifndef NDEBUG // Integrity check.
  assert(mask.size() == this->m_mesh.elements.size());
#endif

#ifndef NVERBOSE
  std::cout << "Refining mesh size: " << std::flush;
#endif

  // Degrees.
  std::vector<std::size_t> degrees;

  for (std::size_t j = 0; j < mask.size(); ++j)
    if (!(mask[j]))
      degrees.emplace_back(this->m_mesh.elements[j].degree);

  // Polygons.
  std::vector<Polygon> diagram;
  std::vector<Polygon> refine;
  std::vector<Polygon> refined;

  for (std::size_t j = 0; j < mask.size(); ++j)
    if (mask[j]) {
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
  this->m_mesh = new_mesh;
};

/**
 * @brief Refines specified elements' degree for a mesh.
 *
 * @param mask Vector to identify which element to refine.
 */
void LaplaceEstimator::mesh_refine_degree(const Mask &mask) {
#ifndef NDEBUG // Integrity check.
  assert(mask.size() == this->m_mesh.elements.size());
#endif

#ifndef NVERBOSE
  std::cout << "Refining mesh degree: " << std::flush;
  std::size_t counter = 0;
#endif

  for (std::size_t j = 0; j < mask.size(); ++j)
    if (mask[j]) {
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
 * @brief hp-Adaptively refine a mesh.
 *
 * @param estimator Error estimator.
 * @param refine Refinement percentage.
 * @param speed Solution's smoothness.
 */
void LaplaceEstimator::mesh_refine(const LaplaceEstimator &estimator,
                                   const Real &refine, const Real &speed) {
#ifndef NDEBUG // Integrity check.
  assert((refine > 0.0L) && (refine < 1.0L));
  assert(speed > 0.0L);
#endif

  // Masks.
  Mask p_mask = this->m_fits > speed;
  Mask h_mask = (this->m_estimates * this->m_estimates) >
                refine * sum(this->m_estimates * this->m_estimates) /
                    this->m_mesh.elements.size();

  // Strategy.
  for (std::size_t j = 0; j < this->m_mesh.elements.size(); ++j) {
    if (!h_mask[j]) // p-Refine only error-marked elements.
      p_mask[j] = false;

    if (p_mask[j] && h_mask[j]) // p > h.
      h_mask[j] = false;
  }

  // Refinements.
  mesh_refine_degree(p_mask);
  mesh_refine_size(h_mask);
};

/**
 * @brief Compute error estimates.
 *
 * @param data Laplace equation data struct.
 * @param laplacian Laplace equation object.
 * @param numerical Numerical solution.
 */
void LaplaceEstimator::computeEstimates(const DataLaplace &data,
                                        const Laplace &laplacian,
                                        const Vector<Real> &numerical) {
#ifndef NVERBOSE
  std::cout << "Evaluating estimates." << std::endl;
#endif

  // Quadrature nodes.
  std::vector<std::size_t> nqn(m_mesh.elements.size(), 0);
  std::transform(m_mesh.elements.begin(), m_mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.reserve(m_mesh.elements.size());
  starts.emplace_back(0);

  for (std::size_t j = 1; j < m_mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + m_mesh.elements[j - 1].dofs());

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = m_mesh.neighbours;

  // Sizes.
  Vector<Real> sizes{m_mesh.elements.size()};

  for (std::size_t j = 0; j < sizes.length; ++j) {
    Element element{m_mesh.elements[j]};

    for (const auto &p : element.element.points)
      for (const auto &q : element.element.points)
        sizes[j] = (distance(p, q) > sizes[j]) ? distance(p, q) : sizes[j];
  }

  // Mass blocks.
  auto blocks = laplacian.block_mass(m_mesh);

  // Coefficients.
  Vector<Real> f_modals = laplacian.modal(m_mesh, data.source_f);
  Vector<Real> g_modals = laplacian.modal(m_mesh, data.DirBC);

  // Loop over the elements.
  for (std::size_t j = 0; j < m_mesh.elements.size(); ++j) {

    // 2D quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

    // Local dofs.
    std::size_t element_dofs = m_mesh.elements[j].dofs();

    // Global matrix indices.
    std::vector<std::size_t> indices;
    indices.reserve(element_dofs);

    for (std::size_t k = 0; k < element_dofs; ++k)
      indices.emplace_back(starts[j] + k);

    // Polygon.
    Polygon polygon = m_mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Loop over the sub-triangulation.
    for (std::size_t k = 0; k < triangles.size(); ++k) {

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangles[k], {nodes_x_2d, nodes_y_2d});

      // Weights scaling.
      Vector<Real> scaled = jacobian_det * weights_2d;

      // Basis functions.
      Matrix<Real> phi = basis_2d(m_mesh, j, {physical_x, physical_y})[0];
      Matrix<Real> lap_phi = lap_basis_2d(m_mesh, j, {physical_x, physical_y});

      // Local numerical laplacian.
      Vector<Real> lap_uh = lap_phi * numerical(indices);

      // Local exact source.
      Vector<Real> f = data.source_f(physical_x, physical_y);

      // Local source approximation.
      Vector<Real> f_bar = phi * f_modals(indices);

      // Local estimator, R_{K, E}^2.
      this->m_estimates[j] += sizes[j] * sizes[j] *
                              dot(scaled, (f_bar + lap_uh) * (f_bar + lap_uh));

      // Local data oscillation, O_{K, E}^2.
      this->m_estimates[j] +=
          sizes[j] * sizes[j] * dot(scaled, (f - f_bar) * (f - f_bar));
    }

    // Element's neighbours.
    std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

    // Penalties.
    Vector<Real> penalties = penalty(m_mesh, j, data.penalty_coeff);

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
          basis_2d(m_mesh, j, {physical_x, physical_y});
      Matrix<Real> scaled_phi{phi};

      for (std::size_t l = 0; l < scaled_phi.columns; ++l)
        scaled_phi.column(l, scaled_phi.column(l) * scaled);

      // Local numerical solution and gradients.
      Vector<Real> uh = phi * numerical(indices);

      Matrix<Real> grad =
          normal_vector[0] * gradx_phi + normal_vector[1] * grady_phi;
      Vector<Real> grad_uh = grad * numerical(indices);

      Matrix<Real> grad_t =
          edge_vector[0] * gradx_phi + edge_vector[1] * grady_phi;
      Vector<Real> grad_uh_t = grad_t * numerical(indices);

      if (neighbour == -1) { // Boundary edge.

        // Local exact Dirichlet and gradient.
        Vector<Real> g = data.DirBC(physical_x, physical_y);
        Vector<Real> grad_g_t =
            edge_vector[0] * data.dc_dx_ex(physical_x, physical_y) +
            edge_vector[1] * data.dc_dy_ex(physical_x, physical_y);

        // Approximate Dirichlet and gradient.
        Vector<Real> g_bar = phi * g_modals(indices);
        Vector<Real> grad_g_t_bar = grad_t * g_modals(indices);

        // Local estimator, R_{K, J}^2.
        this->m_estimates[j] +=
            penalties[k] * dot(scaled, (uh - g_bar) * (uh - g_bar));

        // // Local estimator, R_{K, N}^2.
        // this->m_estimates[j] += sizes[j] * dot(scaled, grad_uh * grad_uh);

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
            basis_2d(m_mesh, neighbour, {physical_x, physical_y});

        std::vector<std::size_t> n_indices;
        std::size_t n_index = element_neighbours[k][1];
        std::size_t n_dofs =
            m_mesh.elements[n_index].dofs(); // Neighbour's dofs.

        n_indices.reserve(n_dofs);
        for (std::size_t h = 0; h < n_dofs; ++h)
          n_indices.emplace_back(starts[n_index] + h);

        // Neighbour's numerical solution and gradients.
        Vector<Real> n_uh = n_phi * numerical(n_indices);

        Matrix<Real> n_grad =
            normal_vector[0] * n_gradx_phi + normal_vector[1] * n_grady_phi;
        Vector<Real> n_grad_uh = n_grad * numerical(n_indices);

        Matrix<Real> n_grad_t =
            edge_vector[0] * n_gradx_phi + edge_vector[1] * n_grady_phi;
        Vector<Real> n_grad_uh_t = n_grad_t * numerical(n_indices);

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

    for (std::size_t i = 0; i < m_mesh.elements[j].degree + 1; ++i)
      for (std::size_t k = 0; k < m_mesh.elements[j].degree + 1 - i; ++k) {
        degrees[counter] = static_cast<Real>(i + k);
        ++counter;
      }

    // Coefficients.
    Vector<Real> coefficients = numerical(indices);

    for (auto &coefficient : coefficients.elements)
      coefficient = std::log(std::abs(coefficient));

    // Fit.
    Vector<Real> fit = polyfit(degrees, coefficients, 1);
    this->m_fits[j] = -fit[1];
  }

  this->m_estimate = std::sqrt(this->m_estimate);
};

/**
 * @brief Outputs the mesh to a polyplot.py readable file.
 *
 * @param filename File name.
 * @param estimates Bool value.
 */
void LaplaceEstimator::write(const std::string &filename, const bool &estimates) {
  // File loading.
  std::ofstream file{filename};

  file << "@ " << filename << "\n";
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
      file << std::setprecision(12) << vertex[0] << " " << std::setprecision(12)
           << vertex[1] << " ";

    if (estimates)
      file << this->m_estimates[count];

    file << "\n";
    count++;
  }
};

/**
 * @brief Compute heat equation error estimates.
 *
 * @param data Heat equation data struct.
 * @param heat Heat equation object.
 * @param ch_old Numerical solution at the previous time step.
 */
void HeatEstimator::computeEstimates(const DataHeat &data, const Heat &heat,
                                     const Vector<Real> &ch_old) {
#ifndef NVERBOSE
  std::cout << "Evaluating estimates." << std::endl;
#endif

  // Quadrature nodes.
  std::vector<std::size_t> nqn(m_mesh.elements.size(), 0);
  std::transform(m_mesh.elements.begin(), m_mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.reserve(m_mesh.elements.size());
  starts.emplace_back(0);

  for (std::size_t j = 1; j < m_mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + m_mesh.elements[j - 1].dofs());

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = m_mesh.neighbours;

  // Sizes.
  Vector<Real> sizes{m_mesh.elements.size()};

  for (std::size_t j = 0; j < sizes.length; ++j) {
    Element element{m_mesh.elements[j]};

    for (const auto &p : element.element.points)
      for (const auto &q : element.element.points)
        sizes[j] = (distance(p, q) > sizes[j]) ? distance(p, q) : sizes[j];
  }

  // Coefficients.
  Vector<Real> f_modals = heat.modal_source(data, m_mesh);
  Vector<Real> g_modals = heat.modal(m_mesh, data.DirBC);

  // Loop over the elements.
  for (std::size_t j = 0; j < m_mesh.elements.size(); ++j) {

    // 2D quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

    // Local dofs.
    std::size_t element_dofs = m_mesh.elements[j].dofs();

    // Global matrix indices.
    std::vector<std::size_t> indices;
    indices.reserve(element_dofs);

    for (std::size_t k = 0; k < element_dofs; ++k)
      indices.emplace_back(starts[j] + k);

    // Polygon.
    Polygon polygon = m_mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Loop over the sub-triangulation.
    for (std::size_t k = 0; k < triangles.size(); ++k) {

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangles[k], {nodes_x_2d, nodes_y_2d});

      // Weights scaling.
      Vector<Real> scaled = jacobian_det * weights_2d;

      // Basis functions.
      Matrix<Real> phi = basis_2d(m_mesh, j, {physical_x, physical_y})[0];
      Matrix<Real> lap_phi = lap_basis_2d(m_mesh, j, {physical_x, physical_y});

      // Local exact source.
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y, heat.t());
      Vector<Real> f = data.source_f(physical_x, physical_y, heat.t(), D_ext);

      // Local numerical heat.
      Vector<Real> lap_uh = D_ext * (lap_phi * heat.ch()(indices));

      // Local time derivative.
      Vector<Real> local_uh = phi * heat.ch()(indices);
      Vector<Real> local_uh_old = phi * ch_old(indices);
      Vector<Real> partial_uh_t = (local_uh - local_uh_old) / data.dt;

      // Local source approximation.
      Vector<Real> f_bar = phi * f_modals(indices);

      // Local estimator, R_{K, E}^2.
      this->m_estimates[j] += sizes[j] * sizes[j] *
                              dot(scaled, (f_bar - partial_uh_t + lap_uh) *
                                              (f_bar - partial_uh_t + lap_uh));

      // Local data oscillation, O_{K, E}^2.
      this->m_estimates[j] +=
          sizes[j] * sizes[j] * dot(scaled, (f - f_bar) * (f - f_bar));
    }

    // Element's neighbours.
    std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

    // Penalties.
    Vector<Real> penalties = penalty(m_mesh, j, data.penalty_coeff);

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
          basis_2d(m_mesh, j, {physical_x, physical_y});
      Matrix<Real> scaled_phi{phi};

      for (std::size_t l = 0; l < scaled_phi.columns; ++l)
        scaled_phi.column(l, scaled_phi.column(l) * scaled);

      // Local numerical solution and gradients.
      Vector<Real> uh = phi * heat.ch()(indices);

      Matrix<Real> grad =
          normal_vector[0] * gradx_phi + normal_vector[1] * grady_phi;
      Vector<Real> grad_uh = grad * heat.ch()(indices);

      Matrix<Real> grad_t =
          edge_vector[0] * gradx_phi + edge_vector[1] * grady_phi;
      Vector<Real> grad_uh_t = grad_t * heat.ch()(indices);

      if (neighbour == -1) { // Boundary edge.

        // Local exact Dirichlet and gradient.
        Vector<Real> g = data.DirBC(physical_x, physical_y, heat.t());
        Vector<Real> grad_g_t =
            edge_vector[0] * data.dc_dx_ex(physical_x, physical_y, heat.t()) +
            edge_vector[1] * data.dc_dy_ex(physical_x, physical_y, heat.t());

        // Approximate Dirichlet and gradient.
        Vector<Real> g_bar = phi * g_modals(indices);
        Vector<Real> grad_g_t_bar = grad_t * g_modals(indices);

        // Local estimator, R_{K, J}^2.
        this->m_estimates[j] +=
            penalties[k] * dot(scaled, (uh - g_bar) * (uh - g_bar));

        // // Local estimator, R_{K, N}^2.
        // this->m_estimates[j] += sizes[j] * dot(scaled, grad_uh * grad_uh);

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
            basis_2d(m_mesh, neighbour, {physical_x, physical_y});

        std::vector<std::size_t> n_indices;
        std::size_t n_index = element_neighbours[k][1];
        std::size_t n_dofs =
            m_mesh.elements[n_index].dofs(); // Neighbour's dofs.

        n_indices.reserve(n_dofs);
        for (std::size_t h = 0; h < n_dofs; ++h)
          n_indices.emplace_back(starts[n_index] + h);

        // Neighbour's numerical solution and gradients.
        Vector<Real> n_uh = n_phi * heat.ch()(n_indices);

        Matrix<Real> n_grad =
            normal_vector[0] * n_gradx_phi + normal_vector[1] * n_grady_phi;
        Vector<Real> n_grad_uh = n_grad * heat.ch()(n_indices);

        Matrix<Real> n_grad_t =
            edge_vector[0] * n_gradx_phi + edge_vector[1] * n_grady_phi;
        Vector<Real> n_grad_uh_t = n_grad_t * heat.ch()(n_indices);

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

    for (std::size_t i = 0; i < m_mesh.elements[j].degree + 1; ++i)
      for (std::size_t k = 0; k < m_mesh.elements[j].degree + 1 - i; ++k) {
        degrees[counter] = static_cast<Real>(i + k);
        ++counter;
      }

    // Coefficients.
    Vector<Real> coefficients = heat.ch()(indices);

    for (auto &coefficient : coefficients.elements)
      coefficient = std::log(std::abs(coefficient));

    // Fit.
    Vector<Real> fit = polyfit(degrees, coefficients, 1);
    this->m_fits[j] = -fit[1];
  }

  this->m_estimate = std::sqrt(this->m_estimate);
};

void HeatEstimator::mesh_refine(Heat &heat, const Mesh &mesh, const Real &refine,
                                const Real &speed) {
#ifndef NDEBUG // Integrity check.
  assert((refine > 0.0L) && (refine < 1.0L));
  assert(speed > 0.0L);
#endif

  // Mesh
  Mesh old_mesh = mesh;

  // Masks.
  Mask p_mask = this->m_fits > speed;
  Mask h_mask = (this->m_estimates * this->m_estimates) >
                refine * sum(this->m_estimates * this->m_estimates) /
                    this->m_mesh.elements.size();

  // Strategy.
  for (std::size_t j = 0; j < this->m_mesh.elements.size(); ++j) {
    if (!h_mask[j]) // p-Refine only error-marked elements.
      p_mask[j] = false;

    if (p_mask[j] && h_mask[j]) // p > h.
      h_mask[j] = false;
  }

  // Refinements and update of solution.
  mesh_refine_degree(p_mask);
  heat.prolong_solution_p(this->m_mesh, old_mesh, p_mask);
  old_mesh = this->m_mesh;
  mesh_refine_size(h_mask);
  heat.prolong_solution_h(this->m_mesh, old_mesh, h_mask);
};

/**
 * @brief Compute fisher equation error estimates.
 *
 * @param data Fisher equation data struct.
 * @param fisher Fisher equation object.
 */
void FisherEstimator::computeEstimates(const DataFKPP &data,
                                       const Fisher &fisher) {
#ifndef NVERBOSE
  std::cout << "Evaluating estimates." << std::endl;
#endif

  // Quadrature nodes.
  std::vector<std::size_t> nqn(m_mesh.elements.size(), 0);
  std::transform(m_mesh.elements.begin(), m_mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.reserve(m_mesh.elements.size());
  starts.emplace_back(0);

  for (std::size_t j = 1; j < m_mesh.elements.size(); ++j)
    starts.emplace_back(starts[j - 1] + m_mesh.elements[j - 1].dofs());

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours = m_mesh.neighbours;

  // Sizes.
  Vector<Real> sizes{m_mesh.elements.size()};

  for (std::size_t j = 0; j < sizes.length; ++j) {
    Element element{m_mesh.elements[j]};

    for (const auto &p : element.element.points)
      for (const auto &q : element.element.points)
        sizes[j] = (distance(p, q) > sizes[j]) ? distance(p, q) : sizes[j];
  }

  // Coefficients.
  Vector<Real> f_modals = fisher.modal_source(data, m_mesh);
  Vector<Real> g_modals = fisher.modal(m_mesh, data.DirBC);

  // Loop over the elements.
  for (std::size_t j = 0; j < m_mesh.elements.size(); ++j) {

    // 2D quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

    // Local dofs.
    std::size_t element_dofs = m_mesh.elements[j].dofs();

    // Global matrix indices.
    std::vector<std::size_t> indices;
    indices.reserve(element_dofs);

    for (std::size_t k = 0; k < element_dofs; ++k)
      indices.emplace_back(starts[j] + k);

    // Polygon.
    Polygon polygon = m_mesh.element(j);

    // Element sub-triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);

    // Loop over the sub-triangulation.
    for (std::size_t k = 0; k < triangles.size(); ++k) {

      // Jacobian's determinant and physical nodes.
      auto [jacobian_det, physical_x, physical_y] =
          get_Jacobian_physical_points(triangles[k], {nodes_x_2d, nodes_y_2d});

      // Weights scaling.
      Vector<Real> scaled = jacobian_det * weights_2d;

      // Basis functions.
      Matrix<Real> phi = basis_2d(m_mesh, j, {physical_x, physical_y})[0];
      Matrix<Real> lap_phi = lap_basis_2d(m_mesh, j, {physical_x, physical_y});

      // Coefficients.
      Vector<Real> D_ext = data.D_ext(physical_x, physical_y, fisher.t());
      Vector<Real> alpha = data.alpha(physical_x, physical_y, fisher.t());

      // Local time derivative.
      Vector<Real> local_uh = phi * fisher.ch()(indices);
      Vector<Real> local_uh_old = phi * fisher.ch_old()(indices);
      Vector<Real> partial_uh_t = (local_uh - local_uh_old) / data.dt;

      // Local numerical laplacian.
      Vector<Real> lap_uh = lap_phi * fisher.ch()(indices);

      // Local numerical non-linear term.
      Vector<Real> nl_uh = alpha * local_uh * (1.0 - local_uh);

      // Local exact source.
      Vector<Real> f =
          data.source_f(physical_x, physical_y, fisher.t(), D_ext, alpha);

      // Local source approximation.
      Vector<Real> f_bar = phi * f_modals(indices);

      // Local estimator, R_{K, E}^2.
      this->m_estimates[j] +=
          sizes[j] * sizes[j] *
          dot(scaled, (f_bar - partial_uh_t + lap_uh + nl_uh) *
                          (f_bar - partial_uh_t + lap_uh + nl_uh));

      // Local data oscillation, O_{K, E}^2.
      this->m_estimates[j] +=
          sizes[j] * sizes[j] * dot(scaled, (f - f_bar) * (f - f_bar));
    }

    // Element's neighbours.
    std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

    // Penalties.
    Vector<Real> penalties = penalty(m_mesh, j, data.penalty_coeff);

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
          basis_2d(m_mesh, j, {physical_x, physical_y});
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

        // // Local estimator, R_{K, N}^2.
        // this->m_estimates[j] += sizes[j] * dot(scaled, grad_uh * grad_uh);

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
            basis_2d(m_mesh, neighbour, {physical_x, physical_y});

        std::vector<std::size_t> n_indices;
        std::size_t n_index = element_neighbours[k][1];
        std::size_t n_dofs =
            m_mesh.elements[n_index].dofs(); // Neighbour's dofs.

        n_indices.reserve(n_dofs);
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

    for (std::size_t i = 0; i < m_mesh.elements[j].degree + 1; ++i)
      for (std::size_t k = 0; k < m_mesh.elements[j].degree + 1 - i; ++k) {
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