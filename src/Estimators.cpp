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
                                       const std::size_t &p) {
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
 * @brief Compute error estimates.
 *
 * @param data Laplace equation data struct.
 * @param mesh Mesh struct.
 * @param laplacian Laplace equation object.
 * @param numerical Numerical solution.
 */
void LaplaceEstimator::computeEstimates(const DataLaplace &data,
                                        const Mesh &mesh,
                                        const Laplace &laplacian,
                                        const Vector<Real> &numerical) {
#ifndef NVERBOSE
  std::cout << "Evaluating estimates." << std::endl;
#endif

  // Quadrature nodes.
  std::vector<std::size_t> nqn(mesh.elements.size(), 0);
  std::transform(mesh.elements.begin(), mesh.elements.end(), nqn.begin(),
                 [](const Element &elem) { return 2 * elem.degree + 1; });

  // Starting indices.
  std::vector<std::size_t> starts;
  starts.reserve(mesh.elements.size());
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

  // Mass blocks.
  auto blocks = laplacian.block_mass(mesh);

  // Coefficients.
  Vector<Real> f_modals = modal(mesh, data.source_f);
  Vector<Real> g_modals = modal(mesh, data.DirBC);

  // Loop over the elements.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {

    // 2D quadrature nodes and weights.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn[j]);

    // Local dofs.
    std::size_t element_dofs = mesh.elements[j].dofs();

    // Global matrix indices.
    std::vector<std::size_t> indices;
    indices.reserve(element_dofs);

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
      Matrix<Real> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];
      Matrix<Real> lap_phi = lap_basis_2d(mesh, j, {physical_x, physical_y});

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
          basis_2d(mesh, j, {physical_x, physical_y});
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
            basis_2d(mesh, neighbour, {physical_x, physical_y});

        std::vector<std::size_t> n_indices;
        std::size_t n_index = element_neighbours[k][1];
        std::size_t n_dofs = mesh.elements[n_index].dofs(); // Neighbour's dofs.

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

    for (std::size_t i = 0; i < mesh.elements[j].degree + 1; ++i)
      for (std::size_t k = 0; k < mesh.elements[j].degree + 1 - i; ++k) {
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
 * @brief Refines specified elements' size for a mesh.
 *
 * @param mesh Mesh.
 * @param mask Mask.
 */
void LaplaceEstimator::mesh_refine_size(Mesh &mesh, const Mask &mask) {
#ifndef NDEBUG // Integrity check.
  assert(mask.size() == mesh.elements.size());
#endif

#ifndef NVERBOSE
  std::cout << "Refining mesh size: " << std::flush;
#endif

  // Degrees.
  std::vector<std::size_t> degrees;

  for (std::size_t j = 0; j < mask.size(); ++j)
    if (!(mask[j]))
      degrees.emplace_back(mesh.elements[j].degree);

  // Polygons.
  std::vector<Polygon> diagram;
  std::vector<Polygon> refine;
  std::vector<Polygon> refined;

  for (std::size_t j = 0; j < mask.size(); ++j)
    if (mask[j]) {
      refine.emplace_back(mesh.element(j));

      for (std::size_t k = 0;
           k < mesh.elements[j].edges.size() +
                   static_cast<std::size_t>(mesh.elements[j].edges.size() > 4);
           ++k)
        degrees.emplace_back(mesh.elements[j].degree);

    } else
      diagram.emplace_back(mesh.element(j));

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
  mesh = Mesh{mesh.domain, diagram, degrees};
};

/**
 * @brief Refines specified elements' degree for a mesh.
 *
 * @param mesh Mesh.
 * @param mask Mask.
 */
void LaplaceEstimator::mesh_refine_degree(Mesh &mesh, const Mask &mask) {
#ifndef NDEBUG // Integrity check.
  assert(mask.size() == mesh.elements.size());
#endif

#ifndef NVERBOSE
  std::cout << "Refining mesh degree: " << std::flush;
  std::size_t counter = 0;
#endif

  for (std::size_t j = 0; j < mask.size(); ++j)
    if (mask[j]) {
      ++mesh.elements[j].degree;

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
 * @param mesh Mesh.
 * @param estimator Estimator.
 * @param refine Refinement percentage.
 * @param speed Solution's smoothness.
 */
void LaplaceEstimator::mesh_refine(Mesh &mesh, const LaplaceEstimator &estimator,
                 const Real &refine, const Real &speed) {
#ifndef NDEBUG // Integrity check.
  assert((refine > 0.0L) && (refine < 1.0L));
  assert(speed > 0.0L);
#endif

  // Masks.
  Mask p_mask = this->m_fits > speed;
  Mask h_mask = (this->m_estimates * this->m_estimates) >
                refine * sum(this->m_estimates * this->m_estimates) /
                    mesh.elements.size();

  // Strategy.
  for (std::size_t j = 0; j < mesh.elements.size(); ++j) {
    if (!h_mask[j]) // p-Refine only error-marked elements.
      p_mask[j] = false;

    if (p_mask[j] && h_mask[j]) // p > h.
      h_mask[j] = false;
  }

  // Refinements.
  mesh_refine_degree(mesh, p_mask);
  mesh_refine_size(mesh, h_mask);
};

} // namespace pacs