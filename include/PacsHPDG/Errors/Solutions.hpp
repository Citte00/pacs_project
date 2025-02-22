/**
 * @file LaplaceSolution.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Equations readable and plottable solution class object.
 * @date 2025-01-12
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_ERRORS_SOLUTIONS_HPP
#define INCLUDE_PACSHPDG_ERRORS_SOLUTIONS_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Geometry.hpp"
#include "../Solvers.hpp"

#include <string>

namespace pacs {

/**
 * @brief Laplace equation solution class.
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class LaplaceSolution {
protected:
  Vector<T> m_x;
  Vector<T> m_y;
  Vector<T> m_numerical;
  Vector<T> m_exact;

public:
  // CONSTRUCTOR.
  LaplaceSolution(const Mesh &mesh_)
      : m_x{mesh_.entries}, m_y{mesh_.entries}, m_numerical{mesh_.entries},
        m_exact{mesh_.entries} {};

  // GETTERS.
  Vector<T> x() const { return this->m_x; };
  Vector<T> y() const { return this->m_y; };
  Vector<T> numerical() const { return this->m_numerical; };
  Vector<T> exact() const { return this->m_exact; };

  // METHODS.
  // Save numerical and exact solution for plotting.
  /**
   * @brief Store Laplace equation numerical and exact solution.
   *
   * @param data_ Laplace equation data structure.
   * @param mesh_ Mesh structure.
   * @param ch_ Numerical solution.
   */
  void solution(const DataLaplace &data_, const Mesh &mesh_,
                const Vector<T> &ch_) {
    // Number of quadrature nodes.
    std::size_t nqn = data_.NqnVisualization;

    // Quadrature nodes.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn);

    // Starting indices.
    std::vector<std::size_t> starts(mesh_.elements.size());
    starts[0] = 0;

    for (std::size_t j = 1; j < mesh_.elements.size(); ++j)
      starts[j] = starts[j - 1] + mesh_.elements[j - 1].dofs();

    // Local vectors indices.
    std::vector<std::size_t> local_indices(nqn * nqn);

    for (std::size_t h = 0; h < nqn * nqn; ++h)
      local_indices[h] = h;

    // Loop over the elements.
    for (std::size_t j = 0; j < mesh_.elements.size(); ++j) {

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
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Basis functions.
        Matrix<T> phi = basis_2d(mesh_, j, {physical_x, physical_y})[0];

        // Numerical solution.
        Vector<T> local_numerical = phi * ch_(indices);

        // Exact solution.
        Vector<T> local_exact = data_.c_ex(physical_x, physical_y);

        // Writing.
        this->m_x(local_indices, physical_x);
        this->m_y(local_indices, physical_y);
        this->m_numerical(local_indices, local_numerical);
        this->m_exact(local_indices, local_exact);

        // Local indices update.
        for (auto &index : local_indices)
          index += nqn * nqn;
      }
    }
  };

  /**
   * @brief Output solution in a .sol file.
   *
   * @param filename_ File name.
   */
  void write(const std::string &filename_) {
    // File loading.
    std::ofstream file{filename_};

    file << "@ " << filename_ << "\n";
    file << "@ contplot.py readable mesh\n";
    file << "@ Structure: [x, y, numerical(x, y), exact(x, y)].\n";

    for (std::size_t j = 0; j < this->m_x.length; ++j) {
      file << this->m_x[j] << ",";
      file << this->m_y[j] << ",";
      file << this->m_numerical[j] << ",";
      file << this->m_exact[j] << "\n";
    }

    file.close();
  };
};

/**
 * @brief
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class HeatSolution : public LaplaceSolution<T> {
public:
  // CONSTRUCTOR.
  HeatSolution(const Mesh &mesh_) : LaplaceSolution<T>(mesh_) {};

  // METHODS.
  /**
   * @brief  Store heat equation numerical and exact solution.
   *
   * @param data_ Heat equation data structure.
   * @param mesh_ Mesh structure.
   * @param heat_ Heat equation object.
   * @param ch_ Numerical solution.
   */
  void solution(const DataHeat &data_, const Mesh &mesh_, const Heat<T> &heat_,
                const Vector<T> &ch_) {
#ifndef NVERBOSE
    std::cout << "Evaluating solutions." << std::endl;
#endif

    // Number of quadrature nodes.
    std::size_t nqn = GAUSS_ORDER;

    // Quadrature nodes.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn);

    // Starting indices.
    std::vector<std::size_t> starts(mesh_.elements.size());
    starts[0] = 0;

    for (std::size_t j = 1; j < mesh_.elements.size(); ++j)
      starts[j] = starts[j - 1] + mesh_.elements[j - 1].dofs();

    // Local vectors indices.
    std::vector<std::size_t> local_indices(nqn * nqn);

    for (std::size_t h = 0; h < nqn * nqn; ++h)
      local_indices[h] = h;

    // Loop over the elements.
    for (std::size_t j = 0; j < mesh_.elements.size(); ++j) {

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
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Basis functions.
        Matrix<T> phi = basis_2d(mesh_, j, {physical_x, physical_y})[0];

        // Numerical solution.
        Vector<T> local_numerical = phi * ch_(indices);

        // Exact solution.
        Vector<T> local_exact = data_.c_ex(physical_x, physical_y, heat_.t());

        // Writing.
        this->m_x(local_indices, physical_x);
        this->m_y(local_indices, physical_y);
        this->m_numerical(local_indices, local_numerical);
        this->m_exact(local_indices, local_exact);

        // Local indices update.
        for (auto &index : local_indices)
          index += nqn * nqn;
      }
    }
  };
};

/**
 * @brief
 *
 * @tparam T Template parameter for numerical variable.
 */
template <NumericType T> class FisherSolution : public HeatSolution<T> {
public:
  // CONSTRUCTOR.
  FisherSolution(const Mesh &mesh_) : HeatSolution<T>(mesh_) {};

  // METHODS.
  /**
   * @brief  Store Fisher-KPP equation numerical and exact solution.
   *
   * @param data_ Fisher-KPP equation data structure.
   * @param mesh_ Mesh structure.
   * @param fisher_ Fisher-KPP equation object.
   * @param ch_ Numerical solution.
   */
  void solution(const DataFKPP &data_, const Mesh &mesh_,
                const Fisher<T> &fisher_, const Vector<T> &ch_) {
    // Number of quadrature nodes.
    std::size_t nqn = data_.NqnVisualization;

    // Quadrature nodes.
    auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn);

    // Starting indices.
    std::vector<std::size_t> starts(mesh_.elements.size());
    starts[0] = 0;

    for (std::size_t j = 1; j < mesh_.elements.size(); ++j)
      starts[j] = starts[j - 1] + mesh_.elements[j - 1].dofs();

    // Local vectors indices.
    std::vector<std::size_t> local_indices(nqn * nqn);

    for (std::size_t h = 0; h < nqn * nqn; ++h)
      local_indices[h] = h;

    // Loop over the elements.
    for (std::size_t j = 0; j < mesh_.elements.size(); ++j) {

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
            get_Jacobian_physical_points(triangle, {nodes_x_2d, nodes_y_2d});

        // Basis functions.
        Matrix<T> phi = basis_2d(mesh_, j, {physical_x, physical_y})[0];

        // Numerical solution.
        Vector<T> local_numerical = phi * ch_(indices);

        // Exact solution.
        Vector<T> local_exact = data_.c_ex(physical_x, physical_y, fisher_.t());

        // Writing.
        this->m_x(local_indices, physical_x);
        this->m_y(local_indices, physical_y);
        this->m_numerical(local_indices, local_numerical);
        this->m_exact(local_indices, local_exact);

        // Local indices update.
        for (auto &index : local_indices)
          index += nqn * nqn;
      }
    }
  };
};

} // namespace pacs

#endif