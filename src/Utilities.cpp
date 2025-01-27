/**
 * @file Penalty.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

namespace pacs {

    /**
     * @brief Returns the penalty coefficients for a given element. "Max" policy.
     * 
     * @param mesh Mesh.
     * @param index Element's index.
     * @return Vector<Real> 
     */
    Vector<Real> penalty(const Mesh &mesh, const std::size_t &index, const Real &coefficient) {
        
        // Element.
        Element element = mesh.elements[index];
        Polygon polygon = mesh.element(index);

        std::vector<std::array<int, 3>> neighbours = mesh.neighbours[index];

        // Sizes.
        std::vector<Segment> edges = polygon.edges();
        Vector<Real> sizes{element.edges.size()};

        for(std::size_t j = 0; j < sizes.length; ++j)
            sizes[j] = std::abs(edges[j]);

        // Element's area.
        Real area = mesh.areas[index];

        // Biggest simplices areas.
        Vector<Real> areas = mesh.max_simplices[index];

        // Inverse constant.
        Vector<Real> inverse = area / areas;

        // Coefficients.
        Real penalty_coefficient = coefficient * (element.degree * element.degree);
        Vector<Real> penalty_dirichlet = penalty_coefficient * inverse * sizes / area;

        // Penalty evaluation.
        Vector<Real> penalties{neighbours.size()};
        Vector<Real> internal{neighbours.size()}; // Element.
        Vector<Real> external{neighbours.size()}; // Neighbour.
        Vector<Real> inverse_external{neighbours.size()};

        for(std::size_t j = 0; j < neighbours.size(); ++j) {
            if(neighbours[j][1] == -1) {
                penalties[j] = penalty_dirichlet[j];
                continue;
            }

            inverse_external[j] = mesh.areas[neighbours[j][1]] / mesh.max_simplices[neighbours[j][1]][neighbours[j][2]];
            internal[j] = penalty_coefficient * inverse[j] * sizes[j] / area;
            external[j] = penalty_coefficient * inverse_external[j] * sizes[j] / mesh.areas[neighbours[j][1]];

            penalties[j] = (internal[j] > external[j]) ? internal[j] : external[j];
        }

        return penalties;
    };

    /**
     * @brief Get the Jacobian physical points object.
     * 
     * @param triangle Triangle.
     * @param points Physical points in R^2
     * @return std::tuple<Real, Vector<Real>, Vector<Real>> 
     */
    std::tuple<Real, Vector<Real>, Vector<Real>>
    get_Jacobian_physical_points(const Polygon &triangle,
                                 const std::array<Vector<Real>, 2> &points) {

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
      auto [nodes_x_2d, nodes_y_2d] = points;
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

      return std::make_tuple(jacobian_det, physical_x, physical_y);
    };

    std::array<Vector<Real>, 4>
    faces_physical_points(const Segment &segment, const Vector<Real> &nodes_1d) {

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

      return {normal_vector, edge_vector, physical_x, physical_y};
    };
}