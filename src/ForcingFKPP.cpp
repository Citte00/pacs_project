/**
 * @file ForcingFKPP.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

namespace pacs {

    /**
     * @brief Assembly of the forcing term.
     * 
     * @param data Data struct.
     * @param mesh Mesh struct.
     * @param t Time step.
     * @return Vector<Real> 
     */
    Vector<Real> forcingFKPP(const DataFKPP &data, const Mesh &mesh, const Real &t) {

        #ifndef NVERBOSE
        std::cout << "Computing the forcing term." << std::endl;
        #endif

        // Number of quadrature nodes.
        std::vector<std::size_t> nqn(mesh.elements.size(), 0);
        std::transform(mesh.elements.begin(), mesh.elements.end(), nqn.begin(),
                       [](const Element &elem) { return 2 * elem.degree + 1; });

        // Degrees of freedom.
        std::size_t dofs = mesh.dofs();

        // Starting indices.
        std::vector<std::size_t> starts;
        starts.emplace_back(0);

        for(std::size_t j = 1; j < mesh.elements.size(); ++j)
            starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

        // Neighbours.
        std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

        // Forcing term.
        Vector<Real> forcing{dofs};

        // Volume integrals.

        // Loop over the elements.
        for(std::size_t j = 0; j < mesh.elements.size(); ++j) {
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
            Real jacobian_det = jacobian(0, 0) * jacobian(1, 1) -
                                jacobian(0, 1) * jacobian(1, 0);

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
            Vector<Real> Dext = data.D_ext(physical_x, physical_y, t);
            Vector<Real> alpha = data.alpha(physical_x, physical_y, t);
            Vector<Real> local_source =
                data.source_f(physical_x, physical_y, t, Dext, alpha);

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
            for(std::size_t k = 0; k < element_neighbours.size(); ++k) {
              // 1D Quadrature nodes and weights.
              auto [nodes_1d, weights_1d] = quadrature_1d(nqn[j]);

              // Neighbour information.
              auto [edge, neighbour, n_edge] = element_neighbours[k];

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
                auto [phi, gradx_phi, grady_phi] = basis_2d(mesh, j, {physical_x, physical_y});

                // Param initialisation.
                Vector<Real> Dext = data.D_ext(physical_x, physical_y, t);

                // Local matrix assembly.
                Matrix<Real> scaled_gradx{gradx_phi};
                Matrix<Real> scaled_grady{grady_phi};
                Matrix<Real> scaled_phi{phi};

                for(std::size_t l = 0; l < scaled_gradx.columns; ++l) {
                    scaled_gradx.column(l, (Dext * scaled_gradx.column(l)) * scaled);
                    scaled_grady.column(l, (Dext * scaled_grady.column(l)) * scaled);
                    scaled_phi.column(l, scaled_phi.column(l) * scaled);
                }

                Matrix<Real> scaled_grad = normal_vector[0] * scaled_gradx + normal_vector[1] * scaled_grady;

                // Boundary conditions.
                Vector<Real> boundary = data.DirBC(physical_x, physical_y, t);

                // Local forcing term.
                local_f -= scaled_grad.transpose() * boundary;
                local_f += penalties[k] * scaled_phi.transpose() * boundary;
            }

            // Global forcing term.
            forcing(indices, local_f);

        }

        return forcing;

    }

}