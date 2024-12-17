/**
 * @file EvaluateCoefficients.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

namespace pacs {

    /**
     * @brief Evaluate the modal coefficients of a function.
     * 
     * @param mesh 
     * @param function Function to evaluate.
     * @param t Time step.
     * @return Vector<Real> 
     */
    Vector<Real> evaluateCoeff(const Mesh &mesh, const TriFunctor &function, const Real &t) {

        if (!function.getFunction())
            throw std::runtime_error("evaluateCoeff received an uninitialized function");

        // Number of quadrature nodes.
        std::size_t degree = GAUSS_ORDER;

        // Coefficients.
        Vector<Real> coefficients{mesh.dofs()};

        // Quadrature nodes.
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(degree);

        // Starting indices.
        std::vector<std::size_t> starts;
        starts.emplace_back(0);

        for(std::size_t j = 1; j < mesh.elements.size(); ++j)
            starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

        // Loop over the elements.
        for(std::size_t j = 0; j < mesh.elements.size(); ++j) {

            // Local dofs.
            std::size_t element_dofs = mesh.elements[j].dofs();

            // Global matrix indices.
            std::vector<std::size_t> indices;

            for(std::size_t k = 0; k < element_dofs; ++k)
                indices.emplace_back(starts[j] + k);
            
            // Polygon.
            Polygon polygon = mesh.element(j);

            // Element sub-triangulation.
            std::vector<Polygon> triangles = triangulate(polygon);

            // Local coefficients.
            Vector<Real> local_coefficients{element_dofs};

            // Loop over the sub-triangulation.
            for(std::size_t k = 0; k < triangles.size(); ++k) {

                // Triangle.
                Polygon triangle = triangles[k];

                // Jacobian.
                Matrix<Real> jacobian{2, 2};

                jacobian(0, 0) = triangle.points[1][0] - triangle.points[0][0];
                jacobian(0, 1) = triangle.points[2][0] - triangle.points[0][0];
                jacobian(1, 0) = triangle.points[1][1] - triangle.points[0][1];
                jacobian(1, 1) = triangle.points[2][1] - triangle.points[0][1];

                // Jacobian's determinant.
                Real jacobian_det = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);

                // Translation.
                Vector<Real> translation{2};

                translation[0] = triangle.points[0][0];
                translation[1] = triangle.points[0][1];

                // Physical nodes.
                Vector<Real> physical_x{nodes_x_2d.length};
                Vector<Real> physical_y{nodes_y_2d.length};

                for(std::size_t l = 0; l < physical_x.length; ++l) {
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
                Matrix<Real> scaled_phi = phi;

                for(std::size_t l = 0; l < scaled_phi.columns; ++l)
                    scaled_phi.column(l, scaled_phi.column(l) * scaled);

                // function solution.
                Vector<Real> local_function = function(physical_x, physical_y, t);

                // Local coefficients.
                local_coefficients += scaled_phi.transpose() * local_function;
            }

            // Update.
            coefficients(indices, local_coefficients);
        }

        return coefficients;
    }

    /**
     * @brief Evaluate the modal coefficient of the source term.
     * 
     * @param mesh 
     * @param function Source function.
     * @param t Time step.
     * @param Alpha 
     * @param D 
     * @return Vector<Real> 
     */
    Vector<Real> evaluateSource(const Mesh &mesh, const FKPPSource &function, const Real &t, const TriFunctor &Alpha, const TriFunctor &D) {

        // Number of quadrature nodes.
        std::size_t degree = GAUSS_ORDER;

        // Coefficients.
        Vector<Real> coefficients{mesh.dofs()};

        // Quadrature nodes.
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(degree);

        // Starting indices.
        std::vector<std::size_t> starts;
        starts.emplace_back(0);

        for(std::size_t j = 1; j < mesh.elements.size(); ++j)
            starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

        // Loop over the elements.
        for(std::size_t j = 0; j < mesh.elements.size(); ++j) {

            // Local dofs.
            std::size_t element_dofs = mesh.elements[j].dofs();

            // Global matrix indices.
            std::vector<std::size_t> indices;

            for(std::size_t k = 0; k < element_dofs; ++k)
                indices.emplace_back(starts[j] + k);
            
            // Polygon.
            Polygon polygon = mesh.element(j);

            // Element sub-triangulation.
            std::vector<Polygon> triangles = triangulate(polygon);

            // Local coefficients.
            Vector<Real> local_coefficients{element_dofs};

            // Loop over the sub-triangulation.
            for(std::size_t k = 0; k < triangles.size(); ++k) {

                // Triangle.
                Polygon triangle = triangles[k];

                // Jacobian.
                Matrix<Real> jacobian{2, 2};

                jacobian(0, 0) = triangle.points[1][0] - triangle.points[0][0];
                jacobian(0, 1) = triangle.points[2][0] - triangle.points[0][0];
                jacobian(1, 0) = triangle.points[1][1] - triangle.points[0][1];
                jacobian(1, 1) = triangle.points[2][1] - triangle.points[0][1];

                // Jacobian's determinant.
                Real jacobian_det = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);

                // Translation.
                Vector<Real> translation{2};

                translation[0] = triangle.points[0][0];
                translation[1] = triangle.points[0][1];

                // Physical nodes.
                Vector<Real> physical_x{nodes_x_2d.length};
                Vector<Real> physical_y{nodes_y_2d.length};

                for(std::size_t l = 0; l < physical_x.length; ++l) {
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
                Matrix<Real> scaled_phi = phi;

                for(std::size_t l = 0; l < scaled_phi.columns; ++l)
                    scaled_phi.column(l, scaled_phi.column(l) * scaled);

                // function solution.
                Vector<Real> Dext = D(physical_x, physical_y, t);
                Vector<Real> alpha = Alpha(physical_x, physical_y, t);
                Vector<Real> local_function = function(physical_x, physical_y, t, alpha, Dext);

                // Local coefficients.
                local_coefficients += scaled_phi.transpose() * local_function;
            }

            // Update.
            coefficients(indices, local_coefficients);
        }

        return coefficients;

    }

    /**
     * @brief Get the initial condition for the Fisher-KPP equation.
     * 
     * @param mesh MEsh.
     * @param mass Projection matrix.
     * @param ch Function to evaluate.
     * @param t Time step.
     * @return std::array<Vector<Real>, 2> 
     */
    std::array<Vector<Real>, 2> EvaluateICFKPP(const Mesh &mesh, const Sparse<Real> &mass, const TriFunctor &ch, const Real &t) {

        // Mass blocks.
        auto blocks = block_mass(mesh);

        // Initial condition.
        Vector<Real> c_h = evaluateCoeff(mesh, ch, 0.0);
        Vector<Real> c_hh = evaluateCoeff(mesh, ch, 0.0-t);

        // Projection for modal coordinates.
        Vector<Real> c_old = (norm(c_h) > TOLERANCE) ? solve(mass, c_h, blocks, DB) : Vector<Real>{mesh.dofs()};
        Vector<Real> c_oold = (norm(c_hh) > TOLERANCE) ? solve(mass, c_hh, blocks, DB) : Vector<Real>{mesh.dofs()};

        return {c_old, c_oold};

    }

    /**
     * @brief Get the initial condition for the Heat equation.
     * 
     * @param mesh Mesh.
     * @param mass Projection matrix.
     * @param ch Function to evaluate.
     * @return Vector<Real> 
     */
    Vector<Real> EvaluateICHeat(const Mesh &mesh, const Sparse<Real> &mass, const TriFunctor &ch) {

        // Mass blocks.
        auto blocks = block_mass(mesh);

        // Initial condition.
        Vector<Real> c_h = evaluateCoeff(mesh, ch, 0.0);

        // Projection for modal coordinates.
        Vector<Real> c_old = (norm(c_h) > TOLERANCE) ? solve(mass, c_h, blocks, DB) : Vector<Real>{mesh.dofs()};

        return c_old;

    }

}