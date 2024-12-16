/**
 * @file Fisher.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-11-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include <PacsHPDG.hpp>

namespace pacs{

    /**
     * @brief Assembly of matrices for Fisher-KPP equation
     * 
     * @param mesh 
     * @param Alpha alpha
     * @param D D_ext
     * @param penalty_coefficient 
     * @return std::array<Sparse<Real>, 4> 
     */
    std::array<Sparse<Real>, 4> fisher(const Mesh &mesh, const size_t &degree, const TriFunctor &Alpha, const TriFunctor &D, const Real &penalty_coefficient)
    {
        #ifndef NVERBOSE
        std::cout << "Computing the Fisher-KPP matrices." << std::endl;
        #endif

        // Number of quadrature nodes.
        std::size_t nqn = 2*degree + 1;

        // Quadrature nodes.
        auto [nodes_1d, weights_1d] = quadrature_1d(nqn);
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn);

        // Degrees of freedom.
        std::size_t dofs = mesh.dofs();

        // Neighbours.
        std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

        // Matrices.
        Sparse<Real> M_prj{dofs, dofs};
        Sparse<Real> M{dofs, dofs};
        Sparse<Real> A{dofs, dofs};
        Sparse<Real> IA{dofs, dofs};
        Sparse<Real> SA{dofs, dofs};

        // Starting indices.
        std::vector<std::size_t> starts;
        starts.emplace_back(0);

        for(std::size_t j = 1; j < mesh.elements.size(); ++j)
            starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

        // Volume integrals.

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

            // Local matrices.
            Matrix<Real> local_M_prj{element_dofs, element_dofs};
            Matrix<Real> local_M{element_dofs, element_dofs};
            Matrix<Real> local_A{element_dofs, element_dofs};

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

                // Param initialization.
                Vector<Real> Dext = D(physical_x, physical_y, 0.0);
                Vector<Real> alpha = Alpha(physical_x, physical_y, 0.0);

                // Basis functions.
                auto [phi, gradx_phi, grady_phi] = basis_2d(mesh, j, {physical_x, physical_y});

                // Some products.
                Matrix<Real> scaled_phi{phi};
                Matrix<Real> scaled_phi_alpha{phi};
                Matrix<Real> scaled_gradx{gradx_phi};
                Matrix<Real> scaled_grady{grady_phi};

                for(std::size_t l = 0; l < scaled_gradx.columns; ++l) {
                    scaled_phi.column(l, scaled_phi.column(l) * scaled);
                    scaled_phi_alpha.column(l, (alpha * scaled_phi_alpha.column(l)) * scaled);
                    scaled_gradx.column(l, (Dext * scaled_gradx.column(l)) * scaled);
                    scaled_grady.column(l, (Dext * scaled_grady.column(l)) * scaled);
                }

                // Local matrix assembly.
                local_M_prj += scaled_phi.transpose() * phi;
                local_M += scaled_phi_alpha.transpose() * phi;
                local_A += scaled_gradx.transpose() * gradx_phi + scaled_grady.transpose() * grady_phi;
            }

            // Global matrix assembly.
            M_prj.insert(indices, indices, local_M_prj);
            M.insert(indices, indices, local_M);
            A.insert(indices, indices, local_A);

            // Face integrals.

            // Local matrices.
            Matrix<Real> local_IA{element_dofs, element_dofs};
            Matrix<Real> local_SA{element_dofs, element_dofs};

            // Element's neighbours.
            std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

            // Local matrices for neighbours.
            std::vector<Matrix<Real>> local_IAN;
            std::vector<Matrix<Real>> local_SAN;

            // Penalties.
            Vector<Real> penalties = penalty(mesh, j, penalty_coefficient);

            // Edges.
            std::vector<Segment> edges{polygon.edges()};

            // Loop over faces.
            for(std::size_t k = 0; k < element_neighbours.size(); ++k) {

                // Neighbour information.
                auto [edge, neighbour, n_edge] = element_neighbours[k];

                // Edge geometry.
                Segment segment{edges[k]};

                // Edge's normal.
                Vector<Real> normal_vector{2};

                normal_vector[0] = segment[1][1] - segment[0][1];
                normal_vector[1] = segment[0][0] - segment[1][0];

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

                for(std::size_t l = 0; l < nodes_1d.length; ++l) {
                    Vector<Real> node{2};

                    node[0] = nodes_1d[l];

                    Vector<Real> transformed = jacobian * node + translation;

                    physical_x[l] = transformed[0];
                    physical_y[l] = transformed[1];
                }

                // Weights scaling.
                Vector<Real> scaled = std::abs(segment) * weights_1d;

                // Param initialization.
                Vector<Real> Dext = D(physical_x, physical_y, 0.0);

                // Basis functions.
                auto [phi, gradx_phi, grady_phi] = basis_2d(mesh, j, {physical_x, physical_y});

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

                if(neighbour == -1) { // Boundary edge.

                    local_IA += scaled_grad.transpose() * phi;
                    local_SA += (penalties[k] * scaled_phi).transpose() * phi;

                    // Empty small matrices.
                    local_IAN.emplace_back(Matrix<Real>{1, 1});
                    local_SAN.emplace_back(Matrix<Real>{1, 1});

                } else {

                    local_IA += 0.5 * scaled_grad.transpose() * phi;
                    local_SA += (penalties[k] * scaled_phi).transpose() * phi;

                    // Neighbour's basis function.
                    Matrix<Real> n_phi = basis_2d(mesh, neighbour, {physical_x, physical_y})[0];

                    // Neighbour's local matrix.
                    local_IAN.emplace_back(- 0.5 * scaled_grad.transpose() * n_phi);
                    local_SAN.emplace_back(- (penalties[k] * scaled_phi).transpose() * n_phi);
                }
            }

            IA.insert(indices, indices, local_IA);
            SA.insert(indices, indices, local_SA);

            // Neighbouring DG matrices assembly.
            for(std::size_t k = 0; k < element_neighbours.size(); ++k) {
                if(element_neighbours[k][1] == -1)
                    continue;

                std::vector<std::size_t> n_indices;
                std::size_t n_index = element_neighbours[k][1];
                std::size_t n_dofs = mesh.elements[n_index].dofs(); // Neighbour's dofs.

                for(std::size_t h = 0; h < n_dofs; ++h)
                    n_indices.emplace_back(starts[n_index] + h);

                IA.add(indices, n_indices, local_IAN[k]);
                SA.add(indices, n_indices, local_SAN[k]);
            }
        }

        // Matrices.
        Sparse<Real> mass_prj = M_prj;
        Sparse<Real> mass = M;
        Sparse<Real> dg_stiffness = A + SA;
        Sparse<Real> stiffness = dg_stiffness - IA - IA.transpose();

        // Compression.
        mass_prj.compress();
        mass.compress();
        dg_stiffness.compress();
        stiffness.compress();
        
        return {mass_prj, mass, stiffness, dg_stiffness};
    }

    /**
     * @brief Assembly of the non-linear matrix for Fisher-KPP equation
     * 
     * @param mesh 
     * @param Alpha 
     * @param uh 
     * @param penalty_coefficient 
     * @return Sparse<Real> 
     */
    Sparse<Real> NLfisher(const Mesh &mesh, const size_t &degree, const TriFunctor &Alpha, const Vector<Real> &uh, const Real &penalty_coefficient)
    {
        #ifndef NVERBOSE
        std::cout << "Computing the Fisher-KPP matrices." << std::endl;
        #endif

        // Number of quadrature nodes.
        std::size_t nqn = 2*degree + 1;

        // Quadrature nodes.
        auto [nodes_1d, weights_1d] = quadrature_1d(nqn);
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(nqn);

        // Degrees of freedom.
        std::size_t dofs = mesh.dofs();

        // Neighbours.
        std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

        // Matrices.
        Sparse<Real> M_star{dofs, dofs};

        // Starting indices.
        std::vector<std::size_t> starts;
        starts.emplace_back(0);

        for(std::size_t j = 1; j < mesh.elements.size(); ++j)
            starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

        // Volume integrals.

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

            // Local matrices.
            Matrix<Real> local_M{element_dofs, element_dofs};

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

                // Param initialization.
                Vector<Real> alpha = Alpha(physical_x, physical_y, 0.0);
                Vector<Real> c_star = uh(indices);

                // Basis functions.
                auto [phi, gradx_phi, grady_phi] = basis_2d(mesh, j, {physical_x, physical_y});

                // Some products.
                Matrix<Real> scaled_phi{phi};
                Vector<Real> c_star_loc = phi * c_star;

                for(std::size_t l = 0; l < scaled_phi.columns; ++l) {
                    scaled_phi.column(l, (alpha * scaled_phi.column(l)) * c_star_loc * scaled);
                }

                // Local matrix assembly.
                local_M += scaled_phi.transpose() * phi;

            }

            // Global matrix assembly.
            M_star.insert(indices, indices, local_M);
            
        }

        // Matrices.
        Sparse<Real> non_lin = M_star;
    
        // Compression.
        non_lin.compress();
        
        return non_lin;
    }

}