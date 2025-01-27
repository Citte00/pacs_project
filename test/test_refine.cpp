/**
 * @file test_refine.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <iostream>
#include <vector>

int main() {

    // Constructs a domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};
    
    // Constructing a mesh.
    pacs::Mesh mesh{domain, pacs::mesh_diagram("meshes/square/square_300.poly")};
    std::cout << "Dofs: " << mesh.dofs() << std::endl;

    // Refining some elements.
    pacs::Mask refinement(300);

    refinement[0] = true;
    refinement[1] = true;
    refinement[2] = true;

    // Refinement.
    pacs::LaplaceEstimator estimator(mesh);
    estimator.mesh_refine_size(mesh, refinement);

    std::cout << "Dofs: " << mesh.dofs() << std::endl;

    // Mesh output.
    mesh.write("output/refined.poly");

}