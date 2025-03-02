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

using namespace pacs;

int main() {

    // Constructs a domain.
    Point a{0.0, 0.0};
    Point b{1.0, 0.0};
    Point c{1.0, 1.0};
    Point d{0.0, 1.0};

    Polygon domain{{a, b, c, d}};
    
    // Constructing a mesh.
    Mesh mesh{domain, mesh_diagram("meshes/square/square_300.poly")};
    std::cout << "Dofs: " << mesh.dofs() << std::endl;

    // Refining some elements.
    pacs::Mask refinement(300);

    refinement[0] = true;
    refinement[1] = true;
    refinement[2] = true;

    // Refinement.
    LaplaceEstimator<Real> estimator(mesh);
    estimator.mesh_refine_size(refinement);

    std::cout << "Dofs: " << mesh.dofs() << std::endl;

    // Mesh output.
    mesh.write("output/refined.poly");
}