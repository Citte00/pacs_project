/**
 * @file square_smooth.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Uniform refinement.
 * @date 2024-07-03
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <PacsHPDG.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono>

int main(int argc, char **argv) {

  // Retrieve problem data from structure.
  pacs::DataLaplace data;

  std::ofstream output{"output/square_s_" + std::to_string(data.degree) +
                       ".error"};

  output << "Square domain - uniform refinement." << "\n";

  std::cout << "Square domain - uniform refinement." << std::endl;
  std::cout << "Output under output/square_s_" + std::to_string(data.degree) +
                   ".error."
            << std::endl;

  // Domain.
  pacs::Polygon domain{data.domain};

  // Diagrams.
  std::vector<std::vector<pacs::Polygon>> diagrams;

  diagrams.emplace_back(pacs::mesh_diagram("meshes/square/square_125.poly"));
  diagrams.emplace_back(pacs::mesh_diagram("meshes/square/square_250.poly"));
  diagrams.emplace_back(pacs::mesh_diagram("meshes/square/square_500.poly"));
  diagrams.emplace_back(pacs::mesh_diagram("meshes/square/square_1000.poly"));
  diagrams.emplace_back(pacs::mesh_diagram("meshes/square/square_2000.poly"));
  diagrams.emplace_back(pacs::mesh_diagram("meshes/square/square_4000.poly"));

  // Test.
  for (std::size_t j = 0; j < diagrams.size(); ++j) {

    // Mesh.
    pacs::Mesh mesh{domain, diagrams[j], data.degree};

    // Matrices.
    pacs::Laplace laplacian(mesh);
    laplacian.assembly(data, mesh);

    // Forcing term.
    pacs::Vector<pacs::Real> forcing = laplacian.assembly_force(data, mesh);

    // Linear system solution.
    pacs::Vector<pacs::Real> numerical = laplacian.solver(mesh, forcing);

    // Errors.
    pacs::LaplaceError error(mesh);

// Solution structure (output).
#ifndef NSOLUTIONS
    pacs::LaplaceSolution solution{mesh};
    solution.computeSolution(data, mesh, numerical);
    std::string solfile = "output/square_s_" + std::to_string(data.degree) + "_" +
                          std::to_string(j) + ".sol";
    solution.write(solfile);
#endif

    // Compute error.
    error.computeErrors(data, mesh, laplacian, numerical);

    // Output.
    output << "\n" << error << "\n";

    output << "Laplacian: " << laplacian.A().rows << " x "
           << laplacian.A().columns << "\n";
    output << "Residual: " << pacs::norm(laplacian.A() * numerical - forcing)
           << std::endl;
  }
}