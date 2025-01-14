/**
 * @file square_hp.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. hp-Adaptive refinement with estimator.
 * @date 2024-06-27
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <PacsHPDG.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char **argv) {

  // Retrieve problem data from structure.
  pacs::DataLaplace data;

  // "Splash".
  std::ofstream output{"output/square_hp_" + std::to_string(data.elements) +
                       "@" + std::to_string(data.degree) + ".error"};
  output << "Square domain - hp-adaptive refinement with estimator." << "\n";

  std::cout << "Square domain - hp-adaptive refinement with estimator."
            << std::endl;
  std::cout << "Output under output/square_hp_" + std::to_string(data.degree) +
                   ".error."
            << std::endl;

  // Domain.
  pacs::Polygon domain{data.domain};

  // Diagram.
  std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(
      "meshes/square/square_" + std::to_string(data.elements) + ".poly");

  // Mesh.
  pacs::Mesh mesh{domain, diagram, data.degree};

  // Sequence of meshes.
  for (std::size_t index = 0; index < TESTS_MAX; ++index) {

    // Verbosity.
    std::cout << "\nSTARTING DEGREE: " << data.degree << "\nINDEX: " << index
              << "\n"
              << std::endl;

    // Mesh output.
    std::string polyfile = "output/square_hp_" + std::to_string(data.elements) +
                           "@" + std::to_string(data.degree) + "_" +
                           std::to_string(index) + ".poly";
    mesh.write(polyfile, true);

    // Matrices.
    pacs::Laplace laplacian(mesh);
    laplacian.assembly(data, mesh);

    // Forcing term.
    pacs::Vector<pacs::Real> forcing = laplacian.forcing(data, mesh);

    // Linear system solution.
    pacs::Vector<pacs::Real> numerical = laplacian.lapsolver(mesh, forcing);

// Solution structure (output).
#ifndef NSOLUTIONS
    pacs::Solution solution{mesh, numerical, exact};
    std::string solfile = "output/square_hp_" + std::to_string(elements) + "@" +
                          std::to_string(degree) + "_" + std::to_string(index) +
                          ".sol";
    solution.write(solfile);
#endif

    // Errors.
    pacs::LaplaceError error(mesh);

    // Output.
    output << "\n" << error << "\n";

    output << "Laplacian: " << laplacian.A().rows << " x "
           << laplacian.A().columns << "\n";
    output << "Residual: " << pacs::norm(laplacian.A() * numerical - forcing)
           << std::endl;

    // Exit.
    if (error.dofs() > DOFS_MAX)
      break;

    // Estimator.
    pacs::LaplaceEstimator estimator(mesh);
    estimator.computeEstimates(data, mesh, laplacian, numerical);

    // Refinement.
    estimator.mesh_refine(mesh, estimator);
  }
}
