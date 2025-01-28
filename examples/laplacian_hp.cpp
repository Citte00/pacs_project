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

  // To save typing the full qualified names.
  using namespace pacs;

  // Retrieve problem data from structure.
  DataLaplace data;

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
  Polygon domain{data.domain};

  // Diagram.
  std::vector<Polygon> diagram = mesh_diagram(
      "meshes/square/square_" + std::to_string(data.elements) + ".poly");

  // Mesh.
  Mesh mesh{domain, diagram, data.degree};

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
    Laplace laplacian(mesh);
    laplacian.assembly(data, mesh);

    // Forcing term.
    Vector<Real> forcing = laplacian.assembly_force(data, mesh);

    // Linear system solution.
    Vector<Real> numerical = laplacian.solver(mesh, forcing);

// Solution structure (output).
#ifndef NSOLUTIONS
    LaplaceSolution solution{mesh};
    solution.computeSolution(data, mesh, numerical);
    std::string solfile = "output/square_hp_" + std::to_string(data.elements) +
                          "@" + std::to_string(data.degree) + "_" +
                          std::to_string(index) + ".sol";
    solution.write(solfile);
#endif

    // Errors.
    LaplaceError error(mesh);
    error.computeErrors(data, mesh, laplacian, numerical);

    // Output.
    output << "\n" << error << "\n";

    output << "Laplacian: " << laplacian.A().rows << " x "
           << laplacian.A().columns << "\n";
    output << "Residual: " << norm(laplacian.A() * numerical - forcing)
           << std::endl;

    // Exit.
    if (error.dofs() > DOFS_MAX)
      break;

    // Estimator.
    LaplaceEstimator estimator(mesh);
    estimator.computeEstimates(data, laplacian, numerical);

    // Refinement.
    estimator.mesh_refine(estimator);

    // Update mesh.
    mesh = estimator.mesh();
  }
}
