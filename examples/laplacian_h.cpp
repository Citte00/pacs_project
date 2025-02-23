/**
 * @file square_h.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Element size adaptive refinement.
 * @date 2024-05-31
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <PacsHPDG.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// To save typing the full qualified names.
using namespace pacs;

int main(int argc, char **argv) {

  // Retrieve problem data from structure.
  DataLaplace data;

  std::ostringstream oss;
  oss << "output/laplacian_h_" << data.elements << "@" << data.degree;
  std::ofstream output(oss.str() + ".error");

  output << "Square domain - element size adaptive refinement." << "\n";

  std::cout << "Square domain - element size adaptive refinement." << std::endl;
  std::cout << "Output under " << oss.str() << ".error." << std::endl;

  // Domain.
  Polygon domain{data.domain};

  // Diagram.
  std::vector<Polygon> diagram = mesh_diagram(
      "meshes/square/square_" + std::to_string(data.elements) + ".poly");

  // Refinement percentage.
  Real refine = 0.75L;

  // Mesh.
  Mesh mesh{domain, std::move(diagram), data.degree};

  // Sequence of meshes.
  for (std::size_t index = 0; index < TESTS_MAX; ++index) {

    // Verbosity.
    std::cout << "\nDEGREE: " << data.degree << "\nINDEX: " << index << "\n"
              << std::endl;

    // Mesh output.
    std::string polyfile = "output/square_h_" + std::to_string(data.elements) +
                           "@" + std::to_string(data.degree) + "_" +
                           std::to_string(index) + ".poly";
    mesh.write(polyfile);

    // Matrices.
    Laplace<Real> laplacian(mesh);
    laplacian.assembly(data, mesh);

    // Forcing term.
    Vector<Real> forcing = laplacian.assembly_force(data, mesh);

    // Linear system solution.
    Vector<Real> numerical = laplacian.solver(mesh, forcing);

// Solution structure (output).
#ifndef NSOLUTIONS
    LaplaceSolution<Real> solution{mesh};
    solution.solution(data, mesh, numerical);
    std::string solfile = "output/square_h_" + std::to_string(data.elements) +
                          "@" + std::to_string(data.degree) + "_" +
                          std::to_string(index) + ".sol";
    solution.write(solfile);
#endif

    // Errors.
    LaplaceError<Real> error(mesh);
    error.error(data, mesh, laplacian, numerical);

    // Output.
    output << "\n"
           << error << "\n"
           << "Laplacian: " << laplacian.A().rows << " x "
           << laplacian.A().columns << "\n"
           << "Residual: " << norm(laplacian.A() * numerical - forcing)
           << std::endl;

    // Exit.
    if (error.dofs() > DOFS_MAX)
      break;

    // Refinement.
    LaplaceEstimator<Real> estimator(mesh);
    estimator.mesh_refine_size(error.L2errors() >
                               refine * max(error.L2errors()));

    // Update mesh.
    mesh = estimator.mesh();
  }
}
