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

// To save typing the full qualified names.
using namespace pacs;

int main(int argc, char **argv) {

  // Retrieve problem data from structure.
  DataLaplace data;

  // "Splash".
  std::ostringstream oss;
  oss << "output/laplacian_hp_" << data.elements << "@" << data.degree;
  std::ofstream output(oss.str() + ".error");

  output << "Square domain - hp-adaptive refinement." << "\n";

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
  Mesh mesh{domain, std::move(diagram), data.degree};

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
    std::string solfile = "output/square_hp_" + std::to_string(data.elements) +
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
           << "Laplacian: " << laplacian.A().m_rows << " x "
           << laplacian.A().m_columns << "\n"
           << "Residual: " << norm(laplacian.A() * numerical - forcing)
           << std::endl;

    // Exit.
    if (error.dofs() > DOFS_MAX)
      break;

    // Estimator.
    LaplaceEstimator<Real> estimator(mesh);
    estimator.computeEstimates(data, laplacian, numerical);

    // Refinement.
    auto [h_mask, p_mask] = estimator.find_elem_to_refine();
    estimator.mesh_refine_degree(p_mask);
    estimator.mesh_refine_size(h_mask);

    // Update mesh.
    mesh = estimator.mesh();
  }
}
