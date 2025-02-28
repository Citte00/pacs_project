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

#include <chrono>
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
  oss << "output/laplacian_" << data.degree;
  std::ofstream output(oss.str() + ".error");

  output << "Square domain - uniform refinement." << "\n";

  std::cout << "Square domain - uniform refinement." << std::endl;
  std::cout << "Output under output/square_s_" + std::to_string(data.degree) +
                   ".error."
            << std::endl;

  // Domain.
  Polygon domain{data.domain};

  // Diagrams.
  std::vector<std::vector<Polygon>> diagrams;

  diagrams.emplace_back(mesh_diagram("meshes/square/square_125.poly"));
  diagrams.emplace_back(mesh_diagram("meshes/square/square_250.poly"));
  diagrams.emplace_back(mesh_diagram("meshes/square/square_500.poly"));
  diagrams.emplace_back(mesh_diagram("meshes/square/square_1000.poly"));
  diagrams.emplace_back(mesh_diagram("meshes/square/square_2000.poly"));
  diagrams.emplace_back(mesh_diagram("meshes/square/square_4000.poly"));

  // Test.
  for (std::size_t j = 0; j < diagrams.size(); ++j) {

    // Mesh.
    Mesh mesh{domain, std::move(diagrams[j]), data.degree};

    // Matrices.
    Laplace<Real> laplacian(mesh);
    laplacian.assembly(data, mesh);

    // Forcing term.
    Vector<Real> forcing = laplacian.assembly_force(data, mesh);

    // Linear system solution.
    Vector<Real> numerical = laplacian.solver(mesh, forcing);

    // Errors.
    LaplaceError<Real> error(mesh);

    // Solution structure (output).
    LaplaceSolution<Real> solution{mesh};
    solution.solution(data, mesh, numerical);
    std::string solfile = "output/square_s_" + std::to_string(data.degree) +
                          "_" + std::to_string(j) + ".sol";
    solution.write(solfile);

    // Compute error.
    error.error(data, mesh, laplacian, numerical);

    // Output.
    output << "\n" << error << "\n";

    output << "Laplacian: " << laplacian.A().m_rows << " x "
           << laplacian.A().m_columns << "\n";
    output << "Residual: " << norm(laplacian.A() * numerical - forcing)
           << std::endl;
  }
}