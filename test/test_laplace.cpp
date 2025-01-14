#include <PacsHPDG.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char **argv) {

  // Retrieve data.
  pacs::DataLap data;

  std::ofstream output{"output/square_s_" + std::to_string(data.N) + "@" +
                       std::to_string(data.degree) + ".error"};

  output << "Square domain - uniform refinement." << "\n";

  std::cout << "Square domain - uniform refinement." << std::endl;
  std::cout << "Output under output/square_" + std::to_string(data.degree) +
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

    std::cout << "INDEX: " << j << std::endl;

    // Mesh.
    pacs::Mesh mesh{domain, diagrams[j], data.degree};

    // Matrices.
    pacs::Laplace laplacian(mesh);
    laplacian.assembly(data, mesh);

    // Forcing term.
    laplacian.assemblyforce(data, mesh);

    // Linear system solution.
    laplacian.solver(mesh);

    // Errors.
    pacs::LaplaceError error(mesh);

    error.computeErrors(data, mesh, laplacian);

    // Output.
    error.print(output, mesh);
  }
}