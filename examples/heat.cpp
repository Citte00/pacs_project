/**
 * @file heat.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Heat equation, uniform refinement.
 * @date 2025-01-20
 * 
 * @copyright Copyright (c) 2025
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
  DataHeat data;

  std::ostringstream oss;
  oss << "output/heat_p_" << data.degree;
  std::ofstream output(oss.str() + ".error");

  output << "Square domain - uniform refinement." << "\n";

  std::cout << "Square domain - uniform refinement." << std::endl;
  std::cout << "Output under " << oss.str() << ".error." << std::endl;

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
    std::cout << "Elements: " << mesh.elements.size() << std::endl;

    // Matrices.
    Heat<Real> heat(mesh);
    heat.assembly(data, mesh);

    // Initial condition.
    Vector<Real> ch_old = heat.modal(mesh, data.c_ex);

    // Forcing term.
    heat.assembly_force(data, mesh);
    Vector<Real> F_old{heat.forcing().length};

    int steps = static_cast<int>(round((data.t_f - data.t_0) / data.dt));
    for (int i = 1; i <= steps; ++i) {

      // Time step.
      heat.t() += data.dt;
      std::cout << "TIME: " << heat.t() << std::endl;

      // Update forcing term.
      F_old = heat.forcing();
      heat.assembly_force(data, mesh);

      // Linear system equation solution.
      Vector<Real> ch = heat.solver(data, mesh, ch_old, F_old);

      // Update solution.
      ch_old = ch;
    }

    // Errors.
    HeatError<Real> error(mesh);

// Solution structure (output).
#ifndef NSOLUTIONS
    HeatSolution<Real> solution{mesh};
    solution.solution(data, mesh, heat, ch_old);
    std::string solfile = "output/square_s_" + std::to_string(data.degree) +
                          "_" + std::to_string(j) + ".sol";
    solution.write(solfile);
#endif

    // Compute error.
    error.error(data, mesh, heat, ch_old);

    // Output.
    output << "\n"
           << error << "\n"
           << "Residual: "
           << norm(heat.M() * ch_old + heat.A() * ch_old - heat.forcing())
           << std::endl;
  }
}