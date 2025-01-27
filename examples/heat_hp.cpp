/**
 * @file heat_h.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Heat equation hp-adaptive refinement.
 * @date 2025-01-26
 *
 * @copyright Copyright (c) 2025
 *
 */
#include <PacsHPDG.hpp>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char **argv) {

  // Retrieve problem data from structure.
  pacs::DataHeat data;

  std::ofstream output{"output/square_heat_hp_" + std::to_string(data.degree) +
                       ".error"};

  std::string outputVTK =
      "output/square_heat_hp_" + std::to_string(data.degree) + ".poly";

  output << "Square domain - hp-adaptive refinement." << "\n";

  std::cout << "Square domain - hp-adaptive refinement." << std::endl;
  std::cout << "Output under output/square_heat_hp_" +
                   std::to_string(data.degree) + ".error."
            << std::endl;

  // Domain.
  pacs::Polygon domain{data.domain};

  // Diagrams.
  std::vector<pacs::Polygon> diagram =
      pacs::mesh_diagram("meshes/square/square_125.poly");

  // Mesh.
  pacs::Mesh mesh{domain, diagram, data.degree};

  // Matrices.
  pacs::Heat heat(mesh);
  heat.assembly(data, mesh);

  // Initial condition.
  pacs::Vector<pacs::Real> ch_old = heat.modal(mesh, data.c_ex);

  // Forcing term.
  heat.assembly_force(data, mesh);

  int counter = 1;

  int steps = static_cast<int>(round(data.t_f / data.dt));
  for (int i = 1; i <= steps; i++) {

    // Time step.
    heat.t() += data.dt;
    std::cout << "TIME: " << heat.t() << std::endl;

    // Update forcing term.
    pacs::Vector<pacs::Real> F_old = heat.forcing();
    heat.assembly_force(data, mesh);

    // Linear system equation solution.
    heat.solver(data, mesh, ch_old, F_old);

    // Errors.
    pacs::HeatError error(mesh);

    // Compute error.
    error.computeErrors(data, mesh, heat);

    if (counter % data.VisualizationStep == 0) {
      // Output.
      output << "\n" << error << "\n";

      output << "Residual: "
             << pacs::norm(heat.M() * heat.ch() + heat.A() * heat.ch() -
                           heat.forcing())
             << std::endl;
    }

    // Compute estimates.
    pacs::HeatEstimator estimates(mesh);
    estimates.computeEstimates(data, heat, ch_old);
    estimates.write(outputVTK, true);

    // hp-Refine and prolong solution.
    estimates.mesh_refine(heat, mesh);

    // Update mesh.
    mesh = estimates.mesh();

    // Update matrices.
    heat.assembly(data, mesh);
    auto blocks = heat.block_mass(mesh);

    // Update solution.
    ch_old.resize(heat.ch().length);
    ch_old = pacs::solve(heat.M(), heat.ch(), blocks, pacs::DB);

    // Update forcing.
    heat.forcing().resize(mesh.dofs());
    heat.assembly_force(data, mesh);

    ++counter;
  }

// Solution structure (output).
#ifndef NSOLUTIONS
  pacs::HeatSolution solution{mesh};
  solution.computeSolution(data, mesh, heat);
  std::string solfile =
      "output/square_s_" + std::to_string(data.degree) + ".sol";
  solution.write(solfile);
#endif
}