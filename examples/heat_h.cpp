/**
 * @file heat_h.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Heat equation element size adaptive refinement.
 * @date 2025-01-22
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

  // To save typing the full qualified names.
  using namespace pacs;

  // Retrieve problem data from structure.
  DataHeat data;

  std::ofstream output{"output/square_heat_h_" + std::to_string(data.degree) +
                       ".error"};

  std::string outputVTK =
      "output/square_heat_h_" + std::to_string(data.degree) + ".poly";

  output << "Square domain - element size adaptive refinement." << "\n";

  std::cout << "Square domain - element size adaptive refinement." << std::endl;
  std::cout << "Output under output/square_heat_h_" +
                   std::to_string(data.degree) + ".error."
            << std::endl;

  // Domain.
  Polygon domain{data.domain};

  // Diagrams.
  std::vector<Polygon> diagram = mesh_diagram(
      "meshes/square/square_" + std::to_string(data.elements) + ".poly");

  // Mesh.
  Mesh mesh{domain, diagram, data.degree};

  // Matrices.
  std::unique_ptr<Heat> heat = std::make_unique<Heat>(mesh);
  heat->assembly(data, mesh);

  // Initial condition.
  Vector<Real> ch_old = heat->modal(mesh, data.c_ex);

  // Forcing term.
  heat->assembly_force(data, mesh);

  int counter = 1;
  Real t = 0.0;

  int steps = static_cast<int>(round(data.t_f / data.dt));
  for (int i = 1; i <= steps; i++) {

    // Time step.
    t += data.dt;
    heat->t() = t;
    std::cout << "TIME: " << heat->t() << std::endl;

    std::cout << "Elements: " << mesh.elements.size() << std::endl;

    // Update forcing term.
    Vector<Real> F_old = heat->forcing();
    heat->assembly_force(data, mesh);

    // Linear system equation solution.
    Vector<Real> ch = heat->solver(data, mesh, ch_old, F_old);

    // Errors.
    HeatError error(mesh);

    // Compute error.
    error.computeErrors(data, mesh, *heat, ch);

    if (counter % data.VisualizationStep == 0) {
      // Output.
      output << "\n" << error << "\n";

      output << "Residual: "
             << norm(heat->M() * ch + heat->A() * ch -
                     heat->forcing())
             << std::endl;
    }

    /*
    HeatSolution solution{mesh};
    solution.computeSolution(data, mesh, *heat, ch);
    std::string solfile =
        "output/square_h_" + std::to_string(mesh.elements.size()) + ".sol";
    solution.write(solfile);
    */

    // Compute estimates.
    HeatEstimator estimates(mesh);
    estimates.computeEstimates(data, *heat, ch, ch_old);

    // Refine.
    Mask h_mask = error.L2errors() > 0.75L * max(error.L2errors());
    estimates.mesh_refine_size(h_mask);

    if (estimates.mesh().dofs() >= 24000)
      continue;

    // Update matrices.
    heat.reset(new Heat(estimates.mesh()));
    heat->t() = t;
    heat->assembly(data, estimates.mesh());

    // Prolong solution.
    ch_old.resize(estimates.mesh().dofs());
    ch_old = heat->prolong_solution_h(estimates.mesh(), mesh, heat->M(), ch, h_mask);

    // Update forcing.
    heat->assembly_force(data, estimates.mesh());

    // Update mesh.
    mesh = std::move(estimates.mesh());
    
    /*
    HeatSolution solution2{mesh};
    solution2.computeSolution(data, mesh, *heat, ch_old);
    solfile =
        "output/square_h_" + std::to_string(mesh.elements.size()) + ".sol";
    solution2.write(solfile);
    */

    ++counter;
  }

// Solution structure (output).
#ifndef NSOLUTIONS
  HeatSolution solution{mesh};
  solution.computeSolution(data, mesh, heat);
  std::string solfile =
      "output/square_h_" + std::to_string(mesh.elements.size()) + ".sol";
  solution.write(solfile);
#endif
}