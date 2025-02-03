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

  // To save typing the full qualified names.
  using namespace pacs;

  // Retrieve problem data from structure.
  DataHeat data;

  std::ostringstream oss;
  oss << "output/square_heat_hp_" << data.degree;
  std::ofstream output(oss.str() + ".error");

  output << "Square domain - hp-adaptive refinement." << "\n";

  std::cout << "Square domain - hp-adaptive refinement." << std::endl;
  std::cout << "Output under " << oss.str() << ".error." << std::endl;

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

  // Parameters.
  int counter = 1;
  Real t = 0.0;
  int steps = static_cast<int>(round(data.t_f / data.dt));

  for (int i = 1; i <= steps; i++) {

    // Time step.
    t += data.dt;
    std::cout << "TIME: " << (heat->t() = t) << std::endl;

    // Update forcing term.
    Vector<Real> F_old = heat->forcing();
    heat->assembly_force(data, mesh);

    // Linear system equation solution.
    Vector<Real> ch = heat->solver(data, mesh, ch_old, F_old);

    if (counter % data.VisualizationStep == 0) {

      // Errors.
      HeatError error(mesh);

      // Compute error.
      error.error(data, mesh, *heat, ch);

      // Output.
      output << "\n"
             << error << "\n"
             << "Residual: "
             << norm(heat->M() * ch + heat->A() * ch - heat->forcing())
             << std::endl;
    }

    // Compute estimator.
    HeatEstimator estimator(mesh);
    estimator.computeEstimates(data, *heat, ch, ch_old);

    // hp-Refine and prolong solution.
    auto [h_mask, p_mask] = estimator.find_elem_to_refine(estimator);

    // p refinement.
    estimator.mesh_refine_degree(p_mask);
    Mesh new_mesh = estimator.mesh();
    
    heat = std::make_unique<Heat>(new_mesh);
    heat->assembly_mass(new_mesh);
    
    ch_old.resize(new_mesh.dofs());
    ch_old = heat->prolong_solution_p(new_mesh, mesh, heat->M(), ch, p_mask);
    
    mesh = std::move(new_mesh);

    // h refinement.
    estimator.mesh_refine_degree(p_mask);
    new_mesh = std::move(estimator.mesh());

    heat = std::make_unique<Heat>(new_mesh);
    heat->assembly(data, new_mesh);
    
    ch_old.resize(new_mesh.dofs());
    ch_old = heat->prolong_solution_p(new_mesh, mesh, heat->M(), ch, p_mask);

    // Update forcing.
    heat->forcing().resize(mesh.dofs());
    heat->assembly_force(data, new_mesh);

    // Update mesh.
    mesh = std::move(new_mesh);

    ++counter;
  }

// Solution structure (output).
#ifndef NSOLUTIONS
  HeatSolution solution{mesh};
  solution.computeSolution(data, mesh, heat, ch_old);
  std::string solfile =
      "output/square_s_" + std::to_string(data.degree) + ".sol";
  solution.write(solfile);
#endif
}