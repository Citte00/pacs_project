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
#include <memory>
#include <iostream>
#include <string>

int main(int argc, char **argv) {

  // To save typing the full qualified names.
  using namespace pacs;

  // Retrieve problem data from structure.
  DataHeat data;

  std::ostringstream oss;
  oss << "output/square_heat_h_" << data.degree;
  std::ofstream output(oss.str() + ".error");

  output << "Square domain - element size adaptive refinement." << "\n";

  std::cout << "Square domain - element size adaptive refinement." << std::endl;
  std::cout << "Output under " << oss.str() << ".error." << std::endl;

  // Domain.
  Polygon domain{data.domain};

  // Diagrams.
  std::vector<Polygon> diagram = mesh_diagram(
      "meshes/square/square_" + std::to_string(data.elements) + ".poly");

  // Mesh.
  Mesh mesh{domain, std::move(diagram), data.degree};

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
  const int dofsLimit = 2 * DOFS_MAX;
  int steps = static_cast<int>(round(data.t_f / data.dt));
  
  for (int i = 1; i <= steps; ++i) {

    // Time step.
    t += data.dt;
    std::cout << "TIME: " << (heat->t() = t) << std::endl;
    // Number of elements.
    std::cout << "Elements: " << mesh.elements.size() << std::endl;

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

    // Compute estimates.
    HeatEstimator estimator(mesh);
    estimator.computeEstimates(data, *heat, ch, ch_old);
    const auto &estimates = estimator.estimates();
    std::string estimate =
        "output/heat_" + std::to_string(mesh.elements.size()) + "@" + std::to_string(data.degree) + ".poly";
    estimator.write(estimate, true);

    // Refine.
    Mask h_mask = estimates > 0.75L * sum(estimates) /
                                              estimator.mesh().elements.size();

    if (std::none_of(h_mask.begin(), h_mask.end(), [](bool v) { return v; }) ||
        estimator.mesh().dofs() >= dofsLimit) {
      ch_old = ch;
      ++counter;
      continue;
    }

    // Refine mesh.
    estimator.mesh_refine_size(h_mask);
    const auto &new_mesh = estimator.mesh();

    // Update matrices.
    heat = std::make_unique<Heat>(new_mesh);
    heat->t() = t;
    heat->assembly(data, new_mesh);

    // Prolong solution.
    ch_old.resize(new_mesh.dofs());
    ch_old = heat->prolong_solution_h(new_mesh, mesh, ch, h_mask);

    // Update forcing.
    heat->assembly_force(data, new_mesh);

    // Update mesh.
    mesh = std::move(new_mesh);

    ++counter;
  }

// Solution structure (output).
#ifndef NSOLUTIONS

#endif

  // Final mesh.
  std::string meshfile =
      "output/square_h_" + std::to_string(mesh.elements.size()) + ".poly";
  mesh.write(meshfile);
}