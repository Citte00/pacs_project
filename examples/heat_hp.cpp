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
#include <memory>
#include <string>

// To save typing the full qualified names.
using namespace pacs;

int main(int argc, char **argv) {

  // Retrieve problem data from structure.
  DataHeat data;

  std::ostringstream oss;
  oss << "output/heat_hp_" << data.elements << "@" << data.degree;
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
  Mesh mesh{domain, std::move(diagram), data.degree};

  // Matrices.
  Heat<Real> heat{mesh};
  heat.assembly(data, mesh);

  // Initial condition.
  Vector<Real> ch_old = heat.modal(mesh, data.c_ex);

  // Forcing term.
  heat.assembly_force(data, mesh);

  // Parameters.
  int counter = 1;
  std::size_t dofsLimit = DOFS_MAX;
  std::size_t degree = 0;
  int steps = static_cast<int>(round(data.t_f / data.dt));

  for (int i = 1; i <= steps; i++) {

    // Time step.
    heat.t() += data.dt;
    std::cout << "TIME: " << heat.t() << std::endl;
    // Number of elements.
    std::cout << "Elements: " << mesh.elements.size() << std::endl;
    // Degree.
    for (const auto &element : mesh.elements)
      degree = (element.degree > degree) ? element.degree : degree;
    std::cout << "Degree: " << degree << std::endl;

    // Update forcing term.
    Vector<Real> F_old = heat.forcing();
    heat.assembly_force(data, mesh);

    // Linear system equation solution.
    Vector<Real> ch = heat.solver(data, mesh, ch_old, F_old);

    if (counter % data.VisualizationStep == 0) {

      // Errors.
      HeatError<Real> error(mesh);

      // Compute error.
      error.error(data, mesh, heat, ch);

      // Output.
      output << "\n" << error << "\n";

      // Compute estimator.
      HeatEstimator<Real> estimator(mesh);
      estimator.computeEstimates(data, heat, ch, ch_old);

      // Determine elements to refine.
      auto [h_mask, p_mask] = estimator.find_elem_to_refine();

      if (mesh.dofs() >= dofsLimit) {
        ch_old = ch;
        ++counter;
        continue;
      }

      // Check if refinement is needed.
      bool refine_p =
          std::any_of(p_mask.begin(), p_mask.end(), [](bool v) { return v; });
      bool refine_h =
          std::any_of(h_mask.begin(), h_mask.end(), [](bool v) { return v; });

      if (refine_p || refine_h) {

        // Initialize ch_old.
        ch_old = ch;

        // Perform p-refinement if necessary
        if (refine_p) {
          // Refine mesh.
          estimator.mesh_refine_degree(p_mask);

          // Prolong solution for p-refinement
          ch_old.resize(estimator.mesh().dofs());
          ch_old =
              heat.prolong_solution_p(estimator.mesh(), mesh, ch_old, p_mask);

          // Update mesh.
          mesh = estimator.mesh();
        }

        // Perform h-refinement if necessary.
        if (refine_h) {
          // Refine mesh.
          estimator.mesh_refine_size(h_mask);

          // Prolong solution for h-refinement.
          ch_old.resize(estimator.mesh().dofs());
          ch_old =
              heat.prolong_solution_h(estimator.mesh(), mesh, ch_old, h_mask);

          // Update mesh.
          mesh = estimator.mesh();
        }

        // Update matrices and forcing term.
        heat.update(data, mesh);

        // Final mesh.
        std::string meshfile =
            "output/heat_hp_" + std::to_string(mesh.elements.size()) + "@" +
            std::to_string(degree) + "_" + std::to_string(heat.t()) + ".poly";
        mesh.write(meshfile, true);
      }
    } else {
      ch_old = ch;
    }
    ++counter;
  }
}