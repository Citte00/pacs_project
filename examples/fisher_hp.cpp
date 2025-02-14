/**
 * @file fisher_hp.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Fisher-KPP equation hp-adaptive refinement.
 * @date 2025-02-06
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
  DataFKPP data;

  std::ostringstream oss;
  oss << "output/fisher_hp_" << data.elements << "@" << data.degree;
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
  std::unique_ptr<Fisher> fisher = std::make_unique<Fisher>(data, mesh);
  fisher->assembly(data, mesh);

  // Initial condition.
  Vector<Real> ch_oold = fisher->modal(mesh, data.c_ex);
  fisher->t() += data.dt;
  fisher->ch_old() = fisher->modal(mesh, data.c_ex);

  // Forcing term.
  fisher->assembly_force(data, mesh);

  // Parameters.
  int counter = 1;
  Real t = 0.0;
  const int dofsLimit = DOFS_MAX;
  int steps = static_cast<int>(round(data.t_f / data.dt));
  std::size_t degree = 0;

  for (int i = 1; i <= steps; ++i) {

    // Time step.
    t += data.dt;
    std::cout << "TIME: " << (fisher->t() = t) << std::endl;
    // Number of elements.
    std::cout << "Elements: " << mesh.elements.size() << std::endl;
    // Degree.
    for (const auto &element : mesh.elements)
      degree = (element.degree > degree) ? element.degree : degree;
    std::cout << "Degree: " << degree << std::endl;

    // Update forcing term.
    Vector<Real> F_old = fisher->forcing();
    fisher->assembly_force(data, mesh);

    // Linear system equation solution.
    Vector<Real> ch = fisher->solver(data, mesh, ch_oold, F_old);

    if (counter % data.VisualizationStep == 0) {

      // Errors.
      FisherError error(mesh);

      // Compute error.
      error.error(data, mesh, *fisher, ch);

      // Output.
      output << "\n" << error << "\n";

      // Compute estimates.
      FisherEstimator estimator(mesh);
      estimator.computeEstimates(data, *fisher, ch);

      // Determine elements to refine.
      auto [h_mask, p_mask] = estimator.find_elem_to_refine();

      if (mesh.dofs() >= dofsLimit) {
        ch_oold = fisher->ch_old();
        fisher->ch_old() = ch;
        ++counter;
        continue;
      }

      // Check if refinement is needed.
      bool refine_p =
          std::any_of(p_mask.begin(), p_mask.end(), [](bool v) { return v; });
      bool refine_h =
          std::any_of(h_mask.begin(), h_mask.end(), [](bool v) { return v; });

      if (refine_p || refine_h) {
        // Iniialize new_mesh.
        Mesh new_mesh = mesh;

        // Initialize ch_old.
        ch_oold = fisher->ch_old();
        Vector<Real> ch_old = ch;

        // Perform p-refinement if necessary
        if (refine_p) {

          estimator.mesh_refine_degree(p_mask);
          new_mesh = estimator.mesh();

          // Prolong solution for p-refinement
          ch_oold.resize(new_mesh.dofs());
          ch_oold = fisher->prolong_solution_p(new_mesh, mesh, ch_oold, p_mask);

          ch_old.resize(new_mesh.dofs());
          ch_old = fisher->prolong_solution_p(new_mesh, mesh, ch_old, p_mask);

          mesh = std::move(new_mesh);
        }

        // Perform h-refinement if necessary
        if (refine_h) {
          estimator.mesh_refine_size(h_mask);
          new_mesh = estimator.mesh();

          // Prolong solution for h-refinement
          ch_oold.resize(new_mesh.dofs());
          ch_oold = fisher->prolong_solution_h(new_mesh, mesh, ch_oold, h_mask);

          ch_old.resize(new_mesh.dofs());
          ch_old = fisher->prolong_solution_h(new_mesh, mesh, ch_old, h_mask);

          mesh = std::move(new_mesh);
        }

        // Update matrices.
        fisher = std::make_unique<Fisher>(data, mesh);
        fisher->t() = t;
        fisher->assembly(data, mesh);

        // Update solution.
        fisher->ch_old() = ch_old;

        // Update forcing term.
        fisher->assembly_force(data, mesh);

        // Write mesh file only if refined
        std::string meshfile =
            "output/fisher_hp_" + std::to_string(mesh.elements.size()) + "@" +
            std::to_string(degree) + "_" + std::to_string(t) + ".poly";
        mesh.write(meshfile, true);
      }
    } else {
      ch_oold = fisher->ch_old();
      fisher->ch_old() = ch;
    }
    ++counter;
  }
}
