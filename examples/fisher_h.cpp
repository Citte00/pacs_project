/**
 * @file fisher_h.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Fisher-KPP equation element size adaptive refinement.
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
  oss << "output/fisher_h_" << data.elements << "@" << data.degree;
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

  for (int i = 1; i <= steps; ++i) {

    // Time step.
    t += data.dt;
    std::cout << "TIME: " << (fisher->t() = t) << std::endl;
    // Number of elements.
    std::cout << "Elements: " << mesh.elements.size() << std::endl;

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
      output << "\n"
             << error << "\n"
             << "Residual: "
             << norm(fisher->M() * ch + fisher->A() * ch - fisher->forcing())
             << std::endl;

      // Compute estimates.
      FisherEstimator estimator(mesh);
      estimator.computeEstimates(data, *fisher, ch);
      const auto &estimates = estimator.estimates();

      // Refine.
      Mask h_mask = estimates > 0.6L * max(estimates);

      bool refine_h =
          std::any_of(h_mask.begin(), h_mask.end(), [](bool v) { return v; });

      if (!refine_h || mesh.dofs() >= dofsLimit) {
        ch_oold = fisher->ch_old();
        fisher->ch_old() = ch;
        ++counter;
        continue;
      }

      // Refine mesh.
      estimator.mesh_refine_size(h_mask);
      const auto &new_mesh = estimator.mesh();

      ch_oold.resize(new_mesh.dofs());
      ch_oold =
          fisher->prolong_solution_h(new_mesh, mesh, fisher->ch_old(), h_mask);

      // Update matrices.
      fisher = std::make_unique<Fisher>(data, new_mesh);
      fisher->t() = t;
      fisher->assembly(data, new_mesh);

      // Prolong solution.
      fisher->ch_old() = fisher->prolong_solution_h(new_mesh, mesh, ch, h_mask);

      // Update forcing.
      fisher->assembly_force(data, new_mesh);

      // Update mesh.
      mesh = std::move(new_mesh);

      // Final mesh.
      std::string meshfile =
          "output/fisher_h_" + std::to_string(mesh.elements.size()) + "@" +
          std::to_string(data.degree) + "_" + std::to_string(t) + ".poly";
      mesh.write(meshfile);

    } else {
      ch_oold = fisher->ch_old();
      fisher->ch_old() = ch;
    }
    ++counter;
  }
}
