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
#include <memory>
#include <string>

// To save typing the full qualified names.
using namespace pacs;

int main(int argc, char **argv) {

  // Retrieve problem data from structure.
  DataHeat data;

  std::ostringstream oss;
  oss << "output/heat_h_" << data.elements << "@" << data.degree;
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

  // Sequence of meshes.
  for (std::size_t index = 0; index < TESTS_MAX; ++index) {

    // Verbosity.
    std::cout << "\nDEGREE: " << data.degree << "\nINDEX: " << index << "\n"
              << std::endl;

    // Matrices.
    Heat<Real> heat(mesh);
    heat.assembly(data, mesh);

    // Initial condition.
    Vector<Real> ch_old = heat.modal(mesh, data.c_ex);

    // Forcing term.
    heat.assembly_force(data, mesh);

    // Parameters.
    Real t = 0.0;
    const int dofsLimit = DOFS_MAX;
    int steps = static_cast<int>(round((data.t_f - data.t_0) / data.dt));

    for (int i = 1; i <= steps; ++i) {

      // Time step.
      heat.t() += data.dt;
      std::cout << "TIME: " << heat.t() << std::endl;
      // Number of elements.
      std::cout << "Elements: " << mesh.elements.size() << std::endl;

      // Update forcing term.
      Vector<Real> F_old = heat.forcing();
      heat.assembly_force(data, mesh);

      // Linear system equation solution.
      Vector<Real> ch = heat.solver(data, mesh, ch_old, F_old);

      // Errors.
      HeatError<Real> error(mesh);

      // Compute error.
      error.error(data, mesh, heat, ch);

      // Output.
      output << "\n" << error << "\n";

      // Exit.
      if (error.dofs() > dofsLimit)
        break;

      // Compute estimates.
      HeatEstimator<Real> estimator(mesh);
      estimator.computeEstimates(data, heat, ch, ch_old);
      const auto &estimates = estimator.estimates();

      // Refine.
      Mask h_mask = estimates > 0.75L * max(estimates);

      bool refine_h =
          std::any_of(h_mask.begin(), h_mask.end(), [](bool v) { return v; });

      // Refine mesh.
      estimator.mesh_refine_size(h_mask);
      const auto &new_mesh = estimator.mesh();

      // Update mesh.
      mesh = std::move(new_mesh);

      // Final mesh.
      std::string meshfile =
          "output/heat_h_" + std::to_string(mesh.elements.size()) + "@" +
          std::to_string(data.degree) + "_" + std::to_string(t) + ".poly";
      mesh.write(meshfile);
    }
  }
}
