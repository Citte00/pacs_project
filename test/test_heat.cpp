/**
 * @file test_heat.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief
 * @date 2024-11-26
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <PacsHPDG.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char **argv) {

  // Retrieve data from struct.
  pacs::DataHeat data;

  // "Splash".
  std::ofstream output{"output/square_hp_" + std::to_string(data.N) + "@" +
                       std::to_string(data.degree) + ".error"};
  output << "Square domain - hp-adaptive refinement with estimator." << "\n";

  std::cout << "Square domain - hp-adaptive refinement with estimator."
            << std::endl;
  std::cout << "Output under output/square_hp_" + std::to_string(data.degree) +
                   ".error."
            << std::endl;

  // Domain.
  pacs::Polygon domain{data.domain};

  // Diagram.
  std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(data.meshFileSeq);

  // Mesh.
  pacs::Mesh mesh{domain, diagram, data.degree};

  // Writes mesh informations to a file.
  std::string polyfile = "output/square_" + std::to_string(data.N) + "@" +
                         std::to_string(data.degree) + ".poly";
  mesh.write(polyfile);

  // Initializing counter for printing the solution.
  int counter = 1;

  // Builds the Fisher-KPP matrices.
  pacs::Heat heat(mesh);

  // Builds the Fisher-KPP matrices.
  heat.assembly(data, mesh);

  // Compute initial forcing.
  heat.forcing(data, mesh);

  // Compute initial condition.
  pacs::Vector<pacs::Real> ch_old = heat.evaluateIC(data, mesh);

  for (int i = 1; i <= data.t_f / data.dt; i++) {

    heat.t() += data.dt;
    std::cout << "TIME: " << heat.t() << std::endl;

    // Update the forcing term.
    pacs::Vector<pacs::Real> F_old = heat.forcing();
    heat.forcing(data, mesh);

    // Compute the solution.
    heat.solver(data, mesh, ch_old, F_old);

// Solution structure (output).
#ifndef NSOLUTIONS
    if (counter % data.VisualizationStep == 0) {
      pacs::Solution solution{data, mesh, ch, t};
      std::string solfile =
          "output/square_" + std::to_string(mesh.elements.size()) + "@" +
          std::to_string(data.degree) + "_" + std::to_string(t) + ".sol";
      solution.write(solfile);
    }
#endif

    // Errors.
    if (counter % data.VisualizationStep == 0){
      
      pacs::HeatError error(mesh);
      error.computeErrors(data, mesh, heat);

      // Output.
      output << "\n" << error << "\n" << std::endl;
    }

    // Estimator.
    pacs::HeatEstimator estimator(mesh);
    estimator.computeEstimates(data, heat, heat.ch());

    // Masks.
    pacs::Mask p_mask = estimator.fits() > 1.0;
    pacs::Mask h_mask = (estimator.estimates() * estimator.estimates()) >
                        0.75 *
                            sum(estimator.estimates() * estimator.estimates()) /
                            estimator.mesh().elements.size();

    // Strategy.
    for (std::size_t j = 0; j < estimator.mesh().elements.size(); ++j) {
      if (!h_mask[j]) // p-Refine only error-marked elements.
        p_mask[j] = false;

      if (p_mask[j] && h_mask[j]) // p > h.
        h_mask[j] = false;
    }

    // Refinement.
    estimator.mesh_refine_size(h_mask);
    pacs::Mesh new_mesh = estimator.mesh();
    heat.prolong_solution_h(new_mesh, mesh, h_mask);

    if (counter % data.VisualizationStep == 0){
      std::string polyfile = "output/square_" + std::to_string(data.N) + "@" +
                             std::to_string(heat.t()) + ".poly";
      new_mesh.write(polyfile);
    }

    // Update matrices.
    heat.assembly(data, new_mesh);

    // Update solution.
    auto blocks = heat.block_mass(new_mesh);
    ch_old.resize(new_mesh.dofs());
    ch_old = pacs::solve(heat.M(), heat.ch(), blocks, pacs::DB);

    // Update forcing term.
    heat.forcing().resize(new_mesh.dofs());
    heat.forcing(data, new_mesh);

    // Update mesh.
    mesh = new_mesh;

    ++counter;
  };
}