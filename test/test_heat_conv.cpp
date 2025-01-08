/**
 * @file test_heat_conv.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2025-01-07
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <PacsHPDG.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

int main(int argc, char **argv) {

  // Retrieve data from struct.
  pacs::DataHeat Data;

  std::ofstream output{"output/square_" + std::to_string(Data.degree) + ".error"};

  output << "Square domain - uniform refinement." << "\n";

  std::cout << "Square domain - uniform refinement." << std::endl;
  std::cout << "Output under output/square_" + std::to_string(Data.degree) +
                   ".error.\n"
            << std::endl;

  // Domain.
  pacs::Polygon domain{Data.domain};

  // Diagrams.
  std::vector<std::vector<pacs::Polygon>> diagrams;

  diagrams.emplace_back(pacs::mesh_diagram("meshes/square/square_125.poly"));
  diagrams.emplace_back(pacs::mesh_diagram("meshes/square/square_250.poly"));
  diagrams.emplace_back(pacs::mesh_diagram("meshes/square/square_500.poly"));
  diagrams.emplace_back(pacs::mesh_diagram("meshes/square/square_1000.poly"));
  diagrams.emplace_back(pacs::mesh_diagram("meshes/square/square_2000.poly"));

  for (std::size_t index = 0; index < diagrams.size(); ++index) {

    // Mesh.
    pacs::Mesh mesh{domain, diagrams[index], Data.degree};

    std::cout << "\nELEMENTS: " << mesh.elements.size() << std::endl;

    // Initializing counter for printing the solution.
    int counter = 1;
    pacs::Real t = 0.0;

    // Builds the Fisher-KPP matrices.
    std::array<pacs::Sparse<pacs::Real>, 3> Matrices = pacs::heat(Data, mesh);

    // Get initial condition.
    pacs::Vector<pacs::Real> ch_old =
        pacs::EvaluateICHeat(mesh, Matrices[0], Data.c_ex);

    // Compute initial forcing.
    pacs::Vector<pacs::Real> F_new = pacs::forcingHeat(Data, mesh, Data.t_0);

    for (int i = 1; i <= 200; i++) {

      t += Data.dt;

      std::cout << "TIME: " << t << std::endl;

      // Update the forcing term.
      pacs::Vector<pacs::Real> F_old = F_new;
      F_new = pacs::forcingHeat(Data, mesh, t);

      pacs::Vector<pacs::Real> ch =
          pacs::HeatSolver(Data, mesh, Matrices, ch_old, {F_old, F_new});

      // Errors.
      pacs::HeatError error{Data, mesh, {Matrices[0], Matrices[2]}, ch, t};

// Solution structure (output).
#ifndef NSOLUTIONS
      if (counter % Data.VisualizationStep == 0) {
        pacs::Solution solution{Data, mesh, ch, t};
        std::string solfile =
            "output/square_" + std::to_string(mesh.elements.size()) + "@" +
            std::to_string(Data.degree) + "_" + std::to_string(t) + ".sol";
        solution.write(solfile);
      }
#endif

      // Output.
      if (i == 200)
        output << "\n" << error << "\n" << std::endl;

      // Update of the solution.
      ch_old = ch;

      ++counter;
    }
  }
}