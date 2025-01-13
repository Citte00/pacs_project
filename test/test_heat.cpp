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
  pacs::DataHeat Data;

  std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(Data.meshFileSeq);

  // "Splash".
  std::ofstream output{"output/square_hp_" + std::to_string(Data.N) + "@" +
                       std::to_string(Data.degree) + ".error"};
  output << "Square domain - hp-adaptive refinement with estimator." << "\n";

  std::cout << "Square domain - hp-adaptive refinement with estimator."
            << std::endl;
  std::cout << "Output under output/square_hp_" + std::to_string(Data.degree) +
                   ".error."
            << std::endl;

  pacs::Polygon domain{Data.domain};

  // Mesh.
  pacs::Mesh mesh{domain, diagram, Data.degree};

  // Writes mesh informations to a file.
  std::string polyfile = "output/square_" + std::to_string(Data.N) + "@" +
                         std::to_string(Data.degree) + ".poly";
  mesh.write(polyfile);

  // Initializing counter for printing the solution.
  int counter = 1;

  // Build the heat equation object.
  pacs::Heat equation(mesh);

  // Build the heat equation errors object.
  pacs::HeatError errors(mesh);

  // Builds the Fisher-KPP matrices.
  equation.assembly(Data, mesh);

  // Get initial condition.
  equation.evaluateIC(Data, mesh);

  // Compute initial forcing.
  equation.assemblyforce(Data, mesh);

  for (int i = 1; i <= 100*Data.t_f; i++) {

    // Builds the Fisher-KPP matrices.
    // std::array<pacs::Sparse<pacs::Real>, 3> Matrices = pacs::heat(Data,
    // mesh);

    // Get initial condition.
    // pacs::Vector<pacs::Real> ch_old = (counter == 1) ?
    // pacs::EvaluateICHeat(mesh, Matrices[0], Data.c_ex) : pacs::refine(mesh,
    // ch);

    // Compute initial forcing.
    // pacs::Vector<pacs::Real> F_old = (counter == 1) ?
    // pacs::forcingHeat(mesh, D_ext, Source, g_D, Data.t_0) :
    // pacs::forcingHeat(mesh, D_ext, Source, g_D, t - Data.dt);

    equation.t() += Data.dt;
    std::cout << "TIME: " << equation.t() << std::endl;

    // Solve the equation.
    equation.solver(Data, mesh);

    // Errors.
    errors.computeError(mesh, equation, Data.c_ex);

    // Output errors.
    errors.print(output, mesh);

    // Update of the solution.
    equation.ch_old() = equation.ch();

// Solution structure (output).
#ifndef NSOLUTIONS
    if (counter % Data.VisualizationStep == 0) {
      //pacs::Solution solution{Data, mesh, ch, t};
      std::string solfile =
          "output/square_" + std::to_string(mesh.elements.size()) + "@" +
          std::to_string(Data.degree) + "_" + std::to_string(equation.t()) + ".sol";
      //solution.write(solfile);
    }
#endif

    // Estimator.
    // pacs::HeatEstimator estimator{mesh, Matrices[0], ch, ch_old, Source, t,
    // D_ext, c_ex, grad_exact};

    // Refinement.
    // pacs::mesh_refine(mesh, estimator);

    ++counter;
  }
}