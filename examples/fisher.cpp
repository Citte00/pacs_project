/**
 * @file fisher.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Fisher equation, uniform refinement.
 * @date 2025-01-22
 *
 * @copyright Copyright (c) 2025
 *
 */
#include <PacsHPDG.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// To save typing the full qualified names.
using namespace pacs;

int main(int argc, char **argv) {

  // Retrieve problem data from structure.
  DataFKPP data;

  std::ostringstream oss;
  oss << "output/fisher_2_" << data.degree;
  std::ofstream output(oss.str() + ".error");

  output << "Square domain - uniform refinement." << "\n";

  std::cout << "Square domain - uniform refinement." << std::endl;
  std::cout << "Output under " << oss.str() << ".error." << std::endl;

  // Domain.
  Polygon domain{data.domain};

  // Diagrams.
  std::vector<std::vector<Polygon>> diagrams;


  diagrams.emplace_back(mesh_diagram("meshes/square/square_4000.poly"));

  // Test.
  for (std::size_t j = 0; j < diagrams.size(); ++j) {

    // Mesh.
    Mesh mesh{domain, std::move(diagrams[j]), data.degree};
    std::cout << "Elements: " << mesh.elements.size() << std::endl;

    // Errors.
    FisherError error(mesh);

    // Matrices.
    Fisher fisher(data, mesh);
    fisher.assembly(data, mesh);

    // Initial condition.
    Vector<Real> ch_oold = fisher.modal(mesh, data.c_ex);
    fisher.t() += data.dt;
    fisher.ch_old() = fisher.modal(mesh, data.c_ex);

    // Forcing term.
    fisher.assembly_force(data, mesh);

    int steps = static_cast<int>(round(data.t_f / data.dt));
    for (int i = 1; i <= steps; i++) {

      // Time step.
      fisher.t() += data.dt;
      std::cout << "TIME: " << fisher.t() << std::endl;

      // Update forcing term.
      Vector<Real> F_old = fisher.forcing();
      fisher.assembly_force(data, mesh);

      // Linear system equation solution.
      Vector<Real> ch = fisher.solver(data, mesh, ch_oold, F_old);

      // Compute error.
      error.error(data, mesh, fisher, ch);

      FisherEstimator estimator(mesh);
      estimator.computeEstimates(data, fisher, ch);
      std::string estimate = "output/square_fisher_" +
                             std::to_string(estimator.mesh().elements.size()) +
                             ".poly";
      estimator.write(estimate, true);

      // Update solution.
      ch_oold = fisher.ch_old();
      fisher.ch_old() = ch;
    }

// Solution structure (output).
#ifndef NSOLUTIONS
    FisherSolution solution{mesh};
    solution.computeSolution(data, mesh, fisher, fisher.ch_old());
    std::string solfile = "output/square_s_" + std::to_string(data.degree) +
                          "_" + std::to_string(j) + ".sol";
    solution.write(solfile);
#endif

    // Output.
    output << "\n"
           << error << "\n"
           << "Residual: "
           << norm(fisher.M() * fisher.ch_old() + fisher.A() * fisher.ch_old() +
                   fisher.M_alpha() * fisher.ch_old() * (1 - fisher.ch_old()) -
                   fisher.forcing())
           << std::endl;
  }
}