/**
 * @file test_fisher.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-11-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#define GAUSS_ORDER

#include <PacsHPDG.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

int main(int argc, char **argv) {

    /* Degree.
    if(argc <= 1) {
        std::cout << "Usage: " << argv[0] << " DEGREE [ELEMENTS] [TIME]." << std::endl;
        std::exit(-1);
    }

    std::size_t degree = static_cast<std::size_t>(std::stoi(argv[1]));

    // Initial diagram.
    std::size_t elements = 300;

    // Fina time.
    pacs::Real Time = 2;

    if(argc == 3)
        elements = static_cast<std::size_t>(std::stoi(argv[2]));

    if(argc == 4) {
        elements = static_cast<std::size_t>(std::stoi(argv[2]));
        Time = static_cast<pacs::Real>(std::stoi(argv[3]));
    }*/

    pacs::DataFKPP Data;

    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(Data.meshFileSeq);

    // "Splash".
    std::ofstream output{"output/square_hp_" + std::to_string(Data.N) + "@" + std::to_string(Data.degree) + ".error"};
    output << "Square domain - hp-adaptive refinement with estimator." << "\n";

    std::cout << "Square domain - hp-adaptive refinement with estimator." << std::endl;
    std::cout << "Output under output/square_hp_" + std::to_string(Data.degree) + ".error." << std::endl;

    // Domain.
    pacs::Polygon domain{Data.domain};

    // Mesh.
    pacs::Mesh mesh{domain, diagram, Data.degree};

    // Writes mesh informations to a file.
    std::string polyfile = "output/square_" + std::to_string(Data.N) + "@" + std::to_string(Data.degree) + ".poly";
    mesh.write(polyfile, true);

    // Builds the Fisher-KPP matrices.
    std::array<pacs::Sparse<pacs::Real>, 4> Matrices = pacs::fisher(Data, mesh);

    // Get initial condition.
    std::array<pacs::Vector<pacs::Real>, 2> ch_old = pacs::EvaluateICFKPP(mesh, Matrices[0], Data.c_ex, Data.dt);

    // Compute initial forcing.
    pacs::Vector<pacs::Real> F_new = pacs::forcingFKPP(Data, mesh, Data.t_0);

    for(pacs::Real t = (Data.t_0 + Data.dt); t <= Data.t_f; t += Data.dt) {

        std::cout << "TIME: " << t << std::endl;

        // Update the forcing term.
        pacs::Vector<pacs::Real> F_old = F_new;
        F_new = pacs::forcingFKPP(Data, mesh, t);

        pacs::Vector<pacs::Real> ch = pacs::FKPPsolver(Data, mesh, Matrices, ch_old, {F_old, F_new}, t);

        // Errors.
        pacs::FKPPError error{Data, mesh, {Matrices[0], Matrices[3]}, ch, t};

        // Solution structure (output).
        pacs::Solution solution{Data, mesh, ch, t};
        std::string solfile = "output/square_" + std::to_string(Data.N) + "@" + std::to_string(Data.degree) + "_" + std::to_string(t) + ".sol";
        solution.write(solfile);

        // Output.
        output << "\n" << error << "\n" << std::endl;

        ch_old[1] = ch_old[0];
        ch_old[0] = ch;

    }
}