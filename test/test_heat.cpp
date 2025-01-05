/**
 * @file test_fisher.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-11-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include <PacsHPDG.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

int main(int argc, char **argv) {

    /*
    // Degree.
    if(argc <= 1) {
        std::cout << "Usage: " << argv[0] << " DEGREE [TIME] [ELEMENTS]." << std::endl;
        std::exit(-1);
    }

    std::size_t degree = static_cast<std::size_t>(std::stoi(argv[1]));

    // Initial diagram.
    std::size_t elements = 300;

    // Fina time.
    pacs::Real Time = 2;

    if(argc == 3)
        Time = static_cast<pacs::Real>(std::stod(argv[2]));;

    if(argc == 4) {
        Time = static_cast<pacs::Real>(std::stod(argv[2]));
        elements = static_cast<std::size_t>(std::stoi(argv[3]));
    }*/
    
    // Retrieve data from struct.
    pacs::DataHeat Data;

    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(Data.meshFileSeq);

    // "Splash".
    std::ofstream output{"output/square_hp_" + std::to_string(Data.N) + "@" + std::to_string(Data.degree) + ".error"};
    output << "Square domain - hp-adaptive refinement with estimator." << "\n";

    std::cout << "Square domain - hp-adaptive refinement with estimator." << std::endl;
    std::cout << "Output under output/square_hp_" + std::to_string(Data.degree) + ".error." << std::endl;

    pacs::Polygon domain{Data.domain};

    // Mesh.
    pacs::Mesh mesh{domain, diagram, Data.degree};

    // Writes mesh informations to a file.
    std::string polyfile = "output/square_" + std::to_string(Data.N) + "@" + std::to_string(Data.degree) + ".poly";
    mesh.write(polyfile);

    // Initializing counter for printing the solution.
    int counter = 1;

    // Builds the Fisher-KPP matrices.
    std::array<pacs::Sparse<pacs::Real>, 3> Matrices = pacs::heat(Data, mesh);
    std::ofstream matrixfile{"output/square_" + std::to_string(Data.N) + "@" + std::to_string(Data.degree) + ".mat"};
    matrixfile << Matrices[1];

    // Get initial condition.
    pacs::Vector<pacs::Real> ch_old = pacs::EvaluateICHeat(mesh, Matrices[0], Data.c_ex);

    // Compute initial forcing.
    pacs::Vector<pacs::Real> F_new = pacs::forcingHeat(Data, mesh, Data.t_0);

    for(pacs::Real t = (Data.t_0 + Data.dt); t <= Data.t_f; t += Data.dt) {

        // Builds the Fisher-KPP matrices.
        //std::array<pacs::Sparse<pacs::Real>, 3> Matrices = pacs::heat(Data, mesh);

        // Get initial condition.
        //pacs::Vector<pacs::Real> ch_old = (counter == 1) ? pacs::EvaluateICHeat(mesh, Matrices[0], Data.c_ex) : pacs::refine(mesh, ch);

        // Compute initial forcing.
        //pacs::Vector<pacs::Real> F_old = (counter == 1) ? pacs::forcingHeat(mesh, D_ext, Source, g_D, Data.t_0) : pacs::forcingHeat(mesh, D_ext, Source, g_D, t - Data.dt);

        std::cout << "TIME: " << t << std::endl;

        // Update the forcing term.
        pacs::Vector<pacs::Real> F_old = F_new;
        F_new = pacs::forcingHeat(Data, mesh, t);

        pacs::Vector<pacs::Real> ch = pacs::HeatSolver(Data, mesh, Matrices, ch_old, {F_old, F_new});

        // Errors.
        pacs::HeatError error{Data, mesh, {Matrices[0], Matrices[2]}, ch, t};

        // Solution structure (output).
        #ifndef NSOLUTIONS
        if (counter % Data.VisualizationStep == 0) { 
            pacs::Solution solution{Data, mesh, ch, t};
            std::string solfile = "output/square_" + std::to_string(mesh.elements.size()) + "@" + std::to_string(Data.degree) + "_" + std::to_string(t) + ".sol";
            solution.write(solfile);
        }
        #endif

        // Output.
        output << "\n" << error << "\n" << std::endl;

        // Estimator.
        //pacs::HeatEstimator estimator{mesh, Matrices[0], ch, ch_old, Source, t, D_ext, c_ex, grad_exact};

        // Refinement.
       // pacs::mesh_refine(mesh, estimator);

        // Update of the solution.
        //ch_old = pacs::refine(ch);
        ch_old = ch;

        ++counter;

    }
    
}