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

pacs::Vector<pacs::Real> exact(const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &, const pacs::Real &);
pacs::Vector<pacs::Real> exact_x(const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &, const pacs::Real &);
pacs::Vector<pacs::Real> exact_y(const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &, const pacs::Real &);
pacs::Vector<pacs::Real> exact_t(const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &);
pacs::Vector<pacs::Real> dirichlet(const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &, const pacs::Real &);

pacs::Vector<pacs::Real> D(const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &, const pacs::Real &);
pacs::Vector<pacs::Real> source(const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &, const pacs::Real &, const pacs::Vector<pacs::Real> &);

int main(int argc, char **argv) {

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
        Time = static_cast<pacs::Real>(std::stoi(argv[2]));;

    if(argc == 4) {
        Time = static_cast<pacs::Real>(std::stoi(argv[2]));
        elements = static_cast<std::size_t>(std::stoi(argv[3]));
    }
    
    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram("meshes/square/square_" + std::to_string(elements) + ".poly");

    // "Splash".
    std::ofstream output{"output/square_hp_" + std::to_string(elements) + "@" + std::to_string(degree) + ".error"};
    output << "Square domain - hp-adaptive refinement with estimator." << "\n";

    std::cout << "Square domain - hp-adaptive refinement with estimator." << std::endl;
    std::cout << "Output under output/square_hp_" + std::to_string(degree) + ".error." << std::endl;

    // Domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Mesh.
    pacs::Mesh mesh{domain, diagram, degree};

    // Writes mesh informations to a file.
    std::string polyfile = "output/square_" + std::to_string(elements) + "@" + std::to_string(degree) + ".poly";
    mesh.write(polyfile, true);

    // Functor definition.
    pacs::TriFunctor c_ex(exact);
    pacs::GeneralTwoFunctor<pacs::Vector<pacs::Real>, pacs::Vector<pacs::Real>, pacs::Vector<pacs::Real>, pacs::Real> grad_exact(exact_x, exact_y);
    pacs::BiFunctor dt_exact(exact_t);
    pacs::TriFunctor g_D(dirichlet);
    pacs::TriFunctor D_ext(D);
    pacs::HeatSource Source(source);

    // Builds the Fisher-KPP matrices.
    std::array<pacs::Sparse<pacs::Real>, 3> Matrices = pacs::heat(mesh, degree, D_ext);

    // Get initial condition.
    pacs::Vector<pacs::Real> ch_old = pacs::EvaluateICHeat(mesh, Matrices[0], c_ex);

    // Compute initial forcing.
    pacs::Vector<pacs::Real> F_new = pacs::forcingHeat(mesh, D_ext, Source, g_D, 0.0, degree);

    // Initializing counter for printing the solution.
    int counter = 1;

    for(pacs::Real t = 0.01; t <= Time; t += 0.01) {

        std::cout << "TIME: " << t << std::endl;

        // Update the forcing term.
        pacs::Vector<pacs::Real> F_old = F_new;
        F_new = pacs::forcingHeat(mesh, D_ext, Source, g_D, t, degree);

        pacs::Vector<pacs::Real> ch = pacs::HeatSolver(mesh, degree, Matrices, ch_old, {F_old, F_new}, 0.01);

        // Errors.
        pacs::GeneralError error{mesh, {Matrices[0], Matrices[2]}, ch, c_ex, grad_exact, t};

        // Solution structure (output).

        if (counter % 10 == 0) { 
            pacs::Solution solution{mesh, ch, c_ex, t};
            std::string solfile = "output/square_" + std::to_string(elements) + "@" + std::to_string(degree) + "_" + std::to_string(t) + ".sol";
            solution.write(solfile);
        }

        // Output.
        output << "\n" << error << "\n" << std::endl;

        // Update of the solution.
        ch_old = ch;

        ++counter;

    }
    
}

/**
 * @brief Exact solution.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @return pacs::Real 
 */
inline pacs::Vector<pacs::Real> exact(const pacs::Vector<pacs::Real> &x, const pacs::Vector<pacs::Real> &y, const pacs::Real &t) {
    
    pacs::Vector<pacs::Real> result{x.size()};

    for (size_t i = 0; i < x.size(); i++)
    {
        result[i] = (std::cos(2.0*M_PI*x[i]) * std::cos(2.0*M_PI*y[i]) + 2.0) * (1-t);
    }

    return result;
}

/**
 * @brief Exact solution, x-derivative.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @return pacs::Real 
 */
inline pacs::Vector<pacs::Real> exact_x(const pacs::Vector<pacs::Real> &x, const pacs::Vector<pacs::Real> &y, const pacs::Real &t) {
    
    pacs::Vector<pacs::Real> result{x.size()};

    for (size_t i = 0; i < x.size(); i++)
    {
        result[i] = -2.0L * M_PI * std::sin(2.0*M_PI*x[i]) * std::cos(2.0*M_PI*y[i]) * (1-t);
    }

    return result;
}

/**
 * @brief Exact solution, y-derivative.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @return pacs::Real 
 */
inline pacs::Vector<pacs::Real> exact_y(const pacs::Vector<pacs::Real> &x, const pacs::Vector<pacs::Real> &y, const pacs::Real &t) {
    
    pacs::Vector<pacs::Real> result{x.size()};

    for (size_t i = 0; i < x.size(); i++)
    {
        result[i] = -2.0L * M_PI * std::cos(2.0*M_PI*x[i]) * std::sin(2.0*M_PI*y[i]) * (1-t);
    }

    return result;
}

/**
 * @brief Exact solution, time-derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Vector<pacs::Real> exact_t(const pacs::Vector<pacs::Real> &x, const pacs::Vector<pacs::Real> &y) {
    
    pacs::Vector<pacs::Real> result{x.size()};

    for (size_t i = 0; i < x.size(); i++)
    {
        result[i] = - (std::cos(2.0*M_PI*x[i]) * std::cos(2.0*M_PI*y[i]) + 2.0);
    }

    return result;
}

/**
 * @brief Dirichlet BC.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @return pacs::Real 
 */
inline pacs::Vector<pacs::Real> dirichlet(const pacs::Vector<pacs::Real> &x, const pacs::Vector<pacs::Real> &y, const pacs::Real &t) {
    return exact(x, y, t);
}

/**
 * @brief D_ext definition.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @return pacs::Real 
 */
inline pacs::Vector<pacs::Real> D(const pacs::Vector<pacs::Real> &x, const pacs::Vector<pacs::Real> &y, const pacs::Real &t){
    return 1.0 + 0.0 * x * y * t; 
}

/**
 * @brief Test source.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @param D 
 * @return pacs::Real 
 */
inline pacs::Vector<pacs::Real> source(const pacs::Vector<pacs::Real> &x, const pacs::Vector<pacs::Real> &y, const pacs::Real &t, const pacs::Vector<pacs::Real> &D) 
{
    pacs::Vector<pacs::Real> result{x.size()};

    for (size_t i = 0; i < x.size(); i++)
    {
        result[i] = - (std::cos(2.0*M_PI*x[i])*std::cos(2*M_PI*y[i])+2) + 8*M_PI*M_PI*D[i]*std::cos(2*M_PI*x[i])*std::cos(2*M_PI*y[i])*(1-t);
    }
    
    return result;
}