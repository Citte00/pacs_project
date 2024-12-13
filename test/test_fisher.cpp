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
pacs::Vector<pacs::Real> alpha(const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &, const pacs::Real &);
pacs::Vector<pacs::Real> source(const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &, const pacs::Real &, const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &);

int main(int argc, char **argv) {

    // Degree.
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
    pacs::TriFunctor exactfunc(exact);
    pacs::GeneralTwoFunctor<pacs::Vector<pacs::Real>, pacs::Vector<pacs::Real>, pacs::Vector<pacs::Real>, pacs::Real> exactfunc2(exact_x, exact_y);
    pacs::BiFunctor exactfunc_t(exact_t);
    pacs::TriFunctor DirBC(dirichlet);
    pacs::TriFunctor alphafunc(alpha);
    pacs::TriFunctor Dextfunc(D);
    pacs::SourceFunctor Source(source);

    // Builds the Fisher-KPP matrices.
    std::array<pacs::Sparse<pacs::Real>, 4> Matrices = pacs::fisher(mesh, alphafunc, Dextfunc);

    // Get initial condition.
    std::array<pacs::Vector<pacs::Real>, 2> ch_old = pacs::getInitialCond(mesh, Matrices[0], exactfunc, 0.0);

    // Compute initial forcing.
    pacs::Vector<pacs::Real> F_new = pacs::forcingFKPP(mesh, alphafunc, Dextfunc, Source, DirBC, 0.0);

    for(pacs::Real t = 0.01; t <= Time; t += 0.01) {

        std::cout << "TIME: " << t << std::endl;

        // Update the forcing term.
        pacs::Vector<pacs::Real> F_old = F_new;
        F_new = pacs::forcingFKPP(mesh, alphafunc, Dextfunc, Source, DirBC, t);

        pacs::Vector<pacs::Real> ch = pacs::FKPPsolver(mesh, alphafunc, Matrices, ch_old, F_old, F_new, t);

        // Errors.
        pacs::GeneralError error{mesh, {Matrices[0], Matrices[3]}, ch, exactfunc, exactfunc2, t};

        // Solution structure (output).
        pacs::Solution solution{mesh, ch, exactfunc, t};
        std::string solfile = "output/square_" + std::to_string(elements) + "@" + std::to_string(degree) + "_" + std::to_string(t) + ".sol";
        solution.write(solfile);

        // Output.
        output << "\n" << error << "\n" << std::endl;

        ch_old[1] = ch_old[0];
        ch_old[0] = ch;

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
 * @brief Definition of parameter in non linear term.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @return pacs::Real 
 */
inline pacs::Vector<pacs::Real> alpha(const pacs::Vector<pacs::Real> &x, const pacs::Vector<pacs::Real> &y, const pacs::Real &t) {
    return 1.0 + 0.0 * x * y * t;
}

/**
 * @brief Test source.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @param D 
 * @param alpha 
 * @return pacs::Real 
 */
inline pacs::Vector<pacs::Real> source(const pacs::Vector<pacs::Real> &x, const pacs::Vector<pacs::Real> &y, const pacs::Real &t, const pacs::Vector<pacs::Real> &D, const pacs::Vector<pacs::Real> &alpha) 
{
    pacs::Vector<pacs::Real> result{x.size()};

    for (size_t i = 0; i < x.size(); i++)
    {
        result[i] = - (std::cos(2.0*M_PI*x[i])*std::cos(2*M_PI*y[i])+2) + 8*M_PI*M_PI*D[i]*std::cos(2*M_PI*x[i])*std::cos(2*M_PI*y[i])*(1-t) - alpha[i]*(1-t)*(std::cos(2*M_PI*x[i])*std::cos(2*M_PI*y[i])+2)*(1-(1-t)*(std::cos(2*M_PI*x[i])*std::cos(2*M_PI*y[i])+2));
    }
    
    return result;
}