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

#include <iostream>

pacs::Real exact(const pacs::Real &, const pacs::Real &, const pacs::Real &);
pacs::Real exact_x(const pacs::Real &, const pacs::Real &, const pacs::Real &);
pacs::Real exact_y(const pacs::Real &, const pacs::Real &, const pacs::Real &);
pacs::Real exact_t(const pacs::Real &, const pacs::Real &);
pacs::Real source(const pacs::Real &, const pacs::Real &, const pacs::Real &);
pacs::Real dirichlet(const pacs::Real &, const pacs::Real &, const pacs::Real &);
pacs::Vector<pacs::Real> D(const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &, const pacs::Real &);
pacs::Vector<pacs::Real> alpha(const pacs::Vector<pacs::Real> &, const pacs::Vector<pacs::Real> &, const pacs::Real &);

int main() {

    // Loads a mesh.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};
    pacs::Mesh mesh{domain, pacs::mesh_diagram("meshes/square/square_100.poly"), 3};

    // Writes mesh informations to a file.
    mesh.write("output/fisher.poly");

    pacs::GeneralFunctor<pacs::Vector<pacs::Real>, pacs::Vector<pacs::Real>, pacs::Vector<pacs::Real>, pacs::Real> alphafunc(alpha);
    pacs::GeneralFunctor<pacs::Vector<pacs::Real>, pacs::Vector<pacs::Real>, pacs::Vector<pacs::Real>, pacs::Real> Dextfunc(D);

    // Builds the Fisher-KPP matrices.
    auto [M_prj, M, K, DG] = pacs::fisher(mesh, alphafunc, Dextfunc);

    // Builds the forcing term.
    //pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, );

    // Linear system solution.
    //pacs::Vector<pacs::Real> numerical = pacs::solve(laplacian, forcing, pacs::BICGSTAB);

    // Errors.
    //pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact, {exact_x, exact_y}};

    // Solution structure (output).
    //pacs::Solution solution{mesh, numerical, exact};
    //solution.write("output/poisson.sol");

    // Output.
    //std::cout << "\n" << error << "\n" << std::endl;

}

/**
 * @brief Exact solution.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @return pacs::Real 
 */
inline pacs::Real exact(const pacs::Real &x, const pacs::Real &y, const pacs::Real &t) {
    return (std::cos(2.0*M_PI*x) * std::cos(2.0*M_PI*y) + 2.0) * (1-t);
}

/**
 * @brief Exact solution, x-derivative.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @return pacs::Real 
 */
inline pacs::Real exact_x(const pacs::Real &x, const pacs::Real &y, const pacs::Real &t) {
    return -2.0 * M_PI * std::sin(2.0*M_PI*x) * std::cos(2.0L*M_PI*y) * (1-t);
}

/**
 * @brief Exact solution, y-derivative.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @return pacs::Real 
 */
inline pacs::Real exact_y(const pacs::Real &x, const pacs::Real &y, const pacs::Real &t) {
    return -2.0 * M_PI * std::cos(2.0*M_PI*x) * std::sin(2.0L*M_PI*y) * (1-t);
}

/**
 * @brief Exact solution, time-derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real exact_t(const pacs::Real &x, const pacs::Real &y) {
    return - (std::cos(2.0*M_PI*x) * std::sin(2.0L*M_PI*y) + 2);
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
inline pacs::Real source(const pacs::Real &x, const pacs::Real &y, const pacs::Real &t) {
    return - (std::cos(2.0*M_PI*x)*std::cos(2*M_PI*y)+2) + 8*M_PI*M_PI*D(x,y,t)*std::cos(2*M_PI*x)*std::cos(2*M_PI*y)*(1-t) - alpha(x,y,t)*(1-t)*(std::cos(2*M_PI*x)*std::cos(2*M_PI*y)+2)*(1-(1-t)*(std::cos(2*M_PI*x)*std::cos(2*M_PI*y)+2));
}

/**
 * @brief Dirichlet BC.
 * 
 * @param x 
 * @param y 
 * @param t 
 * @return pacs::Real 
 */
inline pacs::Real dirichlet(const pacs::Real &x, const pacs::Real &y, const pacs::Real &t) {
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