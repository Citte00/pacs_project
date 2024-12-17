/**
 * @file SolverHeat.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include <PacsHPDG.hpp>

namespace pacs {

    Vector<Real> HeatSolver(const Mesh &mesh, const size_t &degree, const std::array<Sparse<Real>, 3> &Matrices, const Vector<Real> &c_h, const std::array<Vector<Real>, 2> &forcing, const Real &t, const Real &theta, const Real &TOL) {
        // Mass blocks.
        auto blocks = block_mass(mesh);

        // Extract Matrices from tuple.
        auto [M, A, DG] = Matrices;

        // Extract forcing term.
        auto [F_old, F_new] = forcing;

        // Assembling the constant component of the matrices.
        Sparse<Real> LHS = M + t * theta * A;
        Sparse<Real> RHS = M - t * (1.0 - theta) * A;

        // Construction of the complete RHS for the theta method.
        Vector<Real> F = RHS * c_h + t * theta * F_new + t * (1-theta) * F_old;

        // Solves using GMRES.
        return solve(LHS, F);
    }    

}