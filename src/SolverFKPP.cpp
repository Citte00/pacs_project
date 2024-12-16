/**
 * @file SolverFKPP.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

namespace pacs {

    Vector<Real> FKPPsolver(const Mesh &mesh, const size_t &degree, const TriFunctor &alpha, const std::array<Sparse<Real>, 4> &Matrices, const std::array<Vector<Real>, 2> &ch, const Vector<Real> &F_old, const Vector<Real> &F_new, const Real &t, const Real &theta, const Real &TOL) {

        // Mass blocks.
        auto blocks = block_mass(mesh);

        // Extract Matrices from tuple.
        Sparse<Real> M_prj = Matrices[0];
        Sparse<Real> M = Matrices[1];
        Sparse<Real> A = Matrices[2];

        // Extract solutions.
        Vector<Real> c_old = ch[0];
        Vector<Real> c_oold = ch[1];

        // Assembling the constant component of the matrices.
        Sparse<Real> LHS = M_prj + t * theta * A;
        Sparse<Real> RHS = M_prj - t * (1.0 - theta) * A;

        // Assembling the dynamic componenet of the marices.
        if (theta == 0) {
            
            Sparse<Real> MN = NLfisher(mesh, degree, alpha, c_old);
            RHS += t * (M - MN);

        } else if (theta == 1) {
            
            Sparse<Real> MN = NLfisher(mesh, degree, alpha, c_old);
            LHS -= t * (M - MN);
        
        } else if (theta == 0.5) {

            Vector<Real> c_star = 1.5*c_old - 0.5*c_oold;
            Sparse<Real> MN = NLfisher(mesh, degree, alpha, c_star);
            LHS -= 0.5 * t * (M - MN);
            RHS += 0.5 * t * (M - MN);

        }

        // Construction of the complete RHS for the theta method.
        Vector<Real> F = RHS * c_old + t * theta * F_new + t * (1-theta) * F_old;

        // Solves using GMRES.
        return solve(LHS, F);

    }

}