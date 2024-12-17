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

    /**
     * @brief Fisher-KPP equation solver
     * 
     * @param mesh Mesh.
     * @param degree Polynomial degree.
     * @param alpha Non-linear term coefficient. 
     * @param Matrices [M_prj, M, A, DG] system matrices. 
     * @param ch [c_h^n, c_h^(n-1)] solutions at previous time steps.
     * @param F_old forcing term at time step n
     * @param F_new forcing term at time step n+1
     * @param t Time step.
     * @param theta Time discretization parameter (defaulted to 0.5 for CN method).
     * @param TOL tolerance.
     * @return Vector<Real> 
     */
    Vector<Real> FKPPsolver(const Mesh &mesh, const size_t &degree, const TriFunctor &alpha, const std::array<Sparse<Real>, 4> &Matrices, const std::array<Vector<Real>, 2> &ch, const Vector<Real> &F_old, const Vector<Real> &F_new, const Real &t, const Real &theta, const Real &TOL) {

        // Mass blocks.
        auto blocks = block_mass(mesh);

        // Extract Matrices from tuple.
        auto [M_prj, M, A, DG] = Matrices;

        // Extract solutions.
        auto [c_old, c_oold] = ch;

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