/**
 * @file SolverFKPP.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FKPP_SOLVERS_PACS
#define FKPP_SOLVERS_PACS

#include "../Base.hpp"
#include "../Algebra.hpp"

#include "./Fisher.hpp"

namespace pacs {

    // Custom Fisher-KPP solver.
    Vector<Real> FKPPsolver(const Mesh &, const TriFunctor &, const std::array<Sparse<Real>, 4> &, const std::array<Vector<Real>, 2> &, const Vector<Real> &, const Vector<Real> &, const Real &, const Real &theta=0.5, const Real &TOL = 1E-15);

}

#endif