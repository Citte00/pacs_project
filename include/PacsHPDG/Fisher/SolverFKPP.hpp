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

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"

namespace pacs {

    // Custom Fisher-KPP solver.
Vector<Real> FKPPsolver(const DataFKPP &, const Mesh &, const std::array<Sparse<Real>, 4> &,
                        const std::array<Vector<Real>, 2> &,
                        const std::array<Vector<Real>, 2> &,
                        const Real &, const Real &TOL = 1E-15);
}

#endif