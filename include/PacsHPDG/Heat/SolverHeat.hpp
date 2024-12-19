/**
 * @file SolverHeat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef INCLUDE_PACSHPDG_HEAT_SOLVERHEAT_HPP 
#define INCLUDE_PACSHPDG_HEAT_SOLVERHEAT_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Geometry.hpp"
#include "../Fem.hpp"

namespace pacs {

    // Heat equation solver. 
    Vector<Real> HeatSolver(const Mesh &, const std::array<Sparse<Real>, 3> &, const Vector<Real> &, const std::array<Vector<Real>, 2> &, const Real &, const Real &theta = 0.5, const Real &TOL = 1E-15);    

}

#endif