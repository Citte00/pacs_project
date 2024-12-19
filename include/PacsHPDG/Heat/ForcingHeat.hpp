/**
 * @file ForcingHeat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef INCLUDE_PACSHPDG_HEAT_FORCINGHEAT_HPP 
#define INCLUDE_PACSHPDG_HEAT_FORCINGHEAT_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Geometry.hpp"
#include "../Fem.hpp"

namespace pacs {

    // Forcing term of heat equation.
    Vector<Real> forcingHeat(const Mesh&, const TriFunctor&, const HeatSource&, const TriFunctor&, const Real&, const Real &penalty_coefficient = 10.0);
}

#endif