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
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"

namespace pacs {

    // Forcing term of heat equation.
Vector<Real> forcingHeat(const DataHeat &, const Mesh &, const Real &);
}

#endif