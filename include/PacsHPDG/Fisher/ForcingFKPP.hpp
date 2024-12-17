/**
 * @file Forcing.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-11-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef INCLUDE_PACSHPDG_FISHER_FORCINGFKPP_HPP
#define INCLUDE_PACSHPDG_FISHER_FORCINGFKPP_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Geometry.hpp"
#include "../Fem.hpp"

namespace pacs {

    // Fisher-KPP euation solver.
    Vector<Real> forcingFKPP(const Mesh &, const TriFunctor &, const TriFunctor &, const FKPPSource &, const TriFunctor &, const Real &, const Real &penalty_coefficient = 10.0);
}

#endif