/**
 * @file Heat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef INCLUDE_PACSHPDG_HEAT_HEAT_HPP
#define INCLUDE_PACSHPDG_HEAT_HEAT_HPP

#include "../Base.hpp"
#include "../Algebra.hpp"
#include "../Geometry.hpp"
#include "../Fem.hpp"

namespace pacs {

    // Heat matrices.
    std::array<Sparse<Real>, 3> heat(const Mesh &, const TriFunctor&, const Real &penalty_coefficient = 10.0);

}


#endif