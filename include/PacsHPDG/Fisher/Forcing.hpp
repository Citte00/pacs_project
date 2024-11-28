/**
 * @file Forcing.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-11-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FORCING_FISHER_PACS
#define FORCING_FISHER_PACS

#include "../Base.hpp"
#include "../Geometry.hpp"
#include "../Fem.hpp"

namespace pacs {

    // RHS.

    Vector<Real> forcingFisher(const Mesh &, const Functor3 &, const Functor3 &dirichlet = Functor3{}, const Real &penalty_coefficient = 10.0, const Real &);
}

#endif