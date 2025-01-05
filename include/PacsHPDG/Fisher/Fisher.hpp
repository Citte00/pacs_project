/**
 * @file Fisher.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-11-28
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef INCLUDE_PACSHPDG_FISHER_FISHER_HPP
#define INCLUDE_PACSHPDG_FISHER_FISHER_HPP

#include "../Base.hpp"
#include "../Algebra.hpp"
#include "../Geometry.hpp"
#include "../Fem.hpp"
#include "../Data.hpp"

namespace pacs {

    // Fisher matrices.
    std::array<Sparse<Real>, 4> fisher(const DataFKPP &, const Mesh &);

    // Non-linear Fisher matrix. 
    Sparse<Real> NLfisher(const DataFKPP &, const Mesh &, const Vector<Real> &);

}

#endif