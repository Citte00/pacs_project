/**
 * @file Fisher.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-11-28
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FISHER_MATRIX_PACS
#define FISHER_MATRIX_PACS

#include "../Base.hpp"
#include "../Algebra.hpp"
#include "../Geometry.hpp"
#include "../Fem.hpp"

namespace pacs {

    // Fisher matrices.
    std::array<Sparse<Real>, 4> fisher(const Mesh &, const TriFunctor&, const TriFunctor&, const Real &penalty_coefficient = 10.0);

    // Non-linear Fisher matrix. 

    Sparse<Real> NLfisher(const Mesh &, const TriFunctor&, const Vector<Real> &, const Real &penalty_coefficient = 10.0);

}

#endif