/**
 * @file EvaluateSolution.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef COEFF_PACS
#define COEFF_PACS

#include "../Base.hpp"
#include "../Geometry.hpp"
#include "../Algebra.hpp"

#include "./GeneralFunctor.hpp"

namespace pacs {

    // Modal coefficients of a function.
    Vector<Real> evaluateCoeff(const Mesh &, const TriFunctor &, const Real &);

    // Modal coefficients of source function.
    Vector<Real> evaluateSource(const Mesh &, const SourceFunctor &, const Real &, const TriFunctor &, const TriFunctor &);

    // Get Initial Condition.
    std::array<Vector<Real>, 2> getInitialCond(const Mesh &, const Sparse<Real> &, const TriFunctor &, const Real &);

}

#endif