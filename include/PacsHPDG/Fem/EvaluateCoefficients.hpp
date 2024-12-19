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

    // Modal coefficients of Fisher-KPP source function.
    Vector<Real> evaluateSourceFKPP(const Mesh &, const FKPPSource &, const Real &, const TriFunctor &, const TriFunctor &);

    // Modal coefficients of Heat source function.
    Vector<Real> evaluateSourceHeat(const Mesh &, const HeatSource &, const Real &, const TriFunctor &);

    // Get Fisher-KPP equation initial condition.
    std::array<Vector<Real>, 2> EvaluateICFKPP(const Mesh &, const Sparse<Real> &, const TriFunctor &, const Real &);

    // Get Heat equation initial condition.
    Vector<Real> EvaluateICHeat(const Mesh &, const Sparse<Real> &, const TriFunctor &);

}

#endif