/**
 * @file EstimatorsHeat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Error estimator for the heat equation. 
 * @date 2024-12-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef INCLUDE_PACSHPDG_HEAT_ESTIMATORSHEAT_HPP
#define INCLUDE_PACSHPDG_HEAT_ESTIMATORSHEAT_HPP

#include "../Base.hpp"
#include "../Fem.hpp"
#include "../Algebra.hpp"
#include "../Geometry.hpp"

#include <iostream>
#include <array>

namespace pacs {

    struct HeatEstimator {

        // DOFs.
        std::size_t dofs;

        // Estimates.
        Real estimate = 0.0L;
        Vector<Real> estimates;

        // Fits.
        Vector<Real> fits;

        // CONSTRUCTORS.
        HeatEstimator(const Mesh &, const Sparse<Real> &, const Vector<Real> &, const Vector<Real> &, const HeatSource &, const Real &, const TriFunctor &, const TriFunctor &dirichlet = TriFunctor{}, const TriTwoFunctor &dirichlet_gradient = TriTwoFunctor{}, const Real &penalty_coefficient = 10.0);

        // OUTPUT.
        // friend std::ostream &operator <<(std::ostream &, const HeatEstimator &);

    };


}

#endif