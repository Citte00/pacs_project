/**
 * @file Estimators.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-11-27
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ESTIMATORS_FKPP_PACS
#define ESTIMATORS_FKPP_PACS

#include "../Base.hpp"
#include "../Fem.hpp"
#include "../Algebra.hpp"
#include "../Geometry.hpp"

#include <iostream>
#include <array>

namespace pacs {

    /**
     * @brief Error estimator struct for laplacian.
     * 
     */
    struct Estimator {

        // DOFs.
        std::size_t dofs;

        // Estimates.
        Real estimate = 0.0L;
        Vector<Real> estimates;

        // Fits.
        Vector<Real> fits;

        // CONSTRUCTORS.

        Estimator(const Mesh &, const Sparse<Real> &, const Vector<Real> &, const Functor &, const Functor &dirichlet = Functor{}, const TwoFunctor &dirichlet_gradient = TwoFunctor{}, const Real &penalty_coefficient = 10.0);

        // OUTPUT.

        friend std::ostream &operator <<(std::ostream &, const Estimator &);

    };

    /**
     * @brief Error estimator struct for Fisher-KPP.
     * 
     */
    struct EstimatorFKPP {

        // DOFs.
        std::size_t dofs;

        // Estimates.
        Real estimate = 0.0L;
        Vector<Real> estimates;

        // Fits.
        Vector<Real> fits;

        // CONSTRUCTORS.

        EstimatorFKPP(const Mesh &, const Sparse<Real> &, const Vector<Real> &, const Functor3 &, const Functor3 &dirichlet = Functor3{}, const TwoFunctor3 &dirichlet_gradient = TwoFunctor3{}, const Real &penalty_coefficient = 10.0);

        // OUTPUT.

        friend std::ostream &operator <<(std::ostream &, const EstimatorFKPP &);

    };

}

#endif