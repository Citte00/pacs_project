/**
 * @file GeneralErrors.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEN_ERRORS_PACS
#define GEN_ERRORS_PACS

#include "../Base.hpp"
#include "../Fem.hpp"
#include "../Algebra.hpp"
#include "../Geometry.hpp"

#include <iostream>
#include <array>

namespace pacs {

    /**
     * @brief GeneralError struct.
     * 
     */
    struct GeneralError {
        
        std::size_t elements;
        std::size_t dofs;
        std::size_t degree; // p.
        Real size; // h.

        Real dg_error;
        Real l2_error;

        Vector<Real> l2_errors;
        Vector<Real> h1_errors;

        // CONSTRUCTORS.
        GeneralError(const Mesh &, const std::array<Sparse<Real>, 2> &, const Vector<Real> &, const TriFunctor &, const GeneralTwoFunctor<Vector<Real>, Vector<Real>, Vector<Real>, Real> &, const Real &);

        // OUTPUT.
        friend std::ostream &operator <<(std::ostream &, const GeneralError &);

    };

}

#endif