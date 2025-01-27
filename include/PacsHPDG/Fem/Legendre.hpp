/**
 * @file Legendre.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef INCLUDE_PACSHPDG_FEM_LEGENDRE_HPP
#define INCLUDE_PACSHPDG_FEM_LEGENDRE_HPP

#include "../Base.hpp"
#include "../Algebra.hpp"

namespace pacs {

    // Legendre polynomials and their derivatives.

    Vector<Real> legendre(const Vector<Real> &, const std::size_t &);
    Vector<Real> grad_legendre(const Vector<Real> &, const std::size_t &);
    Vector<Real> lap_legendre(const Vector<Real> &, const std::size_t &);

}

#endif