/**
 * @file Quadrature.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef INCLUDE_PACSHPDG_FEM_QUADRATURE_HPP
#define INCLUDE_PACSHPDG_FEM_QUADRATURE_HPP

#include "../Base.hpp"
#include "../Algebra.hpp"
#include "../Geometry.hpp"

#include <vector>

namespace pacs {

    // Gauss-Legendre quadrature nodes.

    std::array<Vector<Real>, 2> gauss_legendre(const Real &, const Real &, const std::size_t &);

    // Reference interval and square quadrature nodes.

    std::array<Vector<Real>, 2> quadrature_1d(const std::size_t &);
    std::array<Vector<Real>, 3> quadrature_2d(const std::size_t &);

}

#endif