/**
 * @file Jacobian.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief The function constructs and evaluates the Jacobian of the
 * transformation in the physical points.
 * @date 2025-01-24
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_FEM_JACOBIAN_HPP
#define INCLUDE_PACSHPDG_FEM_JACOBIAN_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Geometry.hpp"

namespace pacs {

// Penalty coefficients.
Vector<Real> penalty(const Mesh &, const std::size_t &, const Real &);

// Get Jaconbian determinant and physical points.
std::tuple<Real, Vector<Real>, Vector<Real>>
get_Jacobian_physical_points(const Polygon &,
                             const std::array<Vector<Real>, 2> &);

// Get faces physical points.
std::array<Vector<Real>, 4>
faces_physical_points(const Segment &, const Vector<Real> &);

}; // namespace pacs

#endif
