/**
 * @file LaplaceSolution.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Laplace equation readable and plottable solution class object.
 * @date 2025-01-12
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_SOLUTIONS_LAPLACE_SOLUTION_HPP
#define INCLUDE_PACSHPDG_SOLUTIONS_LAPLACE_SOLUTION_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Geometry.hpp"

#include <string>

namespace pacs {
class LaplaceSolution {
protected:
  Vector<Real> m_x;
  Vector<Real> m_y;
  Vector<Real> m_numerical;
  Vector<Real> m_exact;

public:
  // CONSTRUCTOR.
  LaplaceSolution(const Mesh &mesh_)
      : m_x{mesh_.dofs()}, m_y{mesh_.dofs()}, m_numerical{mesh_.dofs()},
        m_exact{mesh_.dofs()} {};
};
} // namespace pacs

#endif