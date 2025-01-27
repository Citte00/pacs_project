/**
 * @file LaplaceSolution.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Equations readable and plottable solution class object.
 * @date 2025-01-12
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_ERRORS_SOLUTIONS_HPP
#define INCLUDE_PACSHPDG_ERRORS_SOLUTIONS_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Geometry.hpp"
#include "../Solvers.hpp"

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
      : m_x{mesh_.entries}, m_y{mesh_.entries}, m_numerical{mesh_.entries},
        m_exact{mesh_.entries} {};

  // GETTERS.
  Vector<Real> x() const { return this->m_x; };
  Vector<Real> y() const { return this->m_y; };
  Vector<Real> numerical() const { return this->m_numerical; };
  Vector<Real> exact() const { return this->m_exact; };

  // METHODS.
  // Save numerical and exact solution for plotting.
  void computeSolution(const DataLaplace &, const Mesh &, const Vector<Real> &);
  // Output solution in a .sol file.
  void write(const std::string &);
};

class HeatSolution : public LaplaceSolution {
  public:
  // CONSTRUCTOR.
    HeatSolution(const Mesh &mesh_) : LaplaceSolution(mesh_) {};

    // METHODS.
    void computeSolution(const DataHeat &, const Mesh &,
                         const Heat &);
};

class FisherSolution : public HeatSolution {
public:
  // CONSTRUCTOR.
  FisherSolution(const Mesh &mesh_) : HeatSolution(mesh_) {};

  // METHODS.
  void computeSolution(const DataFKPP &, const Mesh &, const Fisher &);
};

} // namespace pacs

#endif