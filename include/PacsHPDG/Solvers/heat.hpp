/**
 * @file heat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Class for the heat equation.
 * @date 2024-12-28
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef INCLUDE_PACSHPDG_SOLVERS_HEAT_HPP
#define INCLUDE_PACSHPDG_SOLVERS_HEAT_HPP

#include "./laplacian.hpp"

namespace pacs {

class Heat : public Laplace {
  protected:
    // Time istant.
    Real m_t;

  public:
    // CONSTRUCTOR.
    Heat(const Mesh &mesh_) : Laplace{mesh_}, m_t{0.0} {};

    // GETTER.
    Real t() const { return m_t; };
    Real &t() { return m_t; };

    // METHODS.
    // Assembly of the system matrices.
    void assembly(const DataHeat &, const Mesh &);

    // Assembly of the forcing term.
    void assemblyforce(const DataHeat &, const Mesh &);

    // Solve the heat equation system.
    void solver(const DataHeat &, const Mesh &, const Real &TOL = 1E-15);

    // Get coefficient of a function.
    Vector<Real> evaluateCoeff(const Mesh &, const TriFunctor &, const Real &) const;

    // Get coefficients of source function.
    Vector<Real> evaluateSource(const DataHeat &, const Mesh &, const Real &) const;

    // Get initial condition.
    void evaluateIC(const DataHeat &, const Mesh &);

    // Prolong the solution to the finer mesh.
    // void prolong_solution(const Mesh &, const Mesh &, const Mask &, const Mask &);

    // Compute transition vector.
    // Vector<Real> construct_transition(const Mesh &);
};

}

#endif