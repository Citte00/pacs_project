/**
 * @file SolverFKPP.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

namespace pacs {

/**
 * @brief Fisher-KPP equation solver.
 *
 * @param data Data struct.
 * @param mesh Mesh struct.
 * @param Matrices System matrices, [mass, non_lin_mass, stiff, dg_stiff].
 * @param ch Solution vector at previous time step, [c_old, c_oold].
 * @param forcing Forcing temrs, [F_old, F_new].
 * @param TOL Tolerance.
 * @return Vector<Real>
 */
Vector<Real> FKPPsolver(const DataFKPP &data, const Mesh &mesh,
                        const std::array<Sparse<Real>, 4> &Matrices,
                        const std::array<Vector<Real>, 2> &ch,
                        const std::array<Vector<Real>, 2> &forcing,
                        const Real &TOL) {

  // Mass blocks.
  auto blocks = block_mass(mesh);

  // Extract Matrices from tuple.
  auto [M_prj, M, A, DG] = Matrices;

  // Extract forcing term.
  auto [F_old, F_new] = forcing;

  // Extract solutions.
  auto [c_old, c_oold] = ch;

  // Assembling the constant component of the matrices.
  Sparse<Real> LHS = M_prj + data.dt * data.theta * A;
  Sparse<Real> RHS = M_prj - data.dt * (1.0 - data.theta) * A;

  // Assembling the dynamic componenet of the marices.
  if (data.theta == 0) {

    Sparse<Real> MN = NLfisher(data, mesh, c_old);
    RHS += data.dt * (M - MN);

  } else if (data.theta == 1) {

    Sparse<Real> MN = NLfisher(data, mesh, c_old);
    LHS -= data.dt * (M - MN);

  } else if (data.theta == 0.5) {

    Vector<Real> c_star = 1.5 * c_old - 0.5 * c_oold;
    Sparse<Real> MN = NLfisher(data, mesh, c_star);
    LHS -= 0.5 * data.dt * (M - MN);
    RHS += 0.5 * data.dt * (M - MN);
  }

  LHS.compress();
  RHS.compress();

  // Construction of the complete RHS for the theta method.
  Vector<Real> F = RHS * c_old + data.dt * data.theta * F_new + data.dt * (1 - data.theta) * F_old;

  // Solves using GMRES.
  return solve(LHS, F, blocks, GMRES, DBI, TOL);
}
}