/**
 * @file SolverHeat.cpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief 
 * @date 2024-12-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include <PacsHPDG.hpp>

namespace pacs {

/**
 * @brief Solver of the heat equation matricial system.
 *
 * @param data Data struct.
 * @param mesh Mesh struct.
 * @param Matrices System matrices, [mass, stiff, dg_stiff].
 * @param c_old Soltuion at previous time step.
 * @param forcing Array of forcing terms.
 * @param TOL Tolerance.
 * @return Vector<Real>
 */
Vector<Real> HeatSolver(const DataHeat &data, const Mesh &mesh,
                        const std::array<Sparse<Real>, 3> &Matrices,
                        const Vector<Real> &c_old,
                        const std::array<Vector<Real>, 2> &forcing,
                        const Real &TOL) {
  // Mass blocks.
  Laplace laplacian(mesh);
  auto blocks = laplacian.block_mass(mesh);

  // Extract Matrices from tuple.
  auto [M, A, DG] = Matrices;

  // Extract forcing term.
  auto [F_old, F_new] = forcing;

  // Assembling the constant component of the matrices.
  Sparse<Real> LHS = M + data.dt * data.theta * A;
  Sparse<Real> RHS = M - data.dt * (1.0 - data.theta) * A;

  LHS.compress();
  RHS.compress();

  // Construction of the complete RHS for the theta method.
  Vector<Real> F = RHS * c_old + data.dt * data.theta * F_new +
                   data.dt * (1 - data.theta) * F_old;

  // Solves using GMRES.
  return solve(LHS, F, blocks, GMRES, DBI, TOL);
}
}