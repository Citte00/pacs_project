The `include/PacsHPDG/Solvers/` directory provides solver implementations for the Laplace, Heat, and Fisher-KPP equations. These solvers are crucial for numerical simulations in the hp-adaptive Discontinuous Galerkin (DG) method.

## Classes

The class hierarchy follows an inheritance structure, where more complex PDEs extend the functionality of simpler ones. This structure allows code reuse and facilitates the implementation of time-dependent and nonlinear problems.

### [`include/PacsHPDG/Solvers/Laplacian.hpp`](./Laplacian.hpp)

Base class for solving the Laplace equation using the DG method.

```cpp
class Laplace {};
```

### [`include/PacsHPDG/Solvers/Heat.hpp`](./Heat.hpp)

Time-dependent solver that extends Laplace.

```cpp
class Heat : public Laplace {};
```

### [`include/PacsHPDG/Solvers/Fisher.hpp`](./Fisher.hpp)

Nonlinear solver that extends Heat.

```cpp
class Fisher : public Heat {};
```

## Methods

Provides functionality for solving the Laplace, Heat, and Fisher-KPP equations.

### [`include/PacsHPDG/Laplacian/Laplacian.hpp`](./Laplacian.hpp)

Defines methods for system assembly, solver, modal coefficient extraction and mesh-related operations for the Laplace equation.

```cpp
  // Blocks.
  std::vector<std::array<std::vector<std::size_t>, 2>> block_mass(const Mesh &) const;
  // Assembly the laplace system matrix.
  void assembly(const DataLaplace &, const Mesh &);
  // Assembly the forcing term.
  Vector<Real> assembly_force(const DataLaplace &, const Mesh &) const;
  // Solver of the Laplace equation.
  Vector<Real> solver(const Mesh &, const Vector<Real> &,  const Real &TOL = 1E-15) const;
  // Get functions modal coefficients.
  Vector<Real> modal(const Mesh &, const BiFunctor &) const;
```

### [`include/PacsHPDG/Solvers/Heat.hpp`](./Heat.hpp)

Defines methods for system assembly, solver, modal coefficient extraction and hp-adaptive operations for the Heat equation.

```cpp
  // Initialize the methods.
  void initialize(const Mesh &);
  // Assembly the heat equation system matrices.
  void assembly(const DataHeat &, const Mesh &);
  // Assembly the forcing term.
  void assembly_force(const DataHeat &, const Mesh &);
  // Solver of the Heat equation.
  Vector<Real> solver(const DataHeat &, const Mesh &, const Vector<Real> &, const Vector<Real> &, const Real &TOL = 1E-15);
  // Get functions modal coefficients.
  Vector<Real> modal(const Mesh &, const TriFunctor &) const;
  // Get source function modal coefficient.
  Vector<Real> modal_source(const DataHeat &, const Mesh &) const;
  // hp-adaptive methods.
  // Construct matrix with base indeces for each degree.
  Matrix<int> transition(const std::size_t &) const;
  // prolong solution for p.
  Vector<Real> prolong_solution_p(const Mesh &, const Mesh &, const Vector<Real> &, const Mask &) const;
  // prolong solution for h.
  Vector<Real> prolong_solution_h(const Mesh &, const Mesh &, const Vector<Real> &, const Mask &) const
```

### [`include/PacsHPDG/Solvers/Fisher.hpp`](./Fisher.hpp)

Defines methods for system assembly, nonlinear term handling, solver, modal coefficient extraction for the Fisher-KPP equation.

```cpp
  // Assembly the fisher equation system matrices.
  void assembly(const DataFKPP &, const Mesh &);
  // Assembly the non-linear term.
  Sparse<Real> assembly_nl(const DataFKPP &, const Mesh &, const Vector<Real> &);
  // Assembly the forcing term.
  void assembly_force(const DataFKPP &, const Mesh &);
  // Solver of the Heat equation.
  Vector<Real> solver(const DataFKPP &, const Mesh &, const Vector<Real> &, const Vector<Real> &, const Real &TOL = 1E-15);
  // Get source function modal coefficient.
  Vector<Real> modal_source(const DataFKPP &, const Mesh &) const;
```