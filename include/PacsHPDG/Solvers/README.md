The `include/PacsHPDG/Solvers/` directory provides solver implementations for the Laplace, Heat, and Fisher-KPP equations. These solvers are crucial for numerical simulations in the hp-adaptive Discontinuous Galerkin (DG) method.

## Classes

The class hierarchy follows an inheritance structure, where more complex PDEs extend the functionality of simpler ones. This structure allows code reuse and facilitates the implementation of time-dependent and nonlinear problems.

### [`include/PacsHPDG/Solvers/Laplacian.hpp`](./Laplacian.hpp)

Base class for solving the Laplace equation using the DG method.

```cpp
template <NumericType T> class Laplace {};
```

### [`include/PacsHPDG/Solvers/Heat.hpp`](./Heat.hpp)

Time-dependent solver that extends Laplace.

```cpp
template <NumericType T> class Heat : public Laplace<T> {};
```

### [`include/PacsHPDG/Solvers/Fisher.hpp`](./Fisher.hpp)

Nonlinear solver that extends Heat.

```cpp
template <NumericType T> class Fisher : public Heat<T> {};
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
  Vector<T> assembly_force(const DataLaplace &, const Mesh &) const;
  // Solver of the Laplace equation.
  Vector<T> solver(const Mesh &, const Vector<T> &, const Real &TOL = 1E-15) const;
  // Get functions modal coefficients.
  Vector<T> modal(const Mesh &, const Functor<Vector<T>, Vector<T>, Vector<T>> &) const;
```

### [`include/PacsHPDG/Solvers/Heat.hpp`](./Heat.hpp)

Defines methods for system assembly, solver, modal coefficient extraction and hp-adaptive operations for the Heat equation.

```cpp
 // Update matrices and forcing term in adaptive framework.
  void update(const DataHeat &, const Mesh &);
  // Assembly the heat equation system matrices.
  void assembly(const DataHeat &, const Mesh &);
  // Assembly the forcing term.
  void assembly_force(const DataHeat &, const Mesh &);
  // Solver of the Heat equation.
  Vector<T> solver(const DataHeat &, const Mesh &, const Vector<T> &, const Vector<T> &, const T &TOL = 1E-15);
  // Get functions modal coefficients.
  Vector<T> modal(const Mesh &, const Functor<Vector<T>, Vector<T>, Vector<T> , T> &) const;
  // Get source function modal coefficient.
  Vector<T> modal_source(const DataHeat &, const Mesh &) const;
  // hp-adaptive methods.
  // prolong solution for p.
  Vector<T> prolong_solution_p(const Mesh &, const Mesh &, const Vector<T> &, const Mask &) const;
  // prolong solution for h.
  Vector<T> prolong_solution_h(const Mesh &, const Mesh &, const Vector<T> &, const Mask &) const;
```

### [`include/PacsHPDG/Solvers/Fisher.hpp`](./Fisher.hpp)

Defines methods for system assembly, nonlinear term handling, solver, modal coefficient extraction for the Fisher-KPP equation.

```cpp
  // Update matrices and forcing term in adaptive framework.
  void update(const DataFKPP &, const Mesh &);
  // Assembly the fisher equation system matrices.
  void assembly(const DataFKPP &, const Mesh &);
  // Assembly the non-linear term.
  Sparse<T> assembly_nl(const DataFKPP &, const Mesh &, const Vector<T> &);
  // Assembly the forcing term.
  void assembly_force(const DataFKPP &, const Mesh &);
  // Solver of the Heat equation.
  Vector<T> solver(const DataFKPP &, const Mesh &, const Vector<T> &, const Vector<T> &, const T &TOL = 1E-15);
  // Get source function modal coefficient.
  Vector<T> modal_source(const DataFKPP &, const Mesh &) const;
```