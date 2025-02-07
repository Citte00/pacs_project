# `include/PacsHPDG/Solvers/`

In this folder the equation classes are built. Each of them containes all the method needed to construct the algebraic system and to solve it then. The classes are built through inheritance.

## Classes and structs

### [`include/PacsHPDG/Solvers/Laplacian.hpp`](./Laplacian.hpp)

```cpp
class Laplace {};
```

### [`include/PacsHPDG/Solvers/Heat.hpp`](./Heat.hpp)

```cpp
class Heat : public Laplace {};
```

### [`include/PacsHPDG/Solvers/Fisher.hpp`](./Fisher.hpp)

```cpp
class Fisher : public Heat {};
```

## Methods

### [`include/PacsHPDG/Laplacian/Laplacian.hpp`](./Laplacian.hpp)

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