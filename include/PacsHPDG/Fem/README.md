The `include/PacsHPDG/Fem/` directory contains core finite element method (FEM) functionalities. It provides classes, basis functions, Legendre polynomials, quadrature rules, and utility functions necessary for hp-adaptive Discontinuous Galerkin (DG) methods. 

## Classes

Defines function objects for flexible mathematical operations.

### [`include/PacsHPDG/Fem/Functor.hpp`](./Functor.hpp)

These are variadic template class for wrapping a `std::function<return(values)>`. 

```cpp
template <typename ResultType, typename... Args>  class Functor {};
template <typename ResultType, typename... Args>  class TwoFunctor {};
```

## Methods

Here are presented key mathematical functions used in FEM formulations.

### [`include/PacsHPDG/Fem/Basis.hpp`](./Basis.hpp)

These methods are responsible for the construction of the Lagrange basis functions.

```cpp
std::array<Matrix<Real>, 3> basis_2d(const Mesh &, const std::size_t &, const std::array<Vector<Real>, 2> &);
Matrix<Real> lap_basis_2d(const Mesh &, const std::size_t &, const std::array<Vector<Real>, 2> &);
```

### [`include/PacsHPDG/Fem/Legendre.hpp`](./Legendre.hpp)

These methods compute the Legendre polynomials used to construct the basis functions.

```cpp
Vector<Real> legendre(const Vector<Real> &, const std::size_t &);
Vector<Real> grad_legendre(const Vector<Real> &, const std::size_t &);
Vector<Real> lap_legendre(const Vector<Real> &, const std::size_t &);
```

### [`include/PacsHPDG/Fem/Quadrature.hpp`](./Quadrature.hpp)

These methods are responsible of quadrature rules for numerical integration.

```cpp
// Nodes.
std::array<Vector<Real>, 2> gauss_legendre(const Real &, const Real &, const std::size_t &);

// Reference nodes.
std::array<Vector<Real>, 2> quadrature_1d(const std::size_t &);
std::array<Vector<Real>, 3> quadrature_2d(const std::size_t &);
```

### [`include/PacsHPDG/Fem/Utility.hpp`](./Utility.hpp)

These are utilities methods.

```cpp
// Penalty coefficients.
Vector<Real> penalty(const Mesh &, const std::size_t &, const Real &);

// Get Jaconbian determinant and physical points.
std::tuple<Real, Vector<Real>, Vector<Real>> get_Jacobian_physical_points(const Polygon &, const std::array<Vector<Real>, 2> &);

// Get faces physical points.
std::array<Vector<Real>, 4> faces_physical_points(const Segment &, const Vector<Real> &);
```