# `include/PacsHPDG/Fem/`

## Classes and structs

### [`include/PacsHPDG/Fem/Functor.hpp`](./Functor.hpp)

```cpp
class Functor {};
class TwoFunctor {};
```

## Methods

### [`include/PacsHPDG/Fem/Functor.hpp`](./Functor.hpp)

```cpp
inline Real zero(const Real &, const Real &) { return 0.0; }
```

### [`include/PacsHPDG/Fem/Basis.hpp`](./Basis.hpp)

```cpp
std::array<Matrix<Real>, 3> basis_2d(const Mesh &, const std::size_t &, const std::array<Vector<Real>, 2> &);
Matrix<Real> lap_basis_2d(const Mesh &, const std::size_t &, const std::array<Vector<Real>, 2> &);
```

### [`include/PacsHPDG/Fem/Legendre.hpp`](./Legendre.hpp)

```cpp
Vector<Real> legendre(const Vector<Real> &, const std::size_t &);
Vector<Real> grad_legendre(const Vector<Real> &, const std::size_t &);
Vector<Real> lap_legendre(const Vector<Real> &, const std::size_t &);
```

### [`include/PacsHPDG/Fem/Quadrature.hpp`](./Quadrature.hpp)

```cpp
// Nodes.
std::array<Vector<Real>, 2> gauss_legendre(const Real &, const Real &, const std::size_t &);

// Reference nodes.
std::array<Vector<Real>, 2> quadrature_1d(const std::size_t &);
std::array<Vector<Real>, 3> quadrature_2d(const std::size_t &);
```

### [`include/PacsHPDG/Fem/Utility.hpp`](./Utility.hpp)

```cpp
// Penalty coefficients.
Vector<Real> penalty(const Mesh &, const std::size_t &, const Real &);

// Get Jaconbian determinant and physical points.
std::tuple<Real, Vector<Real>, Vector<Real>> get_Jacobian_physical_points(const Polygon &, const std::array<Vector<Real>, 2> &);

// Get faces physical points.
std::array<Vector<Real>, 4> faces_physical_points(const Segment &, const Vector<Real> &);
```