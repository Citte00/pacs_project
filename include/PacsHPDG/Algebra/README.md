The `include/PacsHPDG/Algebra/` directory contains fundamental algebraic structures and methods used for numerical computations. It provides vector, matrix, and sparse matrix representations, along with solvers and mathematical utilities for efficient linear algebra operations.

## Structs

Here are presented the core algebraic data structures.

### [`include/PacsHPDG/Algebra/Vector.hpp`](./Vector.hpp)

Defines a generic representation of mathematical vectors that support numerical operations.

```cpp
template<NumericType T> struct Vector {};
```

### [`include/PacsHPDG/Algebra/Matrix.hpp`](./Matrix.hpp)

Defines a generic dense matrix representation with built-in algebraic functionalities.

```cpp
template<NumericType T> struct Matrix {};
```

### [`include/PacsHPDG/Algebra/Sparse.hpp`](./Sparse.hpp)

Defines a generic sparse matrix representation for efficient storage and computation in large-scale problems.

```cpp
template<NumericType T> struct Sparse {};
```

## Methods

Here are presented the mathematical operations, solvers, and utility functions for algebraic structures.

### [`include/PacsHPDG/Algebra/Solvers.hpp`](./Solvers.hpp)

Implements dense and sparse solvers, including direct and iterative methods for solving linear systems. It also provides preconditioners to improve the convergence of iterative solvers.

```cpp
// Solvers and preconditioners (enums).
enum DenseSolver {LUD, QRD};
enum SparseSolver {GMRES, CGM, BICGSTAB};
enum SparseDSolver {DB};
enum Preconditioner {DI, DBI};

// Matrix solve wrappers.
template<NumericType T> Vector<T> solve(const Matrix<T> &, const Vector<T> &, const DenseSolver &S = QRD);
template<NumericType T> Matrix<T> solve(const Matrix<T> &, const Matrix<T> &, const DenseSolver &S = QRD);

// Sparse solve wrappers.
template<NumericType T> Vector<T> solve(const Sparse<T> &, const Vector<T> &, const SparseSolver &S = GMRES, const Real &TOL = 1E-8);
template<NumericType T> Vector<T> solve(const Sparse<T> &, const Vector<T> &, std::vector<std::array<std::vector<std::size_t>, 2>> &, const SparseSolver &S = GMRES, const Preconditioner &P = DBI, const Real &TOL = 1E-8);
template<NumericType T> Vector<T> solve(const Sparse<T> &, const Vector<T> &, std::vector<std::array<std::vector<std::size_t>, 2>> &, const SparseDSolver &S = DB);

// Matrix solvers.
template<NumericType T> Vector<T> _lu(const Matrix<T> &, const Vector<T> &);
template<NumericType T> Matrix<T> _lu(const Matrix<T> &, const Matrix<T> &);
template<NumericType T> Vector<T> _qr(const Matrix<T> &, const Vector<T> &);
template<NumericType T> Matrix<T> _qr(const Matrix<T> &, const Matrix<T> &);

// Sparse iterative solvers.
template<NumericType T> Vector<T> _gmres(const Sparse<T> &, const Vector<T> &, const Real &TOL = 1E-8);
template<NumericType T> Vector<T> _gmres(const Sparse<T> &, const Vector<T> &, const Vector<T> &guess, const Real &TOL = 1E-8);
template<NumericType T> Vector<T> _cgm(const Sparse<T> &, const Vector<T> &, const Real &TOL = 1E-8);
template<NumericType T> Vector<T> _bicgstab(const Sparse<T> &, const Vector<T> &, const Real &TOL = 1E-8);

// Sparse direct solvers.
template<NumericType T> Vector<T> _db(const Sparse<T> &, const Vector<T> &, const std::vector<std::array<std::vector<std::size_t>, 2>> &);

// Preconditioners.
template<NumericType T> Sparse<T> _di(const Sparse<T> &);
template<NumericType T> Sparse<T> _dbi(const Sparse<T> &, const std::vector<std::array<std::vector<std::size_t>, 2>> &);
```

### [`include/PacsHPDG/Algebra/Methods/Vector.hpp`](./Methods/Vector.hpp)

Provides mathematical and utility functions for `Vector<T>`.

```cpp
// Math.
template<NumericType T> inline T min(const Vector<T> &);
template<NumericType T> inline T max(const Vector<T> &);
template<NumericType T> inline T sum(const Vector<T> &);
template<NumericType T> inline T dot(const Vector<T> &, const Vector<T> &);
template<NumericType T> inline Real norm(const Vector<T> &);

namespace std { template<pacs::NumericType T> pacs::Vector<T> abs(const pacs::Vector<T> vector); }

// Tools.
template<NumericType T> Vector<T> stack(const Vector<T> &, const Vector<T> &);
template<NumericType T> Vector<T> flip(const Vector<T> &);
template<NumericType T> Mask highest(const Vector<T> &, const std::size_t &);
template<NumericType T> Mask lowest(const Vector<T> &, const std::size_t &);
```

### [`include/PacsHPDG/Algebra/Methods/Matrix.hpp`](./Methods/Matrix.hpp)

Implememnts decomposition techniques and utilities for `Matrix<T>`.

```cpp
// Math.
template<NumericType T> std::array<Matrix<T>, 2> LU(const Matrix<T> &);
template<NumericType T> std::array<Matrix<T>, 2> QR(const Matrix<T> &);
template<NumericType T> inline Matrix<T> identity(const std::size_t &);
template<NumericType T> T mtrace(const Matrix<T> &);

// Tools.
template<NumericType T> Vector<T> squash(const Matrix<T> &);
```

### [`include/PacsHPDG/Algebra/Methods/Sparse.hpp`](./Methods/Sparse.hpp)

Provide mathematical operations for `Sparse<T>`.

```cpp
// Math.
template<NumericType T> T mtrace(const Sparse<T> &);
```