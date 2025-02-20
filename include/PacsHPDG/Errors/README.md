The `include/PacsHPDG/Errors/` directory provides error computation, estimation, and solution handling functionalities. It contains classes and methods for evaluating numerical accuracy and refining the computational mesh.

## Classes

The class hierarchy follows an inheritance structure, where more complex PDEs extend the functionality of simpler ones. This design allows code reuse and hierarchical refinement when working with error analysis. Only the methods for the Laplace equation are introduced below, since those for the other equations are quite similar.

### [`include/PacsHPDG/Error/Errors.hpp`](./Errors.hpp)

These classes represent error computations for different equations.

```cpp
class LaplaceError {};
class HeatError : public LaplaceError {};
class FisherError : public HeatError {};
```

### [`include/PacsHPDG/Error/Estimators.hpp`](./Estimators.hpp)

These classes implement error estimation techniques for adaptive refinement.

```cpp
class LaplaceEstimator {};
class HeatEstimator : public LaplaceEstimator {};
class FisherEstimator : public HeatEstimator {};
```

### [`include/PacsHPDG/Error/Solutions.hpp`](./Solutions.hpp)

These classes store and manage numerical and exact solutions.

```cpp
class LaplaceSolution {};
class HeatSolution : public LaplaceSolution {};
class FisherSolution : public HeatSolution {};
```

## Methods

Provides functionality for error computation, estimation, and solution output.

### [`include/PacsHPDG/Error/Errors.hpp`](./Errors.hpp)

These are the main methods for errors computation.

```cpp
// Compute L2 and DG errors.
void error(const DataLaplace &, const Mesh &, const Laplace<Real> &, const Vector<Real> &);
// Compute L2 and H1 errors.
void errors(const DataLaplace &, const Mesh &, const Laplace<Real> &, const Vector<Real> &);
// Friend operator<< for printing.
friend std::ostream &operator<<(std::ostream &ost, const LaplaceError &error); 
```

### [`include/PacsHPDG/Error/Estimators.hpp`](./Estimators.hpp)

These are the main methods for error estimation and mesh refinment.

```cpp
// Compute error estimates.
void computeEstimates(const DataLaplace &, const Laplace<Real> &, const Vector<Real> &numerical);
// Polynomial fit.
Vector<Real> polyfit(const Vector<Real> &, const Vector<Real> &, const std::size_t &) const;

// Refinement.
void mesh_refine_size(const Mask &);
void mesh_refine_degree(const Mask &);

// Elements to refine refinement.
std::array<Mask, 2> find_elem_to_refine(const Real &refine = 0.75, const Real &speed = 1.0);

// Friend operator<< for output printing.
friend std::ostream &operator<<(std::ostream &ost, const LaplaceEstimator &estimator);

// Outputs the error estimate distribution.
void write(const std::string &, const bool &estimates = false);
```

### [`include/PacsHPDG/Error/Solutions.hpp`](./Solutions.hpp)

These are the main methods to store and plot the numerical and exact solutions.

```cpp
// Save numerical and exact solution for plotting.
void solution(const DataLaplace &, const Mesh &, const Vector<Real> &);
// Output solution in a .sol file.
void write(const std::string &);
```