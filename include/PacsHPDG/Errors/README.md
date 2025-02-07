# `include/PacsHPDG/Errors/`

This folder contains the Error, Estimator and Solution classes. All these classes are built through inheritance.

## Classes and structs

### [`include/PacsHPDG/Error/Errors.hpp`](./Errors.hpp)

```cpp
class LaplaceError {};
class HeatError {};
class FisherError {};
```

### [`include/PacsHPDG/Error/Estimators.hpp`](./Estimators.hpp)

```cpp
class LaplaceEstimator {};
class HeatEstimator {};
class FisherEstimator {};
```

### [`include/PacsHPDG/Error/Solutions.hpp`](./Solutions.hpp)

```cpp
class LaplaceSolution {};
class HeatSolution {};
class FisherSolution {};
```

## Methods

### [`include/PacsHPDG/Error/Errors.hpp`](./Errors.hpp)

```cpp
  // Compute L2, DG and energy errors.
  void error(const DataFKPP &, const Mesh &, const Fisher &,
             const Vector<Real> &);
  // Compute L2 and H1 errors.
  void errors(const DataFKPP &, const Mesh &, const Fisher &,
              const Vector<Real> &);
  // Friend operator<< for polymorphic printing.
  friend std::ostream &operator<<(std::ostream &ost, const FisherError &error) 
```

### [`include/PacsHPDG/Error/Estimators.hpp`](./Estimators.hpp)

```cpp
// Compute error estimates.
  void computeEstimates(const DataLaplace &, const Laplace &,
                        const Vector<Real> &numerical);
  // Polynomial fit.
  Vector<Real> polyfit(const Vector<Real> &, const Vector<Real> &,
                       const std::size_t &) const;

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

```cpp
  // Save numerical and exact solution for plotting.
  void solution(const DataLaplace &, const Mesh &, const Vector<Real> &);
  // Output solution in a .sol file.
  void write(const std::string &);
```