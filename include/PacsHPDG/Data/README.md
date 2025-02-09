The `include/PacsHPDG/Data/` directory contains the data structures that define the problem's parameters. These include the domain definition, source term, Dirichlet boundary conditions, time discretization parameters, and space discretization parameters.

## Structs

### [`include/PacsHPDG/Data/DataLaplace.hpp`](./DataLaplace.hpp)

Defines the structure that defines Laplace problem's parameters.

```cpp
struct DataLaplace {};
```
### [`include/PacsHPDG/Data/DataHeat.hpp`](./DataHeat.hpp)

Defines the structure that defines Heat problem's parameters.

```cpp
struct DataHeat {};
```
### [`include/PacsHPDG/Data/DataFisher.hpp`](./DataFisher.hpp)

Defines the structure that defines Fisher-KPP problem's parameters.

```cpp
struct DataFisher {};
```