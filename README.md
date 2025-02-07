# Advanced Programming for Scientific Computing - Project

_The hp-Adaptive Discontinuous Galërkin Method_

## Introduction

This repository provides an implementation of the _hp-adaptive_ discontinuous Galërkin method for the Poisson problem, the Heat equation and the Fisher-KPP equation.

:warning: Refer to the [report](#compiling-the-report) for a detailed explanation of the background and results related to this project.

## Table of Contents

- [Introduction](#introduction)
- [Table of Contents](#table-of-contents)
- [Overview](#overview)
- [Setup](#setup)
    - [Cloning the Repository](#cloning-the-repository)
    - [**Compilation and Execution**](#compilation-and-execution)
        - [Compilation Flags](#compilation-flags)
    - [Compiling the Code into a Library](#compiling-the-code-into-a-library)
- [Using the Reporitory](#using-the-repository)
    - [Domains](#domains)
    - [**Examples**](#examples)
    - [Scripts](#scripts)
- [Using the Code](#using-the-code)
    - [Generating a Mesh](#generating-a-mesh)
    - [Solving the Poisson Problem](#solving-the-poisson-problem)
    - [Mesh Refinement](#mesh-refinement)
- [**Notes to the Reader**](#notes-to-the-reader)
    - [On the Implementation of Basic Objects](#on-the-implementation-of-basic-objects)
    - [On the Adaptation from **lymph**](#on-the-adaptation-from-lymph)
    - [On the Examples Structure](#on-the-examples-structure)
    - [On the Custom Laplacian Solver](#on-the-custom-laplacian-solver)
    - [On Parallelism](#on-parallelism)

:warning: Make sure to take a look at [Notes to the Reader](#notes-to-the-reader) as they provide insight into some design choices about the code.

## Overview

The key components are as follows:

- `domains/`: Stores sample meshes generation scripts.
- `examples/`: Provides the main examples of using the repository.
- `format/`: Containes scripts that set and control the formatting of the code.
- `include/`: Holds definitions for the structures and methods utilized in the repository.
    - [`include/PacsHPDG/Algebra/`](./include/PacsHPDG/Algebra/): Structures and methods for vectors, matrices and linear solvers.
    - [`include/PacsHPDG/Data/`](./include/PacsHPDG/Data/): Structures with problems data.
    - [`include/PacsHPDG/Errors/`](./include/PacsHPDG/Errors/): Classes and methods for errors, error estimators computation and solution plotting.
    - [`include/PacsHPDG/Geometry/`](./include/PacsHPDG/Geometry/): Tools for working with polygons and meshes.
    - [`include/PacsHPDG/Fem/`](./include/PacsHPDG/Fem/): Finite element structures and methods.
    - [`include/PacsHPDG/Solvers/`](./include/PacsHPDG/Solvers/): Classes and methods for the numerical implementation and solution of each problems.
- `meshes/`: Includes sample meshes for simple domains.
- `presentation/`: Contains the LaTeX presentation on the project.
- `report/`: Contains the LaTeX report on the project.
- `scripts/`: Contains Python scripts for meshes, solutions and errors visualization.
- `snippets/`: Simpler examples for the report.
- `src/`: Contains the primary implementations for the repository’s structures and methods.
- `templates/`: Contains TikZ templates.
- `test/`: Contains scripts for testing fundamental features.

Every directory under `include/PacsHPDG/` has a `README.md` that lists the classes, structures, and methods introduced in that particular category.

## Setup

### Cloning the Repository

Clone the repository from [here](https://github.com/Citte00/pacs_project.git):

```bash
git clone git@github.com:Citte00/pacs_project.git
```

### Compilation and Execution

The code is written to work with the C++ standard library. It has an optional dependency on [_OpenMP_](https://www.openmp.org). The `Makefile` is designed to check for the presence of GCC, modules or a custom installation of _OpenMP_ using the definition of `$(OpenMP)` set to `/path/to/libomp`[^clang].

[^clang]: Tested on Apple's clang.

Compile everything[^compilation] with:

[^compilation]: The `-j` option is advised.

```bash
make
```

Compile examples with:

```bash
make examples
```

Compile tests with:

```bash
make tests
```

Compile a single test, which has to be set in the appropriate [Makefile](./Makefile) variable:

```bash
make single_test
```

Compile mesh generation scripts with:

```bash
make domains
```

Compile with debugging flag for `gdb`:

```bash
make debug
```

Executables are located in `executables/` and their outputs in `output/`.

#### Compilation Flags

:warning: Be sure to review the [Makefile](./Makefile) for additional options to customize the code's behavior.

The code uses the following custom compilation flags:

- `-DNVERBOSE`: Disables verbose output.
- `-DNSOLUTIONS`: Disables solution output.
- `-DNDEBUG`: Disables debug output.

### Compiling the Code into a Library

The code can also be compiled into the static library `PacsHPDG` under `./lib` by:

```bash
make lib
```

which can be installed locally under `~` by:

```bash
make install
```

## Using the Repository

Here is a brief guide on using the executables in this repository. Note that tests were written solely to verify specific parts of the code, so their usage is not covered here.

All output generated by the following codes is located in the `output/` directory.

### Domains

Codes written under `domains/` have the purpose of generating meshes over the domains used for testing the code, specifically the square and L-shaped domains.

This repository provides the following:

- `square_domain.cpp`
- `lshape_domain.cpp`
- `rectangular_domain.cpp`

For example, usage would look like:

```bash
./executables/square_domain.out 250
```

This command generates a square mesh over $[0, 1] \times [0, 1]$ with $N = 250$.

```bash
./executables/lshape_domain.out 125
```

This command generates an L-shaped mesh over $[-1, 1] \times [-1, 1] \setminus [0, 1] \times [-1, 0]$ with $N = 125$.

```bash
./executables/rectangular_domain.out 500
```

This command generates a square mesh over $[0, 5] \times [0, 1]$ with $N = 500$.


### Examples

Examples are divided into the following categories based on their equations and meshes:

1. Uniform meshes:
    - `laplace.cpp` Solve the Laplace equation.
    - `heat.cpp` Solve the Heat equation.
    - `fisher.cpp` Solve the Fisher equation.
2. _h-adaptively_ refined meshes with _a posteriori_ estimates:
    - `laplace_h.cpp` Solve the Laplace equation.
    - `heat_h.cpp` Solve the Heat equation.
    - `fisher_h.cpp` Solve the Fisher equation.
5. _hp-adaptively_ refined meshes with _a posteriori_ estimates:[^hp]
    - `laplace_hp.cpp` Solve the Laplace equation.
    - `heat_hp.cpp` Solve the Heat equation.
    - `fisher_hp.cpp` Solve the Fisher equation.

[^hp]: The polynomial degree specified for _hp-adaptively_ refined meshes is treated as the starting degree.

All the parameters related to the geometry of the problem are defined inside the problem data structures in `include/PacsHPDG/Data`; the `meshes/` folder contains various meshes already built for the different geometries, to use a different mesh one has just to change the mesh file in the desired examples scripts, where present. 

:warning: These examples contribute to the graphs presented in the report.

### Scripts

This repository includes some Python scripts to help visualize error trends, meshes and solutions.

Scripts are divided into the following categories based on their function:

1. Mesh-related:
    - `polyplot.py`: Requires a `.poly` file from which it plots a mesh. Accepts `--degrees` for _hp-adaptively_ refined meshes.
    - `estimatesplot.py`: Requires a `.poly` file from which it plots a mesh, with the relative error estimator distribution.
2. Solution-related:
    - `solplot.py`: Requires a `.sol` file from which it plots a solution.
3. Error trends:
    - `errorvsize.py`: Requires *one* `.error` file from which it plots error trends versus mesh size.
    - `errorvdofs.py`: Requires *one or two* `.error` files from which it plots error trends versus DOFs.
    - `hpvdofs.py`: Requires *one or two* `.error` files from which it plots error trends versus DOFs. Works for _hp-adaptively_ refined meshes.
4. TikZ wrappers:
    - `evs_tikz.py`: `errorvsize.py` TikZ wrapper.
    - `evd_tikz.py`: `errorvdofs.py` TikZ wrapper.
    - `hpvd_tikz.py`: `hpvdofs.py` TikZ wrapper.

## Using the Code

To use the code, include `<PacsHPDG.hpp>` which provides all necessary tools.

### Mesh Generation

Start by defining a domain using points:

```cpp
pacs::Point a{0.0, 0.0};
pacs::Point b{1.0, 0.0};
pacs::Point c{1.0, 1.0};
pacs::Point d{0.0, 1.0};

pacs::Polygon domain{{a, b, c, d}};
```

Next, create a diagram with 100 elements over the square domain:

```cpp
std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(domain, 100);
```

Generate a mesh and save it to a `.poly` file:

```cpp
pacs::Mesh mesh{domain, diagram};
mesh.write("output/square.poly");
```

The diagram can be retrieved using the same `mesh_diagram` function:

```cpp
std::vector<pacs::Polygon> retrieved = pacs::mesh_diagram("output/square.poly");
```

For non-convex domains, enable point reflection:

```cpp
std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(domain, 100, reflect=true);
```

### Solving the Poisson Problem

First, import the data related to the Laplace equation:
```cpp
DataLaplace data;
```

Then, construct the Laplacian objedct from your mesh[^laplacian][^real]:

[^laplacian]: `Laplace laplacian(mesh)` initialize all the members needed to solve the Laplace equation. 
[^real]: `pacs::Real` wraps `long double`.

The matrices are then built with the following method:

```cpp
laplacian.assembly(data, mesh);
laplacian.M();
laplacian.A();
laplacian.DG();
```

Next, construct the forcing term using specified source and Dirichlet boundary conditions, defined in the Lapalce equation data struct:

```cpp
pacs::Vector<pacs::Real> forcing = laplace.assembly_force(data, mesh);
```

Finally, solve the linear system to find the solution vector:

```cpp
pacs::Vector<pacs::Real> solution = laplacian.solve(mesh, forcing);
```

This `solution` vector now contains the computed solution to the Poisson problem on the given mesh with specified boundary conditions and source term.

### Mesh refinement

After solving the Poisson problem, you can evaluate _a posteriori_ error estimates, allowing for adaptive mesh refinement.

First, evaluate the estimates:

```cpp
pacs::LaplaceEstimator estimator(mesh);
estimator.computeEstimates(data, laplacian, numerical);
```

Use the `find_elem_to_refine` method to find which elements need to be refined, then use methods `mesh_refine_degree` and `mesh_refine_size` to refine the mesh:

```cpp
auto [h_mask, p_mask] = estimator.find_elem_to_refine();
estimator.mesh_refine_degree(p_mask);
estimator.mesh_refine_size(h_mask);
```

For Heat equation problem and Fisher-Kolmogorov equation problem the iter is pretty the same, with the add of other methods related to the fact that these are time dependant problems.

## Notes to the Reader

### On the Implementation of Basic Objects

This repository implements basic objects such as vectors and matrices, which are also found in many other libraries. Additionally, it includes a polygonal mesher, though it may not be as powerful as existing alternatives.

### On the Adaptation from **lymph**

The implementation of all the problems and FEM tools in this project was adapted from the [**lymph**](https://lymph.bitbucket.io) library.

### On Parallelism

Parallelism through the STL plays a secondary role, primarily utilized for common operations involving standard containers. In contrast, parallelism via _OpenMP_ is fundamental, significantly boosting the performance of the polygonal mesher, the `DB` solver, the `DBI` preconditioner and all the assembler, when looping over the mesh's elements.