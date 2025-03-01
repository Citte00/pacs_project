/**
 * @file Mesh.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief
 * @date 2024-05-04
 *
 * @copyright Copyright (c) 2024
 *
 */
#ifndef INCLUDE_PACSHPDG_GEOMETRY_MESH_HPP
#define INCLUDE_PACSHPDG_GEOMETRY_MESH_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"

#include "./Shapes.hpp"

#include <array>
#include <string>
#include <vector>

namespace pacs {

// Element, implemented under Element.cpp

/**
 * @brief Element struct. Polygon wrapper.
 *
 */
struct Element {

  // Polygon.
  Polygon element;

  std::vector<Point> nodes;
  std::vector<Segment> edges;

  // Polynomial degree.
  std::size_t degree;

  // CONSTRUCTORS.

  Element(const Polygon &);
  Element(const Polygon &, const std::size_t &);

  // METHODS.

  inline std::size_t dofs() const {
    return (this->degree + 1) * (this->degree + 2) / 2;
  }
};

// Mesh, implemented under Mesh.cpp

/**
 * @brief Mesh struct.
 *
 */
struct Mesh {

  // Geometry.
  Polygon domain; // Polygonal domain.
  std::vector<Point> coord;
  std::vector<std::vector<int>> connectivity;

  // Elements.
  std::vector<Element> elements;

  // Neighbours.
  std::vector<std::vector<std::array<int, 3>>> neighbours;

  // Areas.
  std::vector<Real> areas;
  std::vector<Vector<Real>> max_simplices;

  // Entries for the solution.
  std::size_t entries;

  // CONSTRUCTORS.

  Mesh(const Polygon &, const std::vector<Polygon> &,
       const std::vector<std::size_t> &);
  Mesh(const Polygon &, const std::vector<Polygon> &,
       const std::size_t &degree = 1);
  Mesh(const Mesh &);

  Mesh &operator=(const Mesh &);

  // READ, WRAPPERS.

  Polygon element(const std::size_t &) const;

  // STATS.

  std::size_t dofs() const;

  // OUTPUT.

  void write(const std::string &, const bool &degrees = false);
};

// METHODS.
// Implemented under src/Builder.cpp

// Diagrams.
std::vector<Polygon> mesh_diagram(const Polygon &, const std::size_t &,
                                  const bool &reflect = false,
                                  const bool &uniform = false);
std::vector<Polygon> mesh_diagram(const std::string &);
std::vector<Polygon> mesh_relax(const Polygon &, const std::vector<Polygon> &,
                                const bool &reflect = false);

// Data.
std::vector<Element> mesh_elements(const std::vector<Polygon> &,
                                   const std::vector<std::size_t> &);
std::vector<std::vector<std::array<int, 3>>>
mesh_neighbours(const Polygon &, const std::vector<Element> &);

std::vector<Real> mesh_areas(const std::vector<Polygon> &);
std::vector<Vector<Real>> mesh_max_simplices(const std::vector<Polygon> &);

} // namespace pacs

#endif