/**
 * @file domain.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <iostream>
#include <vector>

int main() {

  using namespace pacs;

  // Constructs a domain.
  Point a{-1.0, -1.0};
  Point b{0.0, -1.0};
  Point c{0.0, 0.0};
  Point d{1.0, 0.0};
  Point e{1.0, 1.0};
  Point f{-1.0, 1.0};

  Polygon domain{{a, b, c, d, e, f}};

  // Constructing a mesh.
  Mesh mesh{domain, mesh_diagram(domain, 100, true)};

  // Mesh output.
  mesh.write("output/mesh.poly");
}