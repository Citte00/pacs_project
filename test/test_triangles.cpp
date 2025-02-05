/**
 * @file triangles.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <iostream>
#include <vector>

using namespace pacs;

int main() {

    // Constructs a domain.
    Point a{-1.5, -1.0};
    Point b{1.0, -2.0};
    Point c{2.0, 0.0};
    Point d{1.0, 1.0};
    Point e{-1.0, 1.0};

    Polygon polygon{{a, b, c, d, e}};

    // Centroid triangulation.
    std::vector<Polygon> triangles = triangulate(polygon);
}