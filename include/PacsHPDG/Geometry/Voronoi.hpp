/**
 * @file Voronoi.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef INCLUDE_PACSHPDG_GEOMETRY_VORONOI_HPP
#define INCLUDE_PACSHPDG_GEOMETRY_VORONOI_HPP

#include "../Base.hpp"
#include "./Shapes.hpp"

#include <vector>

#ifndef GEOMETRY_SAFE
#define GEOMETRY_SAFE 1E-6
#endif

#ifndef GEOMETRY_PADDING
#define GEOMETRY_PADDING 5E-2 // better use constexpr than cpp macros
#endif

namespace pacs {

    // VORONOI.

    std::vector<Polygon> voronoi(const Polygon &, const std::vector<Point> &, const bool &reflect = false);
    std::vector<Polygon> voronoi_random(const Polygon &, const std::size_t &, const bool &reflect = false);
    std::vector<Polygon> voronoi_uniform(const Polygon &, const std::size_t &, const bool &reflect = false);

    // TRIANGLES.

    std::vector<Polygon> triangulate(const Polygon &);
    std::vector<Polygon> triangulate(const std::vector<Polygon> &);

}

#endif