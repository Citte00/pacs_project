/**
 * @file line.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <iostream>

int main() {

  using namespace pacs;

  // Constructing two points.
  Point p{0.0, 1.0};
  Point q{1.0, 0.0};

  // Constructing some Lines.
  Line Bisector = bisector(p, q);
  Line bisector_symmetric = bisector(q, p);
  Line line{1.0, 1.0, 1.0};

  // Bisector output test.
  std::cout << Bisector << std::endl;
  std::cout << bisector_symmetric << std::endl;

  // Intersection output.
  std::cout << intersections(line, Bisector)[0] << std::endl;
    
}