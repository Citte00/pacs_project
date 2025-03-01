/**
 * @file sparse.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <iostream>

using namespace pacs;

int main() {

    // Constructing a matrix.
    Sparse<Real> sparse{2, 3};
    
    // Insert.
    sparse.insert(0, 0, 1);
    sparse.insert(1, 1, -1);
    sparse.insert(0, 2, 1);

    // Compression.
    sparse.compress();

    // Output.
    std::cout << sparse << std::endl;

    // Vector product.
    Vector<Real> vector{3};
    
    vector[0] = 1;
    vector[1] = 2;
    vector[2] = 3;

    // Vector product output.
    std::cout << (sparse * vector) << std::endl;
    
}