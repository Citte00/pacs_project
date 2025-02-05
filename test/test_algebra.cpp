/**
 * @file algebra.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <iostream>

int main() {

    // Sparse.
    using namespace pacs;

    // Constructs a Matrix.
    Sparse<Real> sparse{2, 2};

    sparse.insert(0, 0, 4.0);
    sparse.insert(0, 1, 1.0);
    sparse.insert(1, 0, 1.0);
    sparse.insert(1, 1, 3.0);

    sparse.compress();

    // Constructs a Vector.
    Vector<Real> vector{2};

    vector[0] = 1.0;
    vector[1] = 2.0;

    // Linear system (Ax = b) solution.
    std::cout << solve(sparse, vector, GMRES) << std::endl; // GMRES.

    // Dense.
    
    // Constructs a Matrix.
    pacs::Matrix<pacs::Real> dense{2, 2};

    dense(0, 0) = 4.0;
    dense(0, 1) = 1.0;
    dense(1, 0) = 1.0;
    dense(1, 1) = 3.0;

    // Linear system (Ax = b) solution.
    std::cout << solve(dense, vector, QRD) << std::endl; // QR.
    std::cout << solve(dense, vector, LUD) << std::endl; // LU.
    
}