/**
 * @file HeatError.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Heat equation error computation class object.
 * @date 2025-01-11
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_ERRORS_HEAT_ERROR_HPP
#define INCLUDE_PACSHPDG_ERRORS_HEAT_ERROR_HPP

#include "./LaplaceError.hpp"

namespace pacs {

class HeatError : public LaplaceError {
    public:
    // CONSTRUCTOR.
    HeatError(const Mesh &mesh_) : LaplaceError(mesh_) {};

    // METHODS.
    // Compute L2 and H1 errors.
    void computeErrors(const DataHeat &, const Mesh &, const Heat &);
};

} // namespace pacs

#endif