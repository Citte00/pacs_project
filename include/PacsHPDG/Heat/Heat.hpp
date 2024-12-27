/**
 * @file Heat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief
 * @date 2024-12-16
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef INCLUDE_PACSHPDG_HEAT_HEAT_HPP
#define INCLUDE_PACSHPDG_HEAT_HEAT_HPP

#include "../Data.hpp"
#include "../Fem.hpp"

namespace pacs {

// Heat matrices.
std::array<Sparse<Real>, 3> heat(const DataHeat &, const Mesh &);

} // namespace pacs

#endif