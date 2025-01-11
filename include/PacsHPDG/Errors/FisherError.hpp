/**
 * @file FisherError.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Fisher-KPP equation error computation class object.
 * @date 2025-01-11
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef INCLUDE_PACSHPDG_ERRORS_FISHER_ERROR_HPP
#define INCLUDE_PACSHPDG_ERRORS_FISHER_ERROR_HPP

#include "./HeatError.hpp"

namespace pacs {
class FisherError : public HeatError {
private:
  // Energy error member.
  Real m_energy;
  Real m_int_error;

public:
  // CONSTRUCTOR.
  FisherError(const Mesh &mesh_) : HeatError(mesh_) {
    this->m_energy = 0.0;
    this->m_int_error = 0.0;
  };

  // GETTERS.
  Real energy() const { return this->m_energy; };
  Real interror() const { return this->m_int_error; };

  // METHODS.

  // Compute errors.
  void computeErrors(const DataFKPP &, const Mesh &, const Fisher &);

  // Print to stream file for plotting.
  void print(std::ostream &os, const Mesh &mesh_) const override {
    os << "Elements: " << mesh_.elements.size() << "\n";
    os << "Dofs: " << mesh_.dofs() << "\n";
    os << "Degree (p): " << this->p() << "\n";
    os << "Size (h): " << this->h() << "\n";
    os << "L2 Error: " << this->L2error() << "\n";
    os << "DG Error: " << this->DGerror() << "\n";
    os << "Energy Error: " << this->energy() << std::endl;
  };
};

} // namespace pacs

#endif