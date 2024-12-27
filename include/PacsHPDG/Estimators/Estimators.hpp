/**
 * @file Estimators.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Error estimators polymorphic classes.
 * @date 2024-12-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef INCLUDE_PACSHPDG_ESTIMATORS_ESTIMATORS_HPP
#define INCLUDE_PACSHPDG_ESTIMATORS_ESTIMATORS_HPP

#include <array>
#include <iostream>

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Data.hpp"
#include "../Fem.hpp"
#include "../Geometry.hpp"

namespace pacs {

    /**
     * @brief Polymorphism base class error estimator.
     * 
     */
    class BaseEstimator {
    protected:
      // DOFs.
      std::size_t dofs;

      // Estimates.
      Real estimate = 0.0L;
      Vector<Real> estimates;

      // Fits.
      Vector<Real> fits;

    public:
      BaseEstimator(std::size_t num_elements = 0)
          : estimates{num_elements}, fits{num_elements} {};

      // Virtual constructor for polymorphic behaviour.
      virtual ~BaseEstimator() = default;

      // Virtual print method for polymorphic behavior.
      virtual void print(std::ostream &os) const = 0;

      // Virtual getters.
      virtual std::size_t getDofs() const = 0;
      virtual Real getEstimate() const = 0;
      virtual Vector<Real> getEstimates() const = 0;
      virtual Vector<Real> getFits() const = 0;

      // Virtual setters.
      virtual void setDofs(const size_t &dofs_) = 0;
      virtual void setEstimate(const Real &estimate_) = 0;
      virtual void setEstimates(const Vector<Real> &estimates_) = 0;
      virtual void setFits(const Vector<Real> &fits_) = 0;

      // Friend operator<< for polymorphic printing.
      friend std::ostream &operator<<(std::ostream &os,
                                      const BaseEstimator &estimator) {
        estimator.print(os);
        return os;
      }
    };

    /**
     * @brief Laplace equation error estimator class.
     *
     */
    class Estimator : public BaseEstimator {
    public:
      // CONSTRUCTOR.
      Estimator(const Mesh &, const Sparse<Real> &, const Vector<Real> &,
                const Functor &, const Functor &dirichlet = Functor{},
                const TwoFunctor &dirichlet_gradient = TwoFunctor{},
                const Real &penalty_coefficient = 10.0);

      // Getters.
      std::size_t getDofs() const override { return dofs; };
      Real getEstimate() const override { return estimate; };
      Vector<Real> getEstimates() const override { return estimates; };
      Vector<Real> getFits() const override { return fits; };

      // Setters.
      void setDofs(const size_t &dofs_) override { dofs = dofs; };
      void setEstimate(const Real &estimate_) override {
        estimate = estimate_;
      };
      void setEstimates(const Vector<Real> &estimates_) override {
        estimates = estimates_;
      };
      void setFits(const Vector<Real> &fits_) override { fits = fits_; };

      // Output.
      void print(std::ostream &) const override;
    };

    /**
     * @brief Heat equation error estimator class.
     * 
     */
    class HeatEstimator : public BaseEstimator {
        public:
          // CONSTRUCTOR.
          HeatEstimator(const DataHeat &, const Mesh &, const Sparse<Real> &,
                        const Vector<Real> &, const Vector<Real> &,
                        const Real &);

          // Getters.
          std::size_t getDofs() const override { return dofs; };
          Real getEstimate() const override { return estimate; };
          Vector<Real> getEstimates() const override { return estimates; };
          Vector<Real> getFits() const override { return fits; };

          // Setters.
          void setDofs(const size_t &dofs_) override { dofs = dofs; };
          void setEstimate(const Real &estimate_) override {
            estimate = estimate_;
          };
          void setEstimates(const Vector<Real> &estimates_) override {
            estimates = estimates_;
          };
          void setFits(const Vector<Real> &fits_) override { fits = fits_; };

          // Output.
          void print(std::ostream &) const override;
    };

}

#endif