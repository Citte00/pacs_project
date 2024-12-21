/**
 * @file Estimators.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Error estimators polymorphic classes.
 * @date 2024-12-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ESTIMATORS_PACS
#define ESTIMATORS_PACS

#include "../Base.hpp"
#include "../Fem.hpp"
#include "../Algebra.hpp"
#include "../Geometry.hpp"

#include <iostream>
#include <array>

namespace pacs {

    /**
     * @brief Polymorphism base class error estimator.
     * 
     */
    class BaseEstimator {
        public:
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
        virtual void setDofs(const size_t&) = 0;
        virtual void setEstimate(const Real&) = 0;
        virtual void setEstimates(const Vector<Real>&) = 0;
        virtual void setFits(const Vector<Real>&) = 0;

        // Friend operator<< for polymorphic printing.
        friend std::ostream &operator<<(std::ostream &os, const BaseEstimator &estimator) {
            estimator.print(os);
            return os;
        }
    };

    /**
     * @brief Laplace equation error estimator class.
     * 
     */
    class Estimator : public BaseEstimator {
        private:
            // DOFs.
            std::size_t dofs;

            // Estimates.
            Real estimate = 0.0L;
            Vector<Real> estimates;

            // Fits.
            Vector<Real> fits;

        public:
            // CONSTRUCTORS.
            Estimator(const Mesh &, const Sparse<Real> &, const Vector<Real> &, const Functor &, const Functor &dirichlet = Functor{}, const TwoFunctor &dirichlet_gradient = TwoFunctor{}, const Real &penalty_coefficient = 10.0);
            
            // Getter.
            std::size_t getDofs() const override { return this->dofs; };
            Real getEstimate() const override { return this->estimate; };
            Vector<Real> getEstimates() const override { return this->estimates; };
            Vector<Real> getFits() const override { return this->fits; };

            // Virtual setters.
            void setDofs(const size_t &dofs_) override { this->dofs = dofs; };
            void setEstimate(const Real &estimate_) override { this->estimate = estimate_; };
            void setEstimates(const Vector<Real> &estimates_) override { this->estimates = estimates_; };
            void setFits(const Vector<Real> &fits_) override { this->fits = fits_; };

            // Output.
            void print(std::ostream &) const override;

    };

    /**
     * @brief Heat equation error estimator class.
     * 
     */
    class HeatEstimator : public BaseEstimator {
        private:    
            // DOFs.
            std::size_t dofs;

            // Estimates.
            Real estimate = 0.0L;
            Vector<Real> estimates;

            // Fits.
            Vector<Real> fits;

        public:
            // CONSTRUCTORS.
            HeatEstimator(const Mesh &, const Sparse<Real> &, const Vector<Real> &, const Vector<Real> &, const HeatSource &, const Real &, const TriFunctor &, const TriFunctor &dirichlet = TriFunctor{}, const TriTwoFunctor &dirichlet_gradient = TriTwoFunctor{}, const Real &penalty_coefficient = 10.0);

            // Getter.
            std::size_t getDofs() const override { return this->dofs; };
            Real getEstimate() const override { return this->estimate; };
            Vector<Real> getEstimates() const override { return this->estimates; };
            Vector<Real> getFits() const override { return this->fits; };

            // Virtual setters.
            void setDofs(const size_t &dofs_) override { this->dofs = dofs; };
            void setEstimate(const Real &estimate_) override { this->estimate = estimate_; };
            void setEstimates(const Vector<Real> &estimates_) override { this->estimates = estimates_; };
            void setFits(const Vector<Real> &fits_) override { this->fits = fits_; };

            // Output.
            void print(std::ostream &) const override;

    };

}

#endif