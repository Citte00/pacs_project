/**
 * @file Functor.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Variadic template wrapper class for functions.
 * @date 2025-01-28
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef INCLUDE_PACSHPDG_FEM_FUNCTOR_HPP
#define INCLUDE_PACSHPDG_FEM_FUNCTOR_HPP

#include "../Algebra.hpp"
#include "../Base.hpp"

namespace pacs {

/**
 * @brief Variadic template Functor class.
 *
 */
template <typename ResultType, typename... Args> class Functor {
private:
  // Function.
  Function<ResultType, Args...> m_function;

public:
  // CONSTRUCTORS.
  Functor() : m_function{} {};
  explicit Functor(const Function<ResultType, Args...> &function_)
      : m_function{function_} {};
  template <typename Callable, typename = std::enable_if_t<!std::is_same_v<
                                   std::decay_t<Callable>, Functor>>>
  Functor(Callable &&callable) : m_function{std::forward<Callable>(callable)} {}

  // EVALUATION.
  ResultType operator()(const Args &...args) const {
#ifndef NDEBUG
    if (!this->m_function) {
      throw std::bad_function_call();
    }
#endif

    return this->m_function(args...);
  };
};

/**
 * @brief Variadic template TwoFunctor class.
 *
 */
template <typename ResultType, typename... Args> class TwoFunctor {
private:
  // Functions.
  Function<ResultType, Args...> m_first;
  Function<ResultType, Args...> m_second;

public:
  // CONSTRUCTORS.
  TwoFunctor() : m_first{}, m_second{} {};
  explicit TwoFunctor(const Function<ResultType, Args...> &first_,
                      const Function<ResultType, Args...> &second_)
      : m_first{first_}, m_second{second_} {};
  template <typename Callable, typename = std::enable_if_t<!std::is_same_v<
                                   std::decay_t<Callable>, TwoFunctor>>>
  TwoFunctor(Callable &&callable_f, Callable &&callable_s)
      : m_first{std::forward<Callable>(callable_f)},
        m_second{std::forward<Callable>(callable_s)} {}

  // EVALUATION.
  std::array<ResultType, 2> operator()(const Args &...args) const {
#ifndef NDEBUG
    if (!this->m_first || !this->m_second) {
      throw std::bad_function_call();
    }
#endif

    return {this->m_first(args...), this->m_second(args...)};
  };
};

// Some functor.
using BiFunctor = Functor<Vector<Real>, Vector<Real>, Vector<Real>>;
using TriFunctor = Functor<Vector<Real>, Vector<Real>, Vector<Real>, Real>;
using FKPPSource = Functor<Vector<Real>, Vector<Real>, Vector<Real>, Real,
                           Vector<Real>, Vector<Real>>;
using HeatSource =
    Functor<Vector<Real>, Vector<Real>, Vector<Real>, Real, Vector<Real>>;
using TriTwoFunctor =
    TwoFunctor<pacs::Vector<pacs::Real>, pacs::Vector<pacs::Real>,
               pacs::Vector<pacs::Real>, pacs::Real>;

} // namespace pacs

#endif