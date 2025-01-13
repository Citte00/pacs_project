/**
 * @file GeneralFunctor.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief template Functor class for a more general use
 * @date 2024-11-27
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef GENERAL_FUNCTOR_PACS
#define GENERAL_FUNCTOR_PACS

#include "../Algebra.hpp"
#include "../Base.hpp"
#include "../Geometry.hpp"

namespace pacs {

/**
 * @brief Function template alias.
 *
 * @tparam ResultType
 * @tparam Params
 */
template <typename ResultType, typename... Params>
using GenFunc = std::function<ResultType(Params...)>;

/**
 * @brief General Functor class.
 *
 */
template <typename ResultType, typename... Args> 
class GeneralFunctor {
private:
  // Function.
  GenFunc<ResultType, Args...> m_function;

public:
  // CONSTRUCTORS.
  GeneralFunctor() : m_function{} {};
  explicit GeneralFunctor(const GenFunc<ResultType, Args...> &function_)
      : m_function{function_} {};
  template <typename Callable, typename = std::enable_if_t<!std::is_same_v<
                                   std::decay_t<Callable>, GeneralFunctor>>>
  GeneralFunctor(Callable &&callable)
      : m_function{std::forward<Callable>(callable)} {}

  // Conversion operator.
      operator GenFunc<ResultType, Args...>() const {
        return m_function;
    };

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
 * @brief General TwoFunctor class.
 *
 */
template <typename ResultType, typename... Args> class GeneralTwoFunctor {
private:
  // Functions.
  GenFunc<ResultType, Args...> m_first;
  GenFunc<ResultType, Args...> m_second;

public:
  // CONSTRUCTORS.
  GeneralTwoFunctor() : m_first{}, m_second{} {};
  explicit GeneralTwoFunctor(const GenFunc<ResultType, Args...> &first_,
                             const GenFunc<ResultType, Args...> &second_)
      : m_first{first_}, m_second{second_} {};
  template <typename Callable, typename = std::enable_if_t<!std::is_same_v<
                                   std::decay_t<Callable>, GeneralTwoFunctor>>>
  GeneralTwoFunctor(Callable &&callable_f, Callable &&callable_s)
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
using BiFunctor = GeneralFunctor<Vector<Real>, Vector<Real>, Vector<Real>>;
using TriFunctor =
    GeneralFunctor<Vector<Real>, Vector<Real>, Vector<Real>, Real>;
using FKPPSource = GeneralFunctor<Vector<Real>, Vector<Real>, Vector<Real>,
                                  Real, Vector<Real>, Vector<Real>>;
using HeatSource = GeneralFunctor<Vector<Real>, Vector<Real>, Vector<Real>,
                                  Real, Vector<Real>>;
using BiTwoFunctor =
    GeneralTwoFunctor<Vector<Real>, Vector<Real>, Vector<Real>>;
using TriTwoFunctor =
    GeneralTwoFunctor<Vector<Real>, Vector<Real>, Vector<Real>, Real>;

} // namespace pacs

#endif