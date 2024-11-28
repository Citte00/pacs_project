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

#include "../Base.hpp"
#include "../Algebra.hpp"

namespace pacs {

    /**
     * @brief General Functor class.
     * 
     */
    template<typename VariableType, typename... Args>
    class GeneralFunctor {
        private:
            
            // Function.
            GenFunc<VariableType, Args...> m_function;

        public:

            // CONSTRUCTORS.
            
            GeneralFunctor() = default;
            GeneralFunctor(const GenFunc<VariableType, Args...>& function_) : m_function{function_} {};

            // EVALUATION.
            auto operator()(const Args&... args) const {
                // Determine the return type based on inputs
                if constexpr ((IsVector<Args> || ...)) {
                    // At least one argument is a vector: perform element-wise operation
                    std::size_t length = 0;
                    ((length = IsVector<Args> ? args.length : length), ...);

                    Vector<Real> result(length);

                    // Compute element-wise results
                    for (std::size_t i = 0; i < length; ++i) {
                        result[i] = m_function((IsVector<Args> ? args[i] : args)...);
                    }

                    return result;
                    
                } else {
                    // All arguments are scalars: directly apply the function
                    return m_function(args...);
                }
            };

    };

}

#endif