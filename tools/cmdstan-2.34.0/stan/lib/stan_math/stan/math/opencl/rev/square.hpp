#ifndef STAN_MATH_OPENCL_REV_SQUARE_HPP
#define STAN_MATH_OPENCL_REV_SQUARE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `square()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `square()` of the input.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> square(const var_value<T>& A) {
  return make_callback_var(square(A.val()),
                           [A](vari_value<matrix_cl<double>>& res) mutable {
                             A.adj() += 2.0 * elt_multiply(res.adj(), A.val());
                           });
}

}  // namespace math
}  // namespace stan

#endif
#endif