#ifndef __function_traits_hh__
#define __function_traits_hh__

#include <tuple>

// Source: http://stackoverflow.com/a/7943765/2640636

template <typename T>
struct function_traits
  : public function_traits<decltype(&T::operator())>
{};
// For generic types, directly use the result of the signature of its 'operator()'

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...) const>
// we specialize for pointers to member function
{
  enum { num_args = sizeof...(Args) };
  // arity is the number of arguments.

  typedef ReturnType return_type;

  template <size_t i> struct arg
  {
    typedef typename std::tuple_element<i, std::tuple<Args...>>::type type;
    // the i-th argument is equivalent to the i-th tuple element of a tuple
    // composed of those arguments.
  };
};

#endif
