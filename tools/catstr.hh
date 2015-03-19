#ifndef catstr_hh
#define catstr_hh

#include <string>
#include <sstream>

template<typename T>
inline void cat_impl(std::stringstream& ss, const T& t) {
  ss << t;
}

template<typename T, typename... TT>
inline void cat_impl(std::stringstream& ss, const T& t, const TT&... tt) {
  ss << t;
  cat_impl(ss,tt...);
}

template<typename... TT>
inline std::string cat(const TT&... tt) {
  std::stringstream ss;
  cat_impl(ss,tt...);
  return ss.str();
}

#endif
