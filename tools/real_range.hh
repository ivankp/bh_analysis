#ifndef __real_range_hh__
#define __real_range_hh__

#include <istream>
#include <sstream>
#include <string>
#include <type_traits>
#include <stdexcept>

template<typename T = float>
struct real_range {
  static_assert(std::is_floating_point<T>::value,
                "real_range supports only floating point types");
  typedef typename std::remove_const<T>::type type;

  type min, max;

  real_range(): min(0), max(0) { }
  real_range(type min, type max): min(min), max(max) { }

  type len() const noexcept { return max - min; }
};

template<typename T>
std::istream& operator>> (std::istream& is, real_range<T>& r) {
  std::string str;
  is >> str;
  const size_t sep = str.find(':');
  if (sep==std::string::npos) {
    throw std::runtime_error(
      str+" cannot be converted to real_range. \':\' is missing."
    );
  } else {
    std::stringstream(str.substr(0,sep)) >> r.min;
    std::stringstream(str.substr(sep+1)) >> r.max;
  }
  return is;
}

#endif
