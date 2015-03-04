#ifndef __int_range_hh__
#define __int_range_hh__

#include <istream>
#include <sstream>
#include <string>
#include <type_traits>

template<typename T = size_t>
struct int_range {
  static_assert(std::is_integral<T>::value,
                "int_range supports only integral types");
  typedef typename std::remove_const<T>::type type;

  type first, len;

  int_range(): first(0), len(0) { }
  int_range(type len): first(0), len(len) { }
  int_range(type first, type len): first(first), len(len) { }

  type end() const noexcept { return first + len; }
};

template<typename T>
std::istream& operator>> (std::istream& is, int_range<T>& r) {
  std::string str;
  is >> str;
  const size_t sep = str.find(':');
  if (sep==std::string::npos) {
    r.first = 0;
    std::stringstream(str) >> r.len;
  } else {
    std::stringstream(str.substr(0,sep)) >> r.first;
    std::stringstream(str.substr(sep+1)) >> r.len;
  }
  return is;
}

#endif
