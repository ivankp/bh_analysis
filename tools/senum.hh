#ifndef senum_hh
#define senum_hh

#include <iostream>
#include <string>
#include <stdexcept>

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#define ENUM_STR_EN_CASE(r, data, elem) \
  case elem: return BOOST_PP_STRINGIZE(elem);

#define ENUM_STR_STR_CASE(r, enum_name, elem) \
  if (!str.compare( BOOST_PP_STRINGIZE(elem) )) return elem; else

#define senum(enum_name, values) \
  struct enum_name { \
    enum type { BOOST_PP_SEQ_ENUM(values) }; \
    static type get_enum(const std::string& str) { \
      BOOST_PP_SEQ_FOR_EACH( ENUM_STR_STR_CASE, nil, values ) \
      throw std::runtime_error( \
        "invalid value \""+str+"\" specified for stringized enum " \
        BOOST_PP_STRINGIZE(enum_name) \
      ); \
    } \
    static const char* get_str(type en) { \
      switch (en) { \
        BOOST_PP_SEQ_FOR_EACH( ENUM_STR_EN_CASE, nil, values ) \
      } \
      return 0; \
    } \
  }; \
\
std::istream& operator>> (std::istream& is, enum_name::type& x) { \
  std::string str; \
  is >> str; \
  x = enum_name::get_enum(str); \
  return is; \
} \
\
std::ostream& operator<< (std::ostream& os, const enum_name::type& x) { \
  os << enum_name::get_str(x); \
  return os; \
}
  
#endif
