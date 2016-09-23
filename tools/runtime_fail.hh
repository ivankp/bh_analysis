#ifndef runtime_fail_hh
#define runtime_fail_hh

#include <exception>
#include <sstream>
#include <utility>

class runtime_fail : public std::exception {
  std::string str;

  template <typename Arg1>
  inline void impl(std::stringstream& ss, Arg1&& arg1) {
    ss << std::forward<Arg1>(arg1);
  }
  template <typename Arg1, typename... Args>
  inline void impl(std::stringstream& ss, Arg1&& arg1, Args&&... args) {
    ss << std::forward<Arg1>(arg1);
    impl(ss, std::forward<Args>(args)...);
  }

public:
  runtime_fail(runtime_fail&& o) : str(std::move(o.str)) { }
  runtime_fail(const runtime_fail& o) = delete;

  runtime_fail& operator= (runtime_fail&& o) {
    str = std::move(o.str);
    return *this;
  }
  runtime_fail& operator= (const runtime_fail& o) = delete;

  template <typename... Args>
  runtime_fail(Args&&... args) {
    std::stringstream ss;
    impl(ss, std::forward<Args>(args)...);
    str = ss.str();
  }
  virtual const char* what() const noexcept { return str.c_str(); }
};

#endif
