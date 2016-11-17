// Written by Ivan Pogrebnyak

#ifndef IVANP_TIMED_COUNTER_CHRONO_HH
#define IVANP_TIMED_COUNTER_CHRONO_HH

#include <iostream>
#include <iomanip>
#include <chrono>

namespace ivanp {

template <typename I, typename Compare = std::less<I>>
class timed_counter {
public:
  using value_type    = typename std::enable_if<std::is_integral<I>::value,I>::type;
  using compare_type  = Compare;
  using clock_type    = std::chrono::system_clock;
  using time_type     = std::chrono::time_point<clock_type>;
  using duration_type = std::chrono::duration<double>;

private:
  value_type cnt;
  const value_type cnt_start, cnt_end;
  const time_type start;
  time_type last;
  compare_type cmp;

  void print() {
    using std::cout;
    using std::setw;
    using std::setfill;
    std::ios::fmtflags f(cout.flags());
    auto prec = cout.precision();

    const auto dt = duration_type((last = clock_type::now()) - start).count();
    const int hours   = dt/3600;
    const int minutes = (dt-hours*3600)/60;
    const int seconds = (dt-hours*3600-minutes*60);
    int nb = 28;

    cout << setw(12) << cnt << " | ";
    cout.precision(2);
    cout << std::fixed << setw(6) << (
      cnt==cnt_start ? 0. : 100.*float(cnt-cnt_start)/float(cnt_end-cnt_start)
    ) <<'%'<< " | ";
    if (hours) {
      nb += 8;
      cout << setw(5) << hours << ':'
      << setfill('0') << setw(2) << minutes << ':'
      << setw(2) << seconds << setfill(' ');
    } else if (minutes) {
      nb += 2;
      cout << setw(2) << minutes << ':'
      << setfill('0') << setw(2) << seconds << setfill(' ');
    } else {
      cout << setw(2) << seconds <<'s';
    }

    cout.flush();
    for (int i=0; i<nb; ++i) cout << '\b';
    cout.flags(f);
    cout.precision(prec);
  }
  void print_check() {
    if ( duration_type(clock_type::now()-last).count() > 1 ) print();
  }

public:
  timed_counter(I i, I n)
  : cnt(i), cnt_start(i), cnt_end(n), start(clock_type::now()), last(start)
  { print(); }
  timed_counter(I n)
  : cnt(0), cnt_start(0), cnt_end(n), start(clock_type::now()), last(start)
  { print(); }
  ~timed_counter() { print(); std::cout << std::endl; }

  inline bool ok() const noexcept { return cmp(cnt,cnt_end); }

  // prefix
  inline I operator++() { print_check(); return ++cnt; }
  inline I operator--() { print_check(); return --cnt; }

  // postfix
  inline I operator++(int) { print_check(); return cnt++; }
  inline I operator--(int) { print_check(); return cnt--; }

  template <typename T>
  inline I operator+= (T i) { print_check(); return cnt += i; }
  template <typename T>
  inline I operator-= (T i) { print_check(); return cnt -= i; }

  template <typename T>
  inline bool operator<  (T i) const noexcept { return cnt <  i; }
  template <typename T>
  inline bool operator<= (T i) const noexcept { return cnt <= i; }
  template <typename T>
  inline bool operator>  (T i) const noexcept { return cnt >  i; }
  template <typename T>
  inline bool operator>= (T i) const noexcept { return cnt >= i; }

  // cast to integral type
  template <typename T> inline operator T () {
    static_assert( std::is_integral<T>::value,
      "Cannot cast timed_counter to a non-integral type" );
    return cnt;
  }

  friend inline std::ostream&
  operator<<(std::ostream& os, const timed_counter& tc)
  noexcept(noexcept(os << tc.cnt)) {
    return (os << tc.cnt);
  }
};

} // end namespace

#endif
