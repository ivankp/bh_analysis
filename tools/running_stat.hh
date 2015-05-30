// taken from http://www.johndcook.com/blog/standard_deviation/

#ifndef running_stat_h
#define running_stat_h

#include <cmath>
#include <type_traits>

template<typename FT = double, typename IT = unsigned long,
  typename std::enable_if<std::is_floating_point<FT>::value>::type* = nullptr,
  typename std::enable_if<std::is_integral<IT>::value>::type* = nullptr
>
class running_stat {
public:
  using ft = FT;
  using it = IT;

private:
  it n;
  ft oldM, newM, oldS, newS;

public:
  running_stat(): n(0) { }

  inline void clear() noexcept { n = 0; }

  // ------------------------------------------------------

  inline void push(ft x) noexcept {
    ++n;

    // See Knuth TAOCP vol 2, 3rd edition, page 232
    if (n == 1) {
      oldM = newM = x;
      oldS = 0.;
    } else {
      newM = oldM + (x - oldM)/n;
      newS = oldS + (x - oldM)*(x - newM);

      // set up for next iteration
      oldM = newM; 
      oldS = newS;
    }
  }

  // ------------------------------------------------------

  inline it num  () const noexcept { return n; }
  inline ft mean () const noexcept { return (n > 0) ? newM : 0.; }
  inline ft var  () const noexcept { return ( (n > 1) ? newS/(n - 1) : 0. ); }
  inline ft stdev() const noexcept { return std::sqrt( var() ); }
};

#endif
