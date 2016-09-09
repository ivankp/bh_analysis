#ifndef IVANP_DELTA_PHI_HH
#define IVANP_DELTA_PHI_HH

// return absolute value of phi separation
template <typename T>
[[ gnu::const ]] inline T dphi(T phi1, T phi2) noexcept {
  static constexpr T twopi = M_PI*2.;

  T _dphi = phi1 - phi2;
  if (__builtin_expect(_dphi < 0.,0)) _dphi = -_dphi;
  return ( __builtin_expect(_dphi > M_PI,0) ? twopi-_dphi : _dphi );
}

#endif
