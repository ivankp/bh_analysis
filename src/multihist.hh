#ifndef MultiHist_HH
#define MultiHist_HH

#include <vector>

#include <TH1.h>

class TLegend;

class MultiHist {
  vector<TH1*> hists;
  vector<TColor> colors;

public:
  template <typename... TT>
  MultiHist(TT... args): hists{std::forward<TT...>(args)} { }

  decltype(hists)& hists() { return hists; }
  const decltype(hists)& hists() const { return hists; }

  template <typename... TT>
  void colors(TT... args) {
    colors = decltype(colors)(std::forward<TT...>(args));
  }
  decltype(colors)& colors() { return colors; }
  const decltype(colors)& colors() const { return colors; }

  TLegend* legend(TLegend* leg) const;

};

#endif
