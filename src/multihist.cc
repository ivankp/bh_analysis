#include "multihist.hh"

#include <TLegend.h>

TLegend* MultiHist::legend(TLegend* leg) const {
  for (TH1* h : hists) leg->AddEntry(h,h->GetTitle());
}

