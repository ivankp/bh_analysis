#include "hist_wt.hh"

#include <TH1.h>
#include <TDirectory.h>

#include "weight.hh"

using namespace std;

hist_wt::hist_wt(const string& name) {
  TH1* hist = css->mkhist(name);
  for (auto& wt : weight::all) {
    const weight *w = wt.get();
    dirs[w]->cd();
    h[w] = static_cast<TH1*>( hist->Clone() );
  }
  delete hist;
}

void hist_wt::Fill(Double_t x) noexcept {
  for (auto& _h : h)
    _h.second->Fill(x,_h.first->is_float ? _h.first->w.f : _h.first->w.d);
}

shared_ptr<const csshists> hist_wt::css;
unordered_map<const weight*,TDirectory*> hist_wt::dirs;
