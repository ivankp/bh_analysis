#include "hist_wt.hh"

#include <TH1.h>
#include <TDirectory.h>

#include "weight.hh"

using namespace std;

template<typename T> inline T sq(T x) { return x*x; }

hist_wt::hist_wt(const string& name) {
  if (name.size()) {
    TH1* hist = css->mkhist(name);
    keep_sumw2 = ( hist->GetSumw2N() > 0 );
    for (auto& wt : weight::all) {
      const weight *w = wt.get();
      dirs[w]->cd();
      auto &hw = h[w];
      get<0>(hw) = static_cast<TH1*>( hist->Clone() );
      if (keep_sumw2) {
        const size_t n = hist->GetSumw2N();
        Double_t *w = get<2>(hw) = new Double_t[n];
        for (size_t i=0; i<n; ++i) w[i] = 0.;
      }
    }
    delete hist;

    all.push_back(this); // add to the static vector
  }
}

void hist_wt::Fill(Double_t x) noexcept {
  for (auto& _h : h) {
    auto *hist = get<0>(_h.second);
    const Double_t w = _h.first->get();
    const Int_t  bin = hist->FindFixBin(x);
    
    if (keep_sumw2)
      get<1>(_h.second).emplace(bin,0.).first->second += w;
    hist->Fill(x,w);
  }
}

void hist_wt::FillSumw2() noexcept {
  if (keep_sumw2) for (auto& _h : h) {
    for (auto it =get<1>(_h.second).begin();
              it!=get<1>(_h.second).end();
              it =get<1>(_h.second).erase(it) )
    {
      get<2>(_h.second)[it->first] += sq(it->second);
    }
  }
}

void hist_wt::AdoptSumw2() noexcept {
  if (keep_sumw2) for (auto& _h : h)
    get<0>(_h.second)->GetSumw2()->Adopt(
      get<0>(_h.second)->GetSumw2N(),
      get<2>(_h.second)
    );
}

shared_ptr<const csshists> hist_wt::css;
unordered_map<const weight*,TDirectory*> hist_wt::dirs;
vector<hist_wt*> hist_wt::all;
