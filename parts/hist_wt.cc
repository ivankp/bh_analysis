#include "hist_wt.hh"

#include <TH1.h>
#include <TDirectory.h>

#include "weight.hh"

using namespace std;

template<typename T> inline T sq(T x) { return x*x; }

hist_wt::hist_wt(hist_wt&& o): h(std::move(o.h)), keep_sumw2(o.keep_sumw2) { }
hist_wt& hist_wt::operator=(hist_wt&& o) {
  h = std::move(o.h);
  keep_sumw2 = o.keep_sumw2;
  return *this;
}

hist_wt::hist_wt(const string& name) {
  if (name.size()) init(name);
}

void hist_wt::init(const string& name) {
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

void hist_wt::Fill(Double_t x) noexcept {
  // all histograms in a set have the same binning
  const Int_t bin = get<0>(h.cbegin()->second)->FindFixBin(x);

  for (auto& _h : h) {
    auto *hist = get<0>(_h.second);
    const Double_t w = _h.first->get();

    if (keep_sumw2)
      get<1>(_h.second).emplace(bin,0.).first->second += w;

    hist->Fill(x,w/**event_ptr->ncount*/);
  }
}

void hist_wt::FillIncl(Double_t x) noexcept {
  const TH1* hist = get<0>(h.cbegin()->second);

  for (Int_t i=1, n=hist->FindFixBin(x); i<=n; ++i)
    Fill(hist->GetBinCenter(i));
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
