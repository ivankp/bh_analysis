#include "hist_wt.hh"

#include <TH1.h>
#include <TDirectory.h>

#include "weight.hh"

using namespace std;

hist_wt::hist_wt(const string& name, Char_t* part_ptr): part(part_ptr) {
  if (name.size()) {
    TH1* hist = css->mkhist(name);
    for (auto& wt : weight::all) {
      const weight *w = wt.get();
      dirs[w]->cd();
      h[w] = make_pair( static_cast<TH1*>( hist->Clone() ), 0. );
    }
    delete hist;

    all.push_back(this); // add to the static vector
  }
}

void hist_wt::Fill(Double_t x) noexcept {
  if (*part=='R') {
    for (auto& _h : h)
      _h.second.second += _h.first->get();
  } else {
    for (auto& _h : h)
      _h.second.first->Fill(x,_h.first->get());
  }
}

void hist_wt::FillIndirect() noexcept {
  for (auto& _h : h) {
    _h.second.first->Fill(_h.second.second);
    _h.second.second = 0.;
  }
}

shared_ptr<const csshists> hist_wt::css;
unordered_map<const weight*,TDirectory*> hist_wt::dirs;
vector<hist_wt*> hist_wt::all;
