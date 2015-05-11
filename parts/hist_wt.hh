#ifndef hist_wt_hh
#define hist_wt_hh

#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <memory>

#include <Rtypes.h>

#include "csshists.hh"

class TH1;
class TDirectory;

class weight;

class hist_wt {
  std::unordered_map<const weight*,
    std::tuple< TH1*, std::unordered_map<Int_t,Double_t>, Double_t* >> h;
  bool keep_sumw2;
public:
  hist_wt(const std::string& name);
  void Fill(Double_t x) noexcept;
  void FillIncl(Double_t x) noexcept;
  void FillSumw2() noexcept;
  void AdoptSumw2() noexcept;
  
  static std::shared_ptr<const csshists> css;
  static std::unordered_map<const weight*,TDirectory*> dirs;
  
  static std::vector<hist_wt*> all;
};

#endif
