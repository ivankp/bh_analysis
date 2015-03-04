#ifndef __hist_wt_hh__
#define __hist_wt_hh__

#include <unordered_map>
#include <string>
#include <memory>

#include <Rtypes.h>

#include "csshists.hh"

class TH1;
class TDirectory;

class weight;

class hist_wt {
  std::unordered_map<const weight*,TH1*> h;
public:
  hist_wt(const std::string& name);
  void Fill(Double_t x) noexcept;

  static std::shared_ptr<const csshists> css;
  static std::unordered_map<const weight*,TDirectory*> dirs;
};

#endif
