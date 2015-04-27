#ifndef hist_wt_hh
#define hist_wt_hh

#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

#include <Rtypes.h>

#include "csshists.hh"

class TH1;
class TDirectory;

class weight;

class hist_wt {
  Char_t *part;
  std::unordered_map<const weight*,std::pair<TH1*,Double_t>> h;
public:
  hist_wt(const std::string& name, Char_t* part_ptr);
  void Fill(Double_t x) noexcept;
  void FillIndirect() noexcept;

  static std::shared_ptr<const csshists> css;
  static std::unordered_map<const weight*,TDirectory*> dirs;
  
  static std::vector<hist_wt*> all;
};

#endif
