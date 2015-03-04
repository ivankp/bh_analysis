#ifndef __csshists_hh__
#define __csshists_hh__

#include <string>

class TH1;

class csshists {
  class impl;
  impl *_impl;

public:
  csshists(const std::string& cssfilename);
  ~csshists();

  TH1* mkhist(const std::string& name) const;
};

#endif
