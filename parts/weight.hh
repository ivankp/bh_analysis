#ifndef weight_hh
#define weight_hh

#include <string>
#include <vector>
#include <memory>

#include <TTree.h>

// Weights collector ************************************************

class weight {
public:
  std::string name;

private:
  union {
    Double_t d;
    Float_t f;
  } w;
  bool is_float;

public:
  weight(TTree *tree, const std::string& name, bool is_float=true);

  inline Double_t get() const noexcept {
    return (is_float ? w.f : w.d);
  }

  static std::vector<std::unique_ptr<const weight>> all;
  static void add(TTree* tree, const std::string& name, bool is_float=true) noexcept;
};

#endif
