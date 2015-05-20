#include "weight.hh"

using namespace std;

weight::weight(TTree *tree, const string& name, bool is_float)
: name(name), is_float(is_float)
{
  if ( tree->GetBranch(name.c_str()) ) {

    tree->SetBranchStatus(name.c_str(), true);

    if (is_float) tree->SetBranchAddress(name.c_str(), &w.f);
    else tree->SetBranchAddress(name.c_str(), &w.d);

  } else exit(1);
}

void weight::add(TTree* tree, const string& name, bool is_float) noexcept {
  all.emplace_back( new weight(tree,name,is_float) );
}

vector<unique_ptr<const weight>> weight::all;
