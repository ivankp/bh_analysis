#include "BHEvent.hh"

#include <cmath>

#include <TTree.h>

inline void ActivateBranch(TTree* tree, const char* name, void* addr) noexcept {
  //TBranch *br = tree->GetBranch(name);
  tree->SetBranchAddress(name, addr);
  //br->SetStatus(true);
  tree->SetBranchStatus(name, true);
}

void BHEvent::SetTree(TTree* tree, select_t branches, bool old) {
  this->tree = tree;

  switch (branches) {
    case kinematics: {
    
      tree->SetBranchStatus("*",0);

      ActivateBranch(tree, "id", &eid);
      ActivateBranch(tree, "nparticle", &nparticle);
      ActivateBranch(tree, "px", px);
      ActivateBranch(tree, "py", py);
      ActivateBranch(tree, "pz", pz);
      ActivateBranch(tree, "E" , E );
      ActivateBranch(tree, "kf", kf);

      if (!old) { // for treating RS uncertainties properly
        ActivateBranch(tree, "part", part);
      }

    } break;
    case reweighting: {
    
      tree->SetBranchStatus("*",0);

      ActivateBranch(tree, "nparticle", &nparticle);
      ActivateBranch(tree, "px", &px);
      ActivateBranch(tree, "py", &py);
      ActivateBranch(tree, "kf", &kf);
      ActivateBranch(tree, "alphas", &alphas);
      ActivateBranch(tree, "weight2", &weight2);
      ActivateBranch(tree, "me_wgt", &me_wgt);
      ActivateBranch(tree, "me_wgt2", &me_wgt2);
      ActivateBranch(tree, "x1", &x[0]);
      ActivateBranch(tree, "x2", &x[1]);
      ActivateBranch(tree, "x1p", &xp[0]);
      ActivateBranch(tree, "x2p", &xp[1]);
      ActivateBranch(tree, "id1", &id[0]);
      ActivateBranch(tree, "id2", &id[1]);
      ActivateBranch(tree, "fac_scale", &fac_scale);
      ActivateBranch(tree, "ren_scale", &ren_scale);
      ActivateBranch(tree, "usr_wgts", usr_wgts);
      if (!old) {
        ActivateBranch(tree, "alphasPower", &alphas_power);
        ActivateBranch(tree, "part", part);
      }

    } break;
    case cross_section: {
    
      tree->SetBranchStatus("*",0);

      // ActivateBranch(tree, "nparticle", &nparticle);
      // ActivateBranch(tree, "kf", kf);
      ActivateBranch(tree, "weight", &weight);

    } break;
    default: {

      tree->SetBranchAddress("id", &eid);
      tree->SetBranchAddress("nparticle", &nparticle);
      tree->SetBranchAddress("px", px);
      tree->SetBranchAddress("py", py);
      tree->SetBranchAddress("pz", pz);
      tree->SetBranchAddress("E" , E);
      tree->SetBranchAddress("alphas", &alphas);
      tree->SetBranchAddress("kf", kf);
      tree->SetBranchAddress("weight", &weight);
      tree->SetBranchAddress("weight2", &weight2);
      tree->SetBranchAddress("me_wgt", &me_wgt);
      tree->SetBranchAddress("me_wgt2", &me_wgt2);
      tree->SetBranchAddress("x1", &x[0]);
      tree->SetBranchAddress("x2", &x[1]);
      tree->SetBranchAddress("x1p", &xp[0]);
      tree->SetBranchAddress("x2p", &xp[1]);
      tree->SetBranchAddress("id1", &id[0]);
      tree->SetBranchAddress("id2", &id[1]);
      tree->SetBranchAddress("fac_scale", &fac_scale);
      tree->SetBranchAddress("ren_scale", &ren_scale);
      tree->SetBranchAddress("nuwgt", &nuwgt);
      tree->SetBranchAddress("usr_wgts", usr_wgts);
      tree->SetBranchAddress("alphasPower", &alphas_power);
      tree->SetBranchAddress("part", part);

    }
  }

}

void BHEvent::SetPart(Char_t part) { this->part[0] = part; }
void BHEvent::SetAlphasPower(Char_t n) { this->alphas_power = n; }

template<typename T> constexpr T sq(T x) { return x*x; }

Double_t BHEvent::Ht() const noexcept {
  Double_t _Ht = 0.;
  for (Int_t i=0;i<nparticle;++i)
    _Ht += sqrt( sq(px[i]) + sq(py[i]) );
  return _Ht;
}

Double_t BHEvent::Ht_Higgs() const noexcept {
  Double_t _Ht = 0.;
  for (Int_t i=0;i<nparticle;++i) {
    Double_t pt2 = sq(px[i]) + sq(py[i]);
    if (kf[i]==25) pt2 += sq(125.); // mH^2
    _Ht += sqrt(pt2);
  }
  return _Ht;
}
