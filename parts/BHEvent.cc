#include "BHEvent.hh"

#include <cmath>
#include <stdexcept>

#include <TTree.h>

inline void ActivateBranch(TTree* tree, const char* name, void* addr) noexcept {
  tree->SetBranchStatus (name, true);
  tree->SetBranchAddress(name, addr);
}

void BHEvent::SetTree(TTree* tree, select_t branches, bool old) {
  this->tree = tree;

  if (branches!=all) tree->SetBranchStatus("*",0);

  switch (branches) {
    case kinematics: {

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

      ActivateBranch(tree, "nparticle", &nparticle);
      ActivateBranch(tree, "px", px);
      ActivateBranch(tree, "py", py);
      ActivateBranch(tree, "kf", kf);
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

      ActivateBranch(tree, "id", &eid);
      ActivateBranch(tree, "weight2", &weight2);

    } break;
    case all: {

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
      if (!old) {
        tree->SetBranchAddress("alphasPower", &alphas_power);
        tree->SetBranchAddress("part", part);
      }

    }
  }

}

void BHEvent::SetPart(Char_t part) { this->part[0] = part; }
void BHEvent::SetAlphasPower(Char_t n) { this->alphas_power = n; }

template <typename T> inline T sq(T x) noexcept { return x*x; }
template <typename T, typename... TT>
inline T sq(T x, TT... xx) noexcept { return sq(x)+sq(xx...); }

template <typename T1, typename T2> inline bool in(T1 x, T2 a) { return x==a; }
template <typename T1, typename T2, typename... TT>
inline bool in(T1 x, T2 a, TT... bb) { return in(x,a) || in(x,bb...); }

Double_t BHEvent::Ht() const noexcept {
  Double_t _Ht = 0.;
  for (Int_t i=0; i<nparticle; ++i)
    _Ht += sqrt( sq(px[i],py[i]) );
  return _Ht;
}

Double_t BHEvent::Ht_Higgs() const noexcept {
  Double_t _Ht = 0.;
  for (Int_t i=0; i<nparticle; ++i) {
    Double_t pt2 = sq(px[i],py[i]);
    if (kf[i]==25) pt2 += sq(125.); // mH^2
    _Ht += sqrt(pt2);
  }
  return _Ht;
}

Double_t BHEvent::Ht_lnu() const {
  Double_t pt2 = 0.;
  bool found_l = false, found_nu = false;
  Double_t lnu[4] = {0.,0.,0.,0.};
  for (Int_t i=0; i<nparticle; ++i) {
    if (!found_l && in(abs(kf[i]),11,13,15)) {
      found_l = true;
      lnu[0] += px[i];
      lnu[1] += py[i];
      lnu[2] += pz[i];
      lnu[3] += E [i];
    } else if (!found_nu && in(abs(kf[i]),12,14,16)) {
      found_nu = true;
      lnu[0] += px[i];
      lnu[1] += py[i];
      lnu[2] += pz[i];
      lnu[3] += E [i];
    } else {
      pt2 += sq(px[i],py[i]);
    }
  }
  if (found_l && found_nu) pt2 += sq(lnu[0],lnu[1],lnu[2],lnu[3]);
  else throw std::runtime_error(
    "Cannot compute Ht_lnu for event without lepton and neutrino");
  return sqrt(pt2);
}

Double_t BHEvent::Et_gamma() const {
  // Return Et of the first photon we come across
  for (Int_t i=0; i<nparticle; ++i)
    if (kf[i]==22) return sq(px[i],py[i]);
  // No photon found
  throw std::runtime_error("Cannot compute Et_gamma for event without photon");
}
