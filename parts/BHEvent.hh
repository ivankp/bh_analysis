#ifndef BHEvent_h
#define BHEvent_h

#define BHMAXNP 10 // maximum number of partons

#include <Rtypes.h>

class TTree;

class BHEvent {

  TTree *tree;
  void Init(TTree* tree);

public:
  Int_t           eid;
  Int_t           nparticle;
  Float_t         px[BHMAXNP]; //[nparticle]
  Float_t         py[BHMAXNP]; //[nparticle]
  Float_t         pz[BHMAXNP]; //[nparticle]
  Float_t         E [BHMAXNP]; //[nparticle]
  Int_t           kf[BHMAXNP]; //[nparticle]
  Double_t        alphas;
  Double_t        weight;
  Double_t        weight2;
  Double_t        me_wgt;
  Double_t        me_wgt2;
  Double_t        x[2];
  Double_t        xp[2];
  Int_t           id[2];
  Double_t        fac_scale;
  Double_t        ren_scale;
  Int_t           nuwgt;
  Double_t        usr_wgts[18]; //[nuwgt]
  Char_t          alphas_power;
  Char_t          part[2];

  Double_t Ht() const noexcept;
  Double_t Ht_Higgs() const noexcept;

  enum select_t { all, kinematics, reweighting, cross_section };

  void SetTree(TTree* tree, select_t branches=all, bool old=false);

  void SetPart(Char_t part);
  void SetAlphasPower(Char_t n);

};

#endif
