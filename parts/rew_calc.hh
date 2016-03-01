#ifndef rew_calc_h
#define rew_calc_h

#include <string>
#include <vector>
#include <algorithm>

#include "BHEvent.hh"

extern BHEvent event;

// Function to make PDFs
void usePDFset(const std::string& setname);

//-----------------------------------------------
// Function classes to get scales values
//-----------------------------------------------

class mu_fcn {
public:
  virtual double mu() const noexcept =0;
  virtual ~mu_fcn() noexcept { }
};

class mu_fixed: public mu_fcn {
  double _mu;
public:
  mu_fixed(double mu) noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fixed() noexcept { }
};

class mu_fHt: public mu_fcn {
  double fHt;
public:
  mu_fHt(double fHt) noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fHt() noexcept { }
};

class mu_fHt_Higgs: public mu_fcn {
  double fHt;
public:
  mu_fHt_Higgs(double fHt) noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fHt_Higgs() noexcept { }
};

class mu_fEt_gamma: public mu_fcn {
  double fEt;
public:
  mu_fEt_gamma(double fEt) noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fEt_gamma() noexcept { }
};

class mu_fac_default: public mu_fcn {
public:
  mu_fac_default() noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_fac_default() noexcept { }
};

class mu_ren_default: public mu_fcn {
public:
  mu_ren_default() noexcept;
  virtual double mu() const noexcept;
  virtual ~mu_ren_default() noexcept { }
};

//-----------------------------------------------
// Factorization --------------------------------
//-----------------------------------------------

class reweighter;

struct fac_calc {
  fac_calc(const mu_fcn* mu_f) noexcept;
  void calc() const noexcept;
  ~fac_calc();

  bool pdf_unc;
  bool defaultPDF;

private:
  const mu_fcn* mu_f;
  mutable double f[3][2][5], m[9], lf, si[3];

friend class reweighter;
};

//-----------------------------------------------
// Renormalization ------------------------------
//-----------------------------------------------

enum class alphas_fcn: char { all_ren, two_mH };

struct ren_calc {
  ren_calc(const mu_fcn* mu_r) noexcept;
  void calc() const noexcept;
  ~ren_calc();

  alphas_fcn new_alphas;
  bool defaultPDF;

  static alphas_fcn bh_alphas;

private:
  const mu_fcn* mu_r;
  mutable double ar, lr, m0;

friend class reweighter;
};

//-----------------------------------------------
// Reweighter: combines fac and ren -------------
//-----------------------------------------------

class reweighter {
  const fac_calc *fac;
  const ren_calc *ren;
  bool pdf_unc;
  short nk;

  mutable double s[3];

  mutable Float_t weight[3];

public:
  // Constructor creates branches on tree
  reweighter(const std::pair<const fac_calc*,std::string>& fac,
             const std::pair<const ren_calc*,std::string>& ren,
             TTree* tree, bool pdf_unc=false);
  ~reweighter();
  void stitch() const noexcept;
};

#endif
