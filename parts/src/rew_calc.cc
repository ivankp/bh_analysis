#include "rew_calc.hh"

#include <iostream>
#include <sstream>
#include <valarray>
#include <cmath>

#include <TTree.h>

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

#define debug

using namespace std;

constexpr array<int,10> quarks {
   1, 2, 3, 4, 5,
  -1,-2,-3,-4,-5,
};

template<typename T> inline T sq(T x) noexcept { return x*x; }

// PDF variables
string pdfname("def");
const LHAPDF::PDFSet* pdfset(nullptr);
vector<LHAPDF::PDF*> pdfs;
const LHAPDF::PDF* pdf; // central PDF
size_t npdfs;
double alphas_mH;

// PDF garbage collector
struct PDFgc {
  inline void clear() {
    delete pdfset;
    for ( auto pdf : pdfs ) delete pdf;
  }
  ~PDFgc() { clear(); }
} __pdf;

// Function to make PDFs
void usePDFset(const std::string& setname) {
  __pdf.clear();
  pdfset = new LHAPDF::PDFSet(setname);
  pdfname = pdfset->name();
  pdfs = pdfset->mkPDFs();
  pdf = pdfs[0];
  npdfs = pdfs.size();
  alphas_mH = pdf->alphasQ(125.);
}

//-----------------------------------------------
// Function classes to get scales values
//-----------------------------------------------

mu_fcn::mu_fcn(const string& str) noexcept : str(str) { }

mu_const::mu_const(double mu) noexcept
: mu_fcn( [&mu] () { // lambda
    stringstream ss;
    ss << mu << "GeV";
    return ss.str();
  } () ),
  _mu(mu)
{ }
double mu_const::mu() const noexcept {
  return _mu;
}

mu_fHt::mu_fHt(double fHt) noexcept
: mu_fcn( [&fHt] () { // lambda
    stringstream ss;
    ss << fHt << "Ht";
    return ss.str();
  } () ),
  fHt(fHt)
{ }
double mu_fHt::mu() const noexcept {
  return fHt*event.Ht();
}

mu_fHt_Higgs::mu_fHt_Higgs(double fHt) noexcept
: mu_fcn( [&fHt] () { // lambda
    stringstream ss;
    ss << fHt << "Ht";
    return ss.str();
  } () ),
  fHt(fHt)
{ }
double mu_fHt_Higgs::mu() const noexcept {
  return fHt*event.Ht_Higgs();
}

mu_fac_default::mu_fac_default() noexcept : mu_fcn("def") { }
double mu_fac_default::mu() const noexcept {
  return event.fac_scale;
}

mu_ren_default::mu_ren_default() noexcept : mu_fcn("def") { }
double mu_ren_default::mu() const noexcept {
  return event.ren_scale;
}

//-----------------------------------------------
// Reweighting computation
//-----------------------------------------------

std::vector<std::unique_ptr<const calc_base>> calc_base::all;

void calc_all_scales() noexcept {
  for (auto& c : calc_base::all ) c->calc();
}

// Factorization --------------------------------

const fac_calc*
mk_fac_calc(const mu_fcn* mu_f, bool unc, bool defaultPDF) noexcept {
  return new fac_calc(mu_f,unc,defaultPDF);
}

fac_calc::fac_calc(const mu_fcn* mu_f, bool pdf_unc, bool defaultPDF) noexcept
: calc_base(), mu_f(mu_f), pdf_unc(pdf_unc), defaultPDF(defaultPDF)
{
  all.push_back( unique_ptr<const calc_base>(this) );
}

fac_calc::~fac_calc() {
  delete mu_f;
}

// Renormalization ------------------------------

const ren_calc*
mk_ren_calc(const mu_fcn* mu_r, alphas_fcn asfcn, bool defaultPDF) noexcept {
  return new ren_calc(mu_r,asfcn,defaultPDF);
}

ren_calc::ren_calc(const mu_fcn* mu_r, alphas_fcn asfcn, bool defaultPDF) noexcept
: calc_base(), mu_r(mu_r), asfcn(asfcn), defaultPDF(defaultPDF), ar(1.)
{
  all.push_back( unique_ptr<const calc_base>(this) );
}

ren_calc::~ren_calc() {
  delete mu_r;
}

// Reweighter: combines fac and ren -------------

// Constructor creates branches on tree
reweighter::reweighter(const fac_calc* fac, const ren_calc* ren,
                       TTree* tree, bool pdf_unc)
: fac(fac), ren(ren), pdf_unc(pdf_unc), nk(pdf_unc ? 3 : 1)
{
  if (!pdfset) {
    cerr << "\033[31mNo PDF loaded\033[0m"  << endl;
    exit(1);
  }

  stringstream ss;
  ss <<  "Fac" << fac->mu_f->str
     << "_Ren" << ren->mu_r->str
     << "_PDF" << pdfname;

  string name( ss.str()+"_cent" );
  cout << "Creating branch: " << name << endl;
  tree->Branch(name.c_str(), &weight[0], (name+"/F").c_str());

  if (pdf_unc) {
    if (!fac->pdf_unc) {
      cerr << "PDF uncertanties are not set to be calculated for "
           << fac->mu_f->str << endl;
      exit(1);
    }

    name = ss.str()+"_down";
    cout << "Creating branch: " << name << endl;
    tree->Branch(name.c_str(), &weight[1], (name+"/F").c_str());

    name = ss.str()+"_up";
    cout << "Creating branch: " << name << endl;
    tree->Branch(name.c_str(), &weight[2], (name+"/F").c_str());
  }
}

reweighter::~reweighter() { }

/////////////////////////////////////////////////////////////////////
//
// Reweighting adopted from arXiv:1310.7439v1
// "Ntuples for NLO Events at Hadron Colliders" pp. 20 - 25
//
/////////////////////////////////////////////////////////////////////

//-----------------------------------------------
// Factorization
//-----------------------------------------------

inline double quark_sum(double x, double mu_fac) noexcept {
  double f = 0.;
  for (int q : quarks) f += pdf->xfxQ(q, x, mu_fac);
  return f;
}

// Get PDF lower and upper bound
valarray<double> xfxQ_unc(int id, double x, double q) noexcept {
  static vector<double> xfx(npdfs);
  for (size_t j=0;j<npdfs;++j) xfx[j] = pdfs[j]->xfxQ(id, x, q);

  static LHAPDF::PDFUncertainty xErr;
  pdfset->uncertainty(xErr, xfx); // no 3rd arg => use default 1 sigma

  return {
    xErr.central - xErr.errminus, // down
    xErr.central + xErr.errplus   // up
  };
}

valarray<double> quark_sum_unc(double x, double mu_fac) noexcept {
  valarray<double> f(0.,2);
  for (int q : quarks) f += xfxQ_unc(q, x, mu_fac);
  return f;
}

void unfold(const valarray<double>& a, double& x1, double& x2) noexcept {
  x1 = a[0];
  x2 = a[1];
}

void fac_calc::calc() const noexcept {
  // There is only one global event variable
  // so these references always points to the right place
  static const char& part = event.part[0];
  static Double_t* const& usr_wgts = event.usr_wgts;

  const double mu = mu_f->mu();

  // Born & Real
  if (defaultPDF) m[0] = event.weight;
  else {
    for (short i=0;i<2;++i)
      f[0][i][0] = pdf->xfxQ(event.id[i], event.x[i], mu)/event.x[i];

    // PDF uncertainty
    if (pdf_unc) for (short i=0;i<2;++i)
      unfold( xfxQ_unc(event.id[i], event.x[i], mu)/event.x[i],
              f[1][i][0], f[2][i][0] );

    m[0] = event.me_wgt2;
  }

  // Integrated subtraction
  if (part=='I') {

    // Calculate terms in Eq. (43)
    lf = 2.*log( mu / event.fac_scale );
    for (short i=1;i<9;++i)
      m[i] = usr_wgts[i+1] + usr_wgts[i+9]*lf;

    // Calculate terms in Eq. (44)
    static Int_t id;
    static Double_t x, xp;
    double si_[3][2] = { 0., 0., 0., 0., 0., 0. };
    for (short i=0;i<2;++i) {
      id = event.id[i];
      x  = event.x [i];
      xp = event.xp[i];

      f[0][i][1] = ( id==21 // Eq. (46)
                   ? quark_sum(    x, mu)/x
                   : pdf->xfxQ(id, x, mu)/x
      );
      f[0][i][2] = ( id==21 // Eq. (47)
                   ? quark_sum(    x/xp, mu)/x
                   : pdf->xfxQ(id, x/xp, mu)/x
      );
      f[0][i][3] = pdf->xfxQ(21, x,    mu)/x; // Eq. (48)
      f[0][i][4] = pdf->xfxQ(21, x/xp, mu)/x; // Eq. (49)

      if (pdf_unc) { // PDF uncertainty
        unfold( id==21 ? quark_sum_unc(x,    mu)/x : xfxQ_unc(id, x,    mu)/x,
                f[1][i][1], f[2][i][1]
        ); // Eq. (46)
        unfold( id==21 ? quark_sum_unc(x/xp, mu)/x : xfxQ_unc(id, x/xp, mu)/x,
                f[1][i][2], f[2][i][2]
        ); // Eq. (47)
        unfold( xfxQ_unc(21, x,    mu)/x, f[1][i][3], f[2][i][3] ); // Eq. (48)
        unfold( xfxQ_unc(21, x/xp, mu)/x, f[1][i][4], f[2][i][4] ); // Eq. (49)
      }
    }
    for (short k=0,nk=(pdf_unc?3:1);k<nk;++k) {
      for (short i=0;i<2;++i)
        for (short j=1;j<5;++j)
          si_[k][i] += f[k][i][j]*m[j+(i ? 4 : 0)]; // si_[k][0] = f_1^(j)*ω_j
                                                    // si_[k][1] = f_2^(j)*ω_j+4
      // f_1^(j)*ω_j*f_2 + f_1*f_2^(j)*ω_j+4
      si[k] = f[k][1][0]*si_[k][0] + f[k][0][0]*si_[k][1];
    }

  }
}

//-----------------------------------------------
// Renormalization
//-----------------------------------------------

void ren_calc::calc() const noexcept {
  // There is only one global event variable
  // so these references always points to the right place
  static const char& part = event.part[0];
  static const Double_t& alphas  = event.alphas;
  static const Char_t& n = event.alphas_power;
  static const Double_t& ren_scale = event.ren_scale;
  static Double_t* const& usr_wgts = event.usr_wgts;

  const double mu = mu_r->mu();

  // cout << endl;
  // test(mu_r->str)
  // test(ren_scale)
  // test(mu)

  // test(pdf->alphasQ(mu))
  // test(alphas)
  // test(int(n))

  // Calculate α_s change from renormalization
  if (!defaultPDF) {
    const double to_alphas = pdf->alphasQ(mu);
    switch (asfcn) {
      case alphas_fcn::all_mu:
        ar = pow( to_alphas/alphas, n ); break;
      case alphas_fcn::two_mH:
        ar = sq(to_alphas/alphas_mH)*pow(to_alphas/alphas, n-2); break;
    }
  }

  // test(ar)

  if (part=='V' || part=='I') {
    lr = 2.*log( mu / ren_scale ); // Calculate lr, same as l in Eq (30)
    m0 = lr*usr_wgts[0] + 0.5*lr*lr*usr_wgts[1];
  }

}

//-----------------------------------------------
// Stitch
//-----------------------------------------------

void reweighter::stitch() const noexcept {
  // There is only one global event variable
  // so these references always points to the right place
  const char part = event.part[0];

  for (short k=0;k<nk;++k) {
    s[k] = fac->m[0]; // m0

    if (part=='V' || part=='I') s[k] += ren->m0; // m0

    for (short i=0;i<2;++i) s[k] *= fac->f[k][i][0]; // m0 becomes s

    if (part=='I') s[k] += fac->si[k];

    weight[k] = s[k] * ren->ar;

    if (!isfinite(weight[k]))  {
      cerr << "\033[31mEvent " << event.eid << "\033[0m: "
           << "weight=" << weight[k] << endl;
      weight[k] = 0.;
    }
  }

  // test(weight[0])
}
