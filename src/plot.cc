#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include <memory>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TLatex.h>
#include <TLine.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>

using namespace std;
namespace po = boost::program_options;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

template<class T>
inline T* get(TDirectory* dir, const char* name) {
  T *obj = nullptr;
  dir->GetObject(name,obj);
  if (obj) return obj;
  else {
    cerr << "No object " << name << " in " << dir->GetName() << endl;
    exit(1);
  }
}

template<class T>
void get_all(TDirectory* dir, vector<T*>& hists) noexcept {
  TIter nextkey(dir->GetListOfKeys());
  TKey *key;
  while ((key = static_cast<TKey*>(nextkey()))) {
    TObject *obj = key->ReadObj();
    if (obj->InheritsFrom(T::Class()))
      hists.push_back( static_cast<T*>(obj) );
  }
}

string sigma_prt(Double_t sigma, unsigned prec) {
  stringstream ss;
  ss << "#sigma = " << showpoint << setprecision(prec) << sigma
     << noshowpoint << " pb";
  return ss.str();
}

enum class var : char { other, N, pT, mass, y, phi, HT, tau };

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

template<class T>
string combine(const T& container) {
  string s;
  auto it  = container.begin();
  auto end = container.end();
  if (it==end) return s;
  s += *(it++);
  for (;it!=end;++it) s += ", " + *it;
  return s;
}

int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  string fin_name, title;
  string fout_name;
  float of_lim;
  bool of_lim_set = false;
  unsigned sigma_prec;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<string>(&fin_name)->required(),
     "*input files with histograms")
    ("output,o", po::value<string>(&fout_name),
     "output pdf plots")
    ("title,t", po::value<string>(&title),
     "string appended to each title")
    ("sigma-prec", po::value<unsigned>(&sigma_prec)->default_value(3),
     "number of significant digits in cross section")
    ("overflow", po::value<float>(&of_lim),
     "print under- and overflow messages on histograms above the limit")
    ;

    po::positional_options_description pos;
    pos.add("input",1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(all_opt).positional(pos).run(), vm);
    if (argc == 1 || vm.count("help")) {
      cout << all_opt << endl;
      return 0;
    }
    po::notify(vm);

    if (!fout_name.size()) {
      size_t first = fin_name.rfind('/')+1;
      if (first==string::npos) first = 0;
      fout_name = fin_name.substr(first,fin_name.rfind('.')-first)+".pdf";
    }
    if (vm.count("overflow")) of_lim_set = true;

  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Open files
  TFile *fin = new TFile(fin_name.c_str(),"read");
  if (fin->IsZombie()) exit(1);
  cout << "Input file: " << fin->GetName() << endl;
  cout << "Output file: " << fout_name << endl;
  cout << endl;

  // Number of events
  TH1D *h_N = get<TH1D>(fin,"N");
  const Double_t N_events = h_N->GetAt(1);
  cout << "Events: " << lround(h_N->GetEntries()) << endl;

  // Properties
  set<string> pdf_name, jet_alg;

  // Collect histograms
  vector<TDirectory*> dirs;
  dirs.reserve(16);
  get_all(fin,dirs);
  for (auto d : dirs) {
    const string dir_name(d->GetName());
    cout << dir_name << endl;
    static boost::char_separator<char> sep("_");
    tokenizer tok(dir_name,sep);
    for (auto t : tok) {
      //cout << t << endl;
      if (!t.substr(0,3).compare("PDF")) pdf_name.emplace(t,3);
      if (!t.substr(0,3).compare("Jet")) jet_alg .emplace(t,3);
    }
  }
  cout << endl;
  const size_t ndirs = dirs.size();

  vector<TH1D*> central;
  central.reserve(32);
  get_all(dirs[0],central);

  // Draw histograms
  TCanvas canv;
  canv.SaveAs((fout_name+'[').c_str());

  for (auto h_cent : central) {
    const string h_name(h_cent->GetName());
    cout << h_name << endl;
    const size_t nbins = h_cent->GetNbinsX();

    // Variable type
    var h_var = var::other;

    if (h_name.find("_N") != string::npos) {
      h_var = var::N;
    } else if (h_name.find("_pT") != string::npos) {
      h_var = var::pT;
    } else if (h_name.find("_mass") != string::npos) {
      h_var = var::mass;
    } else if ( (h_name.find("_y") != string::npos) ||
                (h_name.find("_dy") != string::npos) ||
                (h_name.find("_deltay") != string::npos) ) {
      h_var = var::y;
    } else if ( (h_name.find("_phi") != string::npos) ||
                (h_name.find("_dphi") != string::npos) ||
                (h_name.find("_deltaphi") != string::npos) ) {
      h_var = var::phi;
    } else if (h_name.find("_HT") != string::npos) {
      h_var = var::HT;
    } else if (h_name.find("_tau") != string::npos) {
      h_var = var::tau;
    }

    // Cross section
    const Double_t sigma   = ( h_var == var::N && h_name.find("incl") != string::npos
                           ? h_cent->GetAt(1)
                           : h_cent->Integral(0,nbins+1) )/N_events;
    const Double_t sigma_u = h_cent->GetAt(0)/N_events;
    const Double_t sigma_o = h_cent->GetAt(nbins+1)/N_events;

    // Book vectors
    vector<Double_t> bins_edge(nbins,0.),
                     bins_wdth(nbins,0.),
                     cent     (nbins,0.),
                     scales_lo(nbins,0.),
                     scales_hi(nbins,0.),
                     pdf_lo   (nbins,0.),
                     pdf_hi   (nbins,0.);

    // Fill central vectors
    for (size_t i=0;i<nbins;++i) {
      bins_edge[i] = h_cent->GetBinLowEdge(i+1);
      bins_wdth[i] = h_cent->GetBinLowEdge(i+2) - bins_edge[i];
      cent     [i] = h_cent->GetAt(i+1);
    }

    bool draw_scale_unc = false;
    bool draw_pdf_unc = false;

    // Loop over other directories
    for (size_t di=1; di<ndirs; ++di) {
      TDirectory * const dir = dirs[di];
      const string dir_name(dir->GetName());
      TH1D * const h = get<TH1D>(dir,h_name.c_str());

      if (dir_name.find("_up_") != string::npos) { // pdf up
        draw_pdf_unc = true;
        for (size_t i=0;i<nbins;++i)
          pdf_hi[i] = h->GetAt(i+1) - cent[i];

      } else if (dir_name.find("_down_") != string::npos) { // pdf down
        draw_pdf_unc = true;
        for (size_t i=0;i<nbins;++i)
          pdf_lo[i] = cent[i] - h->GetAt(i+1);

      } else { // scale variation
        if(!draw_scale_unc) draw_scale_unc = true;
        for (size_t i=0;i<nbins;++i) {
          Double_t x = h->GetAt(i+1) - cent[i];
          if (x>0.) {
            if (scales_hi[i]<x) scales_hi[i] = x;
          } else {
            x = -x;
            if (scales_lo[i]<x) scales_lo[i] = x;
          }
        }
      }
    }

    for (size_t i=0;i<nbins;++i) {
      const Double_t unit = N_events*bins_wdth[i];
      cent     [i] /= unit;
      h_cent->SetAt(cent[i],i+1);
      pdf_lo   [i] /= unit;
      pdf_hi   [i] /= unit;
      scales_hi[i] /= unit;
      scales_lo[i] /= unit;
    }

    // Draw *********************************************************
    TLine zero;
    zero.SetLineStyle(7);

    TGraphAsymmErrors g_scales(nbins,bins_edge.data(),cent.data(),
                                     0,bins_wdth.data(),
                                     scales_lo.data(),scales_hi.data()),
                      g_pdf   (nbins,bins_edge.data(),cent.data(),
                                     0,bins_wdth.data(),
                                     pdf_lo.data(),pdf_hi.data());

    TAxis *xa = g_scales.GetXaxis();
    TAxis *ya = g_scales.GetYaxis();

    xa->SetRangeUser(bins_edge[0],bins_edge.back()+bins_wdth.back());
    g_scales.SetTitle( title.size() ? (h_name+' '+title).c_str()
                                    : h_name.c_str() );

/*
    if (draw_pdf_unc) {
      ya->SetRangeUser(
        min(g_scales.GetMinimum(),g_pdf.GetMinimum())*0.95,
        max(g_scales.GetMaximum(),g_pdf.GetMaximum())*1.1
      );
    }
*/

    ya->SetTitleOffset(1.3);
    switch (h_var) {
      case var::pT:
        xa->SetTitle("p_{T}, GeV");
        ya->SetTitle("d#sigma/dp_{T}, pb/GeV");
        break;
      case var::HT:
        xa->SetTitle("HT, GeV");
        ya->SetTitle("d#sigma/dHT, pb/GeV");
        break;
      case var::mass:
        xa->SetTitle("m, GeV");
        ya->SetTitle("d#sigma/dm, pb/GeV");
        break;
      case var::y:
        xa->SetTitle("y");
        ya->SetTitle("d#sigma/dy, pb");
        break;
      case var::phi:
        xa->SetTitle("#phi, rad");
        ya->SetTitle("d#sigma/d#phi, pb/rad");
        break;
      case var::tau:
        xa->SetTitle("#tau, GeV");
        ya->SetTitle("d#sigma/d#tau, pb/GeV");
        break;
      default:
        ya->SetTitle("#sigma, pb");
    }

    g_scales.SetFillColorAlpha(2,0.5);
    // g_scales .SetLineColor(10);
    // g_scales .SetFillStyle(3004);
    g_scales.SetLineWidth(0);
    // g_scales .SetMarkerColor(4);
    // g_scales .SetMarkerStyle(21);
    g_scales.Draw("a2");

    if (draw_pdf_unc) {
      g_pdf.SetFillColorAlpha(4,0.5);
      // g_pdf.SetLineColor(10);
      // g_pdf.SetFillStyle(3005);
      g_pdf.SetLineWidth(0);
      g_pdf.Draw("2");
    }

    h_cent->Sumw2(false);
    h_cent->SetLineWidth(2);
    if (draw_scale_unc || draw_pdf_unc)
      h_cent->SetLineColor(1);
    else
      h_cent->SetLineColor(602);
    h_cent->Draw("same");

    TLegend leg( draw_pdf_unc ? 0.69 : 0.72, 0.835,0.91,0.89);
    leg.SetNColumns(2);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetMargin(0.4);
    leg.SetColumnSeparation(-0.1);
    if (draw_scale_unc) leg.AddEntry(&g_scales,"Scale");
    if (draw_pdf_unc)   leg.AddEntry(&g_pdf,   "PDF");
    leg.Draw();

    const Double_t lbl_y2 =
      ( draw_scale_unc || draw_pdf_unc ? 0.82 : 0.88 );

    TLatex cs_lbl(0.73,lbl_y2, sigma_prt(sigma,sigma_prec).c_str());
    cs_lbl.SetNDC();
    cs_lbl.SetTextAlign(13);
    cs_lbl.SetTextFont(42);
    cs_lbl.SetTextSize(0.035);
    cs_lbl.Draw();

    TLatex N_lbl(0.73,lbl_y2-0.04, Form("Ent: %ld",lround(h_cent->GetEntries())));
    N_lbl.SetNDC();
    N_lbl.SetTextAlign(13);
    N_lbl.SetTextFont(42);
    N_lbl.SetTextSize(0.035);
    N_lbl.Draw();

    TText jet_alg_lbl(0.73,lbl_y2-0.08,("Jet alg: "+combine(jet_alg)).c_str());
    jet_alg_lbl.SetNDC();
    jet_alg_lbl.SetTextAlign(13);
    jet_alg_lbl.SetTextFont(42);
    jet_alg_lbl.SetTextSize(0.035);
    jet_alg_lbl.Draw();

    unique_ptr<TLatex> pdf_lbl;
    if (pdf_name.size()) {
      pdf_lbl.reset( new TLatex(0.73,lbl_y2-0.12, ("PDF: "+combine(pdf_name)).c_str()) );
      pdf_lbl->SetNDC();
      pdf_lbl->SetTextAlign(13);
      pdf_lbl->SetTextFont(42);
      pdf_lbl->SetTextSize(0.035);
      pdf_lbl->Draw();
    }

    unique_ptr<TText> u_lbl, o_lbl;
    if (of_lim_set) {
      if (sigma_u/sigma > of_lim) {
        stringstream ss;
        ss << setprecision(2) << 100.*sigma_u/sigma << "% underflow";
        cout << "\033[33m" << ss.str() << "\033[0m" << endl;
        u_lbl.reset( new TText(0.3,0.45,ss.str().c_str()) );
        u_lbl->SetNDC();
        u_lbl->SetTextAlign(13);
        u_lbl->SetTextColor(2);
        u_lbl->SetTextSize(0.1);
        u_lbl->Draw();
      }
      if (sigma_o/sigma > of_lim) {
        stringstream ss;
        ss << setprecision(2) << 100.*sigma_o/sigma << "% overflow";
        cout << "\033[33m" << ss.str() << "\033[0m" << endl;
        o_lbl.reset( new TText(0.3,0.55,ss.str().c_str()) );
        o_lbl->SetNDC();
        o_lbl->SetTextAlign(13);
        o_lbl->SetTextColor(2);
        o_lbl->SetTextSize(0.1);
        o_lbl->Draw();
      }
    }

    if (ya->GetXmin()<0. && 0.<ya->GetXmax())
      zero.DrawLine(bins_edge.front(),0,bins_edge.back()+bins_wdth.back(),0);

    // **************************************************************

    canv.SaveAs(fout_name.c_str());
  }

  canv.SaveAs((fout_name+']').c_str());

  delete fin;

  return 0;
}
