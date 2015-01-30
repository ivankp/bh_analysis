#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
//#include <algorithm>

#include <boost/program_options.hpp>
#include<boost/tokenizer.hpp>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TLatex.h>
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
  TH1 * const h_N = get<TH1>(fin,"N");
  const Double_t N = h_N->GetBinContent(1);
  cout << "Events: " << h_N->GetEntries() << endl;

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
      cout << t << endl;
      if (!t.substr(0,3).compare("PDF")) pdf_name.emplace(t,3);
      // TODO: Add Jet algorithm label
    }
  }
  cout << endl;
  const size_t ndirs = dirs.size();
  
  vector<TH1*> central;
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
    } else if ( (h_name.find("_deltay") != string::npos) ||
                (h_name.find("_y") != string::npos) ) {
      h_var = var::y;
    } else if ( (h_name.find("_deltaphi") != string::npos) ||
                (h_name.find("_phi") != string::npos) ) {
      h_var = var::phi;
    } else if (h_name.find("_HT") != string::npos) {
      h_var = var::HT;
    } else if (h_name.find("_tau") != string::npos) {
      h_var = var::tau;
    }

    // Cross section
    const Double_t sigma   = ( h_name.find("_N_incl") != string::npos
                           ? h_cent->GetBinContent(1)
                           : h_cent->Integral(0,nbins+1) );
    const Double_t sigma_u = h_cent->GetBinContent(0);
    const Double_t sigma_o = h_cent->GetBinContent(nbins+1);
    
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
      cent  [i]    = h_cent->GetBinContent(i+1);
    }
    
    // Loop over other directories
    for (size_t di=1; di<ndirs; ++di) {
      TDirectory * const dir = dirs[di];
      const string dir_name(dir->GetName());
      TH1 * const h = get<TH1>(dir,h_name.c_str());
      
      if (dir_name.substr(dir_name.rfind('_')+1)=="up") { // pdf up
        for (size_t i=0;i<nbins;++i) 
          pdf_hi[i] = h->GetBinContent(i+1) - cent[i];
        
      } else if (dir_name.substr(dir_name.rfind('_')+1)=="down") { // pdf down
        for (size_t i=0;i<nbins;++i) 
          pdf_lo[i] = cent[i] - h->GetBinContent(i+1);
        
      } else { // scale variation
        for (size_t i=0;i<nbins;++i) {
          Double_t x = h->GetBinContent(i+1) - cent[i];
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
      const Double_t unit = ( h_var==var::N ? N : N*bins_wdth[i]);
      cent     [i] /= unit;
      h_cent->SetBinContent(i+1,cent[i]);
      pdf_lo   [i] /= unit;
      pdf_hi   [i] /= unit;
      scales_hi[i] /= unit;
      scales_lo[i] /= unit;
    }
    
    // Draw *********************************************************
    
    TGraphAsymmErrors g_scales (nbins,bins_edge.data(),cent.data(),
                                      0,bins_wdth.data(),
                                      scales_lo.data(),scales_hi.data()),
                      g_pdf_unc(nbins,bins_edge.data(),cent.data(),
                                      0,bins_wdth.data(),
                                      pdf_lo.data(),pdf_hi.data());

    g_scales.GetXaxis()
      ->SetRangeUser(bins_edge[0],bins_edge.back()+bins_wdth.back());
    g_scales.SetTitle( title.size() ? (h_name+' '+title).c_str()
                                    : h_name.c_str() );
    TAxis *xa = g_scales.GetXaxis();
    TAxis *ya = g_scales.GetYaxis();
    ya->SetTitleOffset(1.3);
    switch (h_var) {
      case var::pT:
        xa->SetTitle("pT, GeV");
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

    g_pdf_unc.SetFillColorAlpha(4,0.5);
    // g_pdf_unc.SetLineColor(10);
    // g_pdf_unc.SetFillStyle(3005);
    g_pdf_unc.SetLineWidth(0);
    g_pdf_unc.Draw("2");

    h_cent->Sumw2(false);
    h_cent->SetLineWidth(2);
    h_cent->SetLineColor(1);
    h_cent->Draw("same");

    TLegend leg(0.72,0.75,0.89,0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(&g_scales, "Scale unc");
    leg.AddEntry(&g_pdf_unc,"PDF unc");
    leg.Draw();

    TLatex cs_lbl(0.73,0.74, sigma_prt(sigma,3).c_str());
    cs_lbl.SetNDC();
    cs_lbl.SetTextAlign(13);
    cs_lbl.SetTextFont(42);
    cs_lbl.SetTextSize(0.035);
    cs_lbl.Draw();

    TLatex N_lbl(0.73,0.70, Form("Entries: %.2e",h_cent->GetEntries()));
    N_lbl.SetNDC();
    N_lbl.SetTextAlign(13);
    N_lbl.SetTextFont(42);
    N_lbl.SetTextSize(0.035);
    N_lbl.Draw();

    TText jet_alg_lbl(0.73,0.66,("Jet alg: "+combine(jet_alg)).c_str());
    jet_alg_lbl.SetNDC();
    jet_alg_lbl.SetTextAlign(13);
    jet_alg_lbl.SetTextFont(42);
    jet_alg_lbl.SetTextSize(0.035);
    jet_alg_lbl.Draw();

    TLatex pdf_lbl(0.73,0.62, ("PDF: "+combine(pdf_name)).c_str() );
    pdf_lbl.SetNDC();
    pdf_lbl.SetTextAlign(13);
    pdf_lbl.SetTextFont(42);
    pdf_lbl.SetTextSize(0.035);
    pdf_lbl.Draw();

    if (of_lim_set) {
      if (sigma_u/sigma > of_lim) {
        stringstream ss;
        ss << setprecision(2) << 100.*sigma_u/sigma << "% underflow";
        cout << "\033[33m" << ss.str() << "\033[0m" << endl;
        TText *u_lbl = new TText(0.3,0.45,ss.str().c_str());
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
        TText *o_lbl = new TText(0.3,0.55,ss.str().c_str());
        o_lbl->SetNDC();
        o_lbl->SetTextAlign(13);
        o_lbl->SetTextColor(2);
        o_lbl->SetTextSize(0.1);
        o_lbl->Draw();
      }
    }
    
    // **************************************************************
    
    canv.SaveAs(fout_name.c_str());
  }

  canv.SaveAs((fout_name+']').c_str());
  
  delete fin;

  return 0;
}
