#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <algorithm>
#include <initializer_list>
#include <memory>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TText.h>
#include <TLatex.h>
#include <TAxis.h>

#include "hist_range.hh"

using namespace std;
namespace po = boost::program_options;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

namespace std {
  template<typename A, typename B>
  istream& operator>>(istream& is, pair<A,B>& p) {
    /*
    stringstream a, b;
    bool first = true;
    char c;
    while (is >> c) {
      if (first && c==':') {
        first = false;
        continue;
      }
      (first ? a : b) << c;
    }
    if (first) {
      is.setstate(ios::failbit);
    } else {
      a >> p.first;
      b >> p.second;
      cout << p.first << "  " << p.second << endl;
    }
    */

    string str;
    is >> str;
    const size_t sep = str.find(':');
    if (sep==string::npos || sep==0 || sep==str.size()-1) {
      is.setstate(ios::failbit);
    } else {
      stringstream(str.substr(0,sep)) >> p.first;
      stringstream(str.substr(sep+1)) >> p.second;
    }

    return is;
  }
}

template<class T = TObject>
inline T* get(TDirectory* dir, const char* name) {
  T *obj = nullptr;
  dir->GetObject(name,obj);
  if (obj) return obj;
  else {
    cerr << "No object " << name << " in " << dir->GetName() << endl;
    exit(1);
  }
}

bool ior_contains(const string& in, const initializer_list<string>& l) {
  for (auto &s : l) {
    auto it = std::search(
      in.begin(), in.end(), s.begin(), s.end(),
      [](char a, char b){ return std::tolower(a) == std::tolower(b); }
    );
    if (it != in.end()) return true;
  }
  return false;
}

template<class T> string cat(const T& container, const string& delim) {
  string s;
  auto it = container.begin();
  const auto end = container.end();
  if (it==end) return s;
  s += *(it++);
  for (;it!=end;++it) (s += delim) += *it;
  return s;
}

template<typename T> void backadder(string& str, const T& add) { str += add; }
template<> void backadder<set<string>>(string& str, const set<string>& add) {
  str += cat(add, ", ");
}

template<typename T>
string substitute(const string& str, const vector<T>& parts) {
  string result;
  size_t sep, n = 0;
  const char *first = &str[0], *begin = first;
  char *end;
  while ((sep = str.find('\\',n))!=string::npos) {
    sep -= n;
    result.append(begin,sep);
    begin += sep + 1;
    if (!isdigit(*begin))
      throw runtime_error("--group: \\ not followed by a number");
    backadder<T>(result, parts.at( strtol(begin,&end,10) ));
    begin = end;
    n = begin - first;
  }
  result.append(begin,str.size()-n);
  return result;
}

class comma_numpunct : public std::numpunct<char> {
protected:
  virtual char do_thousands_sep()   const { return ','; }
  virtual std::string do_grouping() const { return "\03"; }
};

// VARS ---------------------------------------------------
bool noN;
unsigned ngroups = 1; // 0 don't sort, 1 sort
map<unsigned,boost::regex> rex, ign;
using tokenizer = boost::tokenizer<boost::char_separator<char>>;
map<unsigned,boost::char_separator<char>> tok;
using group_map_t = unordered_map<
  string, pair<unsigned,vector<pair<TH1D*,vector<string>>>>
>;
group_map_t groups;
string group_str;
Double_t N_scale;
// --------------------------------------------------------

void read_hists(TDirectory* dir, vector<string>& dir_names) {
  dir_names.emplace_back(dir->GetName());
  const size_t nd = dir_names.size();

  TIter nextkey(dir->GetListOfKeys());
  TKey *key;
  while ((key = static_cast<TKey*>(nextkey()))) {
    TObject *obj = key->ReadObj();
    if (obj->InheritsFrom(TH1D::Class())) {

      TH1D *h = static_cast<TH1D*>(obj);

      if (!noN) {
        // Skip the N histogram
        if (dir->InheritsFrom(TFile::Class()) && !strcmp(h->GetName(),"N"))
          continue;

        h->Scale(N_scale);
      }

      // Parse names into properties ----------------------
      vector<string> props;
      for (size_t i=0; i<=nd; ++i) {
        auto name = [&]() {
          return (i ? dir_names[i-1] : h->GetName());
        };

        if (rex.count(i)) {
          // boost::smatch result;
          // if ( boost::regex_search(name(), result, rex[i]) )
          //   for (auto &r : result) props.emplace_back(r.first, r.second);
          // else continue;
          const string &s = name();
          boost::sregex_token_iterator it1(s.begin(), s.end(), rex[i], -1);
          static const boost::sregex_token_iterator end;
          while (it1 != end) props.emplace_back(*it1++);
        } else if (tok.count(i)) {
          for (auto &t : tokenizer(name(), tok[i]))
            props.emplace_back(t);
        } else {
          props.emplace_back(name());
        }
      }

      if (ign.size()) {
        bool ignore = false;
        for (size_t i=0, n=props.size(); i<n; ++i)
          if (ign.count(i)) {
            static boost::smatch result;
            ignore = boost::regex_search(props[i], result, ign[i]);
          }
        if (ignore) continue;
      }

      // Add to the map -----------------------------------
      string group_prop = substitute(group_str,props);
      const bool exists = groups.count(group_prop);
      auto &g = groups[move(group_prop)];
      g.second.emplace_back(h,move(props));
      if (!exists) g.first = (ngroups ? ++ngroups : 0);

    // Go into the next level of directories --------------
    } else if (obj->InheritsFrom(TDirectory::Class())) {
      read_hists(static_cast<TDirectory*>(obj), dir_names);
    }
  }
  dir_names.erase(dir_names.begin()+nd-1,dir_names.end());
}

inline void read_hists(TDirectory* dir) {
  vector<string> dir_names;
  read_hists(dir,dir_names);
}

unsigned sigma_prec;
string sigma_prt(Double_t sigma) {
  stringstream ss;
  ss << "#sigma = " << showpoint << setprecision(sigma_prec) << sigma
     << noshowpoint << " pb";
  return ss.str();
}

// ******************************************************************
constexpr Color_t color[] = { 2, 3, 4, 6, 7, 8, 9, 28, 30, 46 };

int main(int argc, char **argv)
{
  // START OPTIONS **************************************************
  vector<string> fin_name, labels;
  string fout_name, draw_opt;
  bool prt_props, prt_saved;

  try {
    bool sort_groups;
    vector<pair<unsigned,string>> rex_str;
    vector<pair<unsigned,string>> tok_str;
    vector<pair<unsigned,string>> ign_str;

    // General Options ------------------------------------
    po::options_description opt("Options");
    opt.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<vector<string>>(&fin_name)->required(),
     "*input files with histograms")
    ("output,o", po::value<string>(&fout_name),
     "output pdf plots")
    ("group,g", po::value<string>(&group_str)->default_value("\\0"),
     "group by propery\n\\0: hist, \\1: file, \\2+: dir")
    ("tokenize,t", po::value<vector<pair<unsigned,string>>>(&tok_str),
     "(i:delim) split name")
    ("regex,r", po::value<vector<pair<unsigned,string>>>(&rex_str),
     "(i:regex) apply regular expression to name")
    ("ignore,n", po::value<vector<pair<unsigned,string>>>(&ign_str),
     "(i:regex) ignore matching ith name value")
    ("sort,s", po::bool_switch(&sort_groups),
     "alphabetically sort order of groups")
    ("label,l", po::value<vector<string>>(&labels),
     "add a label")
    ("no-N", po::bool_switch(&noN), "no N histogram")
    ("draw", po::value<string>(&draw_opt), "TH1 Draw options")
    ("sigma-prec", po::value<unsigned>(&sigma_prec)->default_value(3),
     "number of significant digits in cross section")
    ("prt-props", po::bool_switch(&prt_props),
     "print available properties")
    ("prt-saved", po::bool_switch(&prt_saved),
     "print groups as they are saved")
    ;

    po::positional_options_description pos;
    pos.add("input",-1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(opt).positional(pos).run(), vm);
    if (argc == 1 || vm.count("help")) {
      cout << opt << endl;
      return 0;
    }
    po::notify(vm);

    if (!fout_name.size()) {
      if (fin_name.size()==1) {
        size_t first = fin_name.front().rfind('/')+1;
        if (first==string::npos) first = 0;
        fout_name = fin_name.front().substr(
          first, fin_name.front().rfind('.') - first
        ) + ".pdf";
      } else throw runtime_error("--output is required for multiple inputs");
    }
    // for (const auto& x : rex_str) rex.emplace(x);
    for (const auto& x : rex_str) rex.emplace(x.first, boost::regex(x.second));
    for (const auto& x : tok_str)
      tok.emplace(x.first, boost::char_separator<char>(x.second.c_str()));
    for (const auto& x : ign_str) ign.emplace(x.first, boost::regex(x.second));
    if (sort_groups) ngroups = 0;
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    return 1;
  }
  // END OPTIONS ****************************************************

  // Read histograms ************************************************
  vector<unique_ptr<TFile>> fin;
  fin.reserve(fin_name.size());
  for (const auto &name : fin_name) {
    TFile *f = new TFile(name.c_str(),"read");
    if (f->IsZombie()) return 1;
    fin.emplace_back(f);
  }
  for (auto &f : fin) {
    if (!noN) N_scale = 1./get<TH1D>(f.get(),"N")->GetAt(1);
    read_hists(f.get());
  }

  vector<group_map_t::iterator> order;
  order.reserve(groups.size());
  for (auto it=groups.begin(), end=groups.end(); it!=end; ++it) {
    order.push_back(it);
  }
  if (ngroups) {
    sort(order.begin(), order.end(),
      [](group_map_t::iterator a, group_map_t::iterator b) {
        return (a->second.first < b->second.first);
      }
    );
  } else {
    sort(order.begin(), order.end(),
      [](group_map_t::iterator a, group_map_t::iterator b) {
        return (a->first < b->first);
      }
    );
  }

  // fill sets of properties
  vector<set<string>> uniq;
  for (auto it : order) {
    // cout << it->first << ' ' << it->second.second.size() << endl;
    for (auto &s : it->second.second) {
      for (size_t i=0, n=s.second.size(); i<n; ++i) {
        if (uniq.size() < n) uniq.resize(n);
        uniq[i].insert(s.second[i]);
      }
    }
  }

  // print sets of properties
  if (prt_props) {
    for (size_t i=0, n=uniq.size(); i<n; ++i) {
      bool first = true;
      for (auto &x : uniq[i]) {
        if (first) {
          cout << setw(3) << i << ": ";
          first = false;
        } else cout << "     ";
        cout << x << endl;
      }
    }
    return 0;
  }

  set<size_t> rm_props;

  if (prt_props) {
    // list available properties
    for (size_t i=0, n=uniq.size(); i<n; ++i) {
      bool first = true;
      for (auto &x : uniq[i]) {
        if (first) {
          cout << setw(3) << i << ": ";
          first = false;
        } else cout << "     ";
        cout << x << endl;
      }
    }
    return 0;
  } else {
    // properties used for groups
    size_t sep, n = 0;
    const char *first = &group_str[0], *begin = first;
    char *end;
    while ((sep = group_str.find('\\',n))!=string::npos) {
      begin += (sep -= n) + 1;
      rm_props.insert( strtol(begin,&end,10) );
      n = (begin = end) - first;
    }
    // common properties
    for (size_t i=0, n=uniq.size(); i<n; ++i)
      if (uniq[i].size()==1) rm_props.insert(i);
  }

  // Draw ***********************************************************
  gStyle->SetOptStat(0);

  TLine zero;
  zero.SetLineStyle(7);

  TCanvas canv;
  canv.SaveAs((fout_name+'[').c_str());

  // string, pair<unsigned,vector<pair<TH1D*,vector<string>>>>
  for (auto it : order) {
    const auto& hists = it->second.second;
    hist_range y_range;

    TLegend leg(0.72,0.92-0.0325*hists.size(),0.95,0.92);
    TLegend sig(leg.GetX1(),2*leg.GetY1()-leg.GetY2(),
                leg.GetX2(),leg.GetY1());
    leg.SetFillColorAlpha(0,0.65);
    sig.SetFillColorAlpha(0,0.65);

    int i = 0;
    TH1D *first = nullptr;
    for (const auto& hp : hists) {
      TH1D *h = hp.first;

      const Int_t nbins = h->GetNbinsX()+1;

      // Cross section
      const Double_t sigma = (
        ior_contains(h->GetName(), {"_N_incl","inclusive"})
        ? h->GetAt(1) : h->Integral(0,nbins)
      );

      for (Int_t i=0; i<nbins; ++i)
        h->SetAt(h->GetAt(i)/(h->GetBinWidth(i)),i);

      y_range(h->GetMinimum(),h->GetMaximum());

      if (i==0) {
        first = h;
        h->Draw(draw_opt.c_str());
      } else {
        static string draw_opt_same = draw_opt + "same";
        h->Draw(draw_opt_same.c_str());
      }

      h->SetLineWidth(2);
      h->SetLineColor(color[i]);
      h->SetMarkerColor(color[i]);
      h->Sumw2(false);
      string leg_name;
      for (size_t i=0, n=hp.second.size(); i<n; ++i) {
        if (rm_props.count(i)) continue;
        if (leg_name.size()) leg_name += '_';
        leg_name += hp.second[i];
      }
      leg.AddEntry(h,leg_name.c_str());
      sig.AddEntry(h,sigma_prt(sigma).c_str());

      ++i;
    }

    TLine top_corner_cover(leg.GetX1(),0.9,0.91,0.9);
    TLine right_corner_cover(0.9,sig.GetY1(),0.9,0.91);
    top_corner_cover.SetNDC();
    right_corner_cover.SetNDC();
    top_corner_cover.SetLineColor(0);
    right_corner_cover.SetLineColor(0);
    top_corner_cover.SetLineWidth(3);
    right_corner_cover.SetLineWidth(3);

    const string& title = it->first;
    const string first_name(first->GetName());
    y_range(first)->SetTitle(title.c_str());

    if (first->GetMinimum()<0. && 0.<first->GetMaximum())
      zero.DrawLine(first->GetBinLowEdge(1),0,
                    first->GetBinLowEdge(first->GetNbinsX()+1),0);

    TAxis *xa = first->GetXaxis();
    TAxis *ya = first->GetYaxis();
    ya->SetTitleOffset(1.3);

    if (ior_contains(first_name, {"_pT"})) {
      xa->SetTitle("pT, GeV");
      ya->SetTitle("d#sigma/dpT, pb/GeV");
    } else if (ior_contains(first_name, {"_mass"})) {
      xa->SetTitle("m, GeV");
      ya->SetTitle("d#sigma/dm, pb/GeV");
    } else if (ior_contains(first_name, {"_y","_dy","_deltay"})) {
      xa->SetTitle("y");
      ya->SetTitle("d#sigma/dy, pb");
    } else if (ior_contains(first_name, {"_eta","_deta","_deltaeta"})) {
      xa->SetTitle("#eta");
      ya->SetTitle("d#sigma/d#eta, pb");
    } else if (ior_contains(first_name, {"_phi","_dphi","_deltaphi"})) {
      xa->SetTitle("#phi, rad");
      ya->SetTitle("d#sigma/d#phi, pb/rad");
    } else if (ior_contains(first_name, {"_dR"})) {
      xa->SetTitle("#DeltaR");
      ya->SetTitle("d#sigma/dR, pb");
    } else if (ior_contains(first_name, {"_mT"})) {
      xa->SetTitle("mT, GeV");
      ya->SetTitle("d#sigma/dmT, pb/GeV");
    } else if (ior_contains(first_name, {"_HT"})) {
      xa->SetTitle("HT, GeV");
      ya->SetTitle("d#sigma/dHT, pb/GeV");
    } else if (ior_contains(first_name, {"_tau"})) {
      xa->SetTitle("#tau, GeV");
      ya->SetTitle("d#sigma/d#tau, pb/GeV");
    } else ya->SetTitle("d#sigma, pb");

    top_corner_cover.Draw();
    right_corner_cover.Draw();
    leg.Draw();
    sig.Draw();

    Double_t lbl_x = 0.73, lbl_y = sig.GetY1()+0.03;

    vector<unique_ptr<TLatex>> lbl;
    lbl.reserve(labels.size());
    for (const string& lbl_name : labels) {
      lbl.emplace_back( new TLatex( lbl_x, lbl_y-=0.045,
        substitute(lbl_name,uniq).c_str()
      ) );
      lbl.back()->SetNDC();
      lbl.back()->SetTextAlign(13);
      lbl.back()->SetTextFont(42);
      lbl.back()->SetTextSize(0.035);
      lbl.back()->Draw();
    }

    if (prt_saved) cout << title << endl;
    canv.SaveAs(fout_name.c_str());
  }

  canv.SaveAs((fout_name+']').c_str());

  return 0;
}
