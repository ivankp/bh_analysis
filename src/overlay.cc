#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <list>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <tuple>
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
#include <TAxis.h>

#include "hist_range.hh"

using namespace std;
namespace po = boost::program_options;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

namespace std {
  template<typename A, typename B>
  istream& operator>>(istream &is, pair<A,B> &p) {
    stringstream a, b;
    bool first = true;
    char c;
    while (is >> c) {
      if (c==':') first = false;
      (first ? a : b) << c;
    }
    if (first) {
      is.setstate(ios::failbit);
    } else {
      a >> p.first;
      b >> p.second;
    }
    return is;
  }
}

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

unsigned ngroups = 0; // 0 don't sort, 1 sort
unsigned group_by = 0;
map<unsigned,boost::regex> rex;
using group_map_t = unordered_map<
  string, pair<unsigned,vector<pair<TH1D*,vector<string>>>>
>;
group_map_t groups;

void read_hists(TDirectory* dir, vector<string>* dir_names = new vector<string>()) {
  dir_names->emplace_back(dir->GetName());
  const size_t nd = dir_names->size();

  TIter nextkey(dir->GetListOfKeys());
  TKey *key;
  while ((key = static_cast<TKey*>(nextkey()))) {
    TObject *obj = key->ReadObj();
    if (obj->InheritsFrom(TH1D::Class())) {

      TH1D *h = static_cast<TH1D*>(obj);

      vector<string> props;
      if (rex.count(0)) {
        // process regex
      } else {
        props.emplace_back(h->GetName());
      }
      for (size_t i=0; i<nd; ++i) {
        if (rex.count(i+1)) {
          // process regex
        } else {
          props.emplace_back((*dir_names)[i]);
        }
      }

      if (props.size() <= group_by)
        throw runtime_error("not enough properties");
      string group_prop = props[group_by];
      props.erase(props.begin()+group_by);
      const bool exists = groups.count(group_prop);
      auto &g = groups[group_prop];
      g.second.emplace_back(h,move(props));
      if (!exists) g.first = (ngroups ? ++ngroups : 0);

    } else if (obj->InheritsFrom(TDirectory::Class())) {
      read_hists(static_cast<TDirectory*>(obj), dir_names);
    }
  }
  dir_names->erase(dir_names->begin()+nd-1,dir_names->end());
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
  vector<string> fin_name;
  string fout_name;
  bool sort_groups;

  try {
    vector<pair<unsigned,string>> rex_str;

    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<vector<string>>(&fin_name)->required(),
     "*input files with histograms")
    ("output,o", po::value<string>(&fout_name),
     "output pdf plots")
    ("regex,r", po::value<vector<pair<unsigned,string>>>(&rex_str),
     "regular expressions")
    ("sort,s", po::bool_switch(&sort_groups),
     "alphabetically sort order of groups")
    ("sigma-prec", po::value<unsigned>(&sigma_prec)->default_value(3),
     "number of significant digits in cross section")
    ;

    po::positional_options_description pos;
    pos.add("input",-1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(all_opt).positional(pos).run(), vm);
    if (argc == 1 || vm.count("help")) {
      cout << all_opt << endl;
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
    if (sort_groups) ngroups = 1;
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    return 1;
  }
  // END OPTIONS ****************************************************

  // Read histograms ************************************************
  vector<pair<TFile* const,Double_t>> fin;
  for (const auto& name : fin_name) {
    TFile *f = new TFile(name.c_str(),"read");
    if (f->IsZombie()) return 1;
    fin.emplace_back(f,0);
  }
  for (auto& f : fin) {
    f.second = get<TH1D>(f.first,"N")->GetAt(1);
    read_hists(f.first);
  }

  vector<group_map_t::iterator> order;
  order.reserve(groups.size());
  for (auto it=groups.begin(), end=groups.end(); it!=end; ++it) {
    order.push_back(it);
  }
  if (sort_groups) {
    sort(order.begin(), order.end(),
      [](group_map_t::iterator a, group_map_t::iterator b) {
        return (a->first < b->first);
      }
    );
  } else {
    sort(order.begin(), order.end(),
      [](group_map_t::iterator a, group_map_t::iterator b) {
        return (a->second.first < b->second.first);
      }
    );
  }

  for (auto &it : order) {
    cout << it->first << ' ' << it->second.second.size() << endl;
    for (auto &s : it->second.second) {
      for (auto &v : s.second)
        cout << "    " << v << endl;
      cout << endl;
    }
  }

  for (auto &f : fin) delete f.first;

  return 0;
}
