#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <cctype>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <tuple>
#include <memory>
#include <functional>

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

class splicer {
  vector<function<void(string&, const vector<string>* props)>> parts;

public:
  void operator()(string& str) { // consumes the input string
    size_t sep, i;
    while ((sep = str.find('\\'))!=string::npos) {
      const string sub = str.substr(0,sep);
      if (sub.size())
        parts.emplace_back([sub](string& s, const vector<string>* props){
          s += sub;
        });
      str.erase(0,sep+1);

      int id;
      try { id = std::stoi(str,&i); }
      catch (...) {
        throw runtime_error("--group: \\ not followed by a number");
      }
      parts.emplace_back([id](string& s, const vector<string>* props){
        s += props->at(id < 0 ? props->size()+id : id );
      });
      str.erase(0,i);
    }
    if (str.size())
      parts.emplace_back([str](string& s, const vector<string>* props){
        s += str;
      });
  }
  string operator()(const vector<string>* props) const {
    string str;
    for (auto &f : parts) f(str,props);
    return str;
  }
};

bool noN;
unsigned ngroups = 1; // 0 don't sort, 1 sort
splicer group_name;
map<unsigned,boost::regex> rex;
using tokenizer = boost::tokenizer<boost::char_separator<char>>;
map<unsigned,boost::char_separator<char>> tok;
using group_map_t = unordered_map<
  string, pair<unsigned,vector<pair<TH1D*,vector<string>>>>
>;
group_map_t groups;

void read_hists(TDirectory* dir, vector<string>& dir_names) {
  dir_names.emplace_back(dir->GetName());
  const size_t nd = dir_names.size();

  TIter nextkey(dir->GetListOfKeys());
  TKey *key;
  while ((key = static_cast<TKey*>(nextkey()))) {
    TObject *obj = key->ReadObj();
    if (obj->InheritsFrom(TH1D::Class())) {

      TH1D *h = static_cast<TH1D*>(obj);

      if (!noN) // Skip the N histogram
        if (dir->InheritsFrom(TFile::Class()) && !strcmp(h->GetName(),"N"))
          continue;

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

      // Add to the map -----------------------------------
      // if (props.size() <= group_by)
      //   throw runtime_error("not enough properties");
      // string group_prop = props[group_by];
      // props.erase(props.begin()+group_by);
      string group_prop = group_name(&props);
      const bool exists = groups.count(group_prop);
      auto &g = groups[group_prop];
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
  vector<string> fin_name;
  string fout_name;

  try {
    bool sort_groups;
    vector<pair<unsigned,string>> rex_str;
    vector<pair<unsigned,string>> tok_str;
    // vector<pair<unsigned,string>> ign_str;
    string group_str;

    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<vector<string>>(&fin_name)->required(),
     "*input files with histograms")
    ("output,o", po::value<string>(&fout_name),
     "output pdf plots")
    ("group,g", po::value<string>(&group_str)->default_value("\\0"),
     "group by propery\n"
     "\\0: hist, \\1: file, \\2-: dir")
    ("tokenize,t", po::value<vector<pair<unsigned,string>>>(&tok_str),
     "(i:delim) split name")
    ("regex,r", po::value<vector<pair<unsigned,string>>>(&rex_str),
     "(i:regex) apply regular expression to name")
    ("sort,s", po::bool_switch(&sort_groups),
     "alphabetically sort order of groups")
    ("no-N", po::bool_switch(&noN),
     "no N histogram")
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
    for (const auto& x : tok_str)
      tok.emplace(x.first, boost::char_separator<char>(x.second.c_str()));
    if (sort_groups) ngroups = 0;
    group_name(group_str);
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
    if (!noN) f.second = get<TH1D>(f.first,"N")->GetAt(1);
    read_hists(f.first);
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

  vector<set<string>> uniq;
  for (auto &it : order) {
    cout << it->first << ' ' << it->second.second.size() << endl;
    for (auto &s : it->second.second) {
      for (size_t i=0, n=s.second.size(); i<n; ++i) {
        if (uniq.size() < n) uniq.resize(n);
        uniq[i].insert(s.second[i]);
      }
    }
  }

  for (size_t i=0, n=uniq.size(); i<n; ++i)
    for (auto &x : uniq[i])
      cout << i << ": " << x << endl;

  for (auto &f : fin) delete f.first;
  return 0;
}
