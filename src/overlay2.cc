#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <algorithm>

#include <boost/program_options.hpp>
//#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
/*
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TText.h>
#include <TLatex.h>
#include <TAxis.h>
*/

using namespace std;
namespace po = boost::program_options;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

struct hist {
  TH1* h;
  string name;
  vector<vector<string>> tokens; // 0 - hist, 1 - file, 2... - dir
  hist(TH1* h, const char* name, decltype(tokens)&& tokens)
  : h(h), name(name), tokens(tokens) { }
};

list<hist> hists; // container of all histograms
using hiter = decltype(hists)::iterator;

class regex {
  enum { match, ignore, substitute, tokenize } action;
  int position; // 0 - hist, 1 - file, 2... - dir, -1 - composed name
  int modifiers; // 1 - g

  boost::regex *re1, *re2;

public:
  regex(): action(match), position(0), modifiers(0),
           re1(nullptr), re2(nullptr) { }
  ~regex() {
    delete re1;
    delete re2;
  }

  friend istream& operator>>(istream& s, regex& re) {
    int i = 0;
    char c;
    bool ok;
    string str;
    // read
    for (;;) {
      if (!( ok = (bool)s.get(c) )) c = '/';
      if (c=='/') {
        if (str.size()>0 && str.back()!='\\') {
          switch (i) {
            case 0: {
              if (str=="m") re.action = match;
              else if (str=="i") re.action = ignore;
              else if (str=="s") re.action = substitute;
              else if (str=="t") re.action = tokenize;
              else throw runtime_error("Unexpected regex action "+str);
            } break;
            case 1: {
              if (str=="f") re.position = 1;
              else if (str=="h") re.position = 0;
              else if (str[0]=='d') {
                if (str.size()>1) re.position = atoi(str.c_str()+1)+1;
                else re.position = 2;
              } else if (str=="n") re.position = -1;
              else throw runtime_error("Unexpected regex position "+str);
            } break;
            case 2: re.re1 = new boost::regex(str); break;
            case 3: re.re2 = new boost::regex(str); break;
            case 4: {
              if (str=="g") re.modifiers = 1;
              else throw runtime_error("Unexpected regex modifier "+str);
            } break;
            default: break;
          }
        }
        ++i;
        test(str)
        str = "";
      } else str += c;
      if (!ok) break;
    }
    s.clear(); // std::basic_istream::get() sets failbit

    // validate
    switch (re.action) {
      case match: {
        if (!re.re1) throw runtime_error("Matching uses <re1>");
        if ( re.re2) throw runtime_error("Matching does not use <re2>");

      } break;
      case ignore: {


      } break;
      case substitute: {


      } break;
      case tokenize: {


      } break;
    }

    return s;
  }

  bool operator()(hiter h) const {
    cout << '.';
    switch (action) {
      case match:
      case ignore: {
        //boost::cmatch what;
        const bool matched = boost::regex_match(h->name, *re1);
        if ( (action==match) ^ matched ) {
          hists.erase(h);
          return true;
        }

      } break;
      case substitute: {


      } break;
      case tokenize: {


      } break;
    }
    return false;
  }
};

void get_hists(TDirectory* dir) {
  TIter nextkey(dir->GetListOfKeys());
  TKey *key;
  while ((key = static_cast<TKey*>(nextkey()))) {
    TObject *obj = key->ReadObj();
    if (obj->InheritsFrom(TH1::Class())) {

      TH1* hist = static_cast<TH1*>(obj);
      vector<vector<string>> tokens;
      for (TDirectory* d=dir; d; d=d->GetMotherDir())
        tokens.emplace_back(vector<string>{d->GetName()});
      tokens.emplace_back(vector<string>{hist->GetName()});
      reverse(tokens.begin(),tokens.end());
      hists.emplace_back(hist,hist->GetName(),move(tokens));

    } else if (obj->InheritsFrom(TDirectory::Class())) {
      get_hists(static_cast<TDirectory*>(obj));
    }
  }
}

int main(int argc, char **argv)
{
  vector<string> ifname;
  string ofname, cfname;
  vector<regex> re;

  // options ---------------------------------------------------
  try {
    po::options_description desc("Options");
    desc.add_options()
      ("input,i", po::value(&ifname)->multitoken()->required(),
       "input root files with TH1 histograms")
      ("output,o", po::value(&ofname),
       "output pdf file")
      ("conf,c", po::value(&cfname),
       "configuration file")

      ("regex,r", po::value(&re),
       "regular expression to modify names\n"
       "format: <action>/<position>/<re1>/<re2>/<modifiers>\n"
       "actions:\n"
       "\tm - match (default)\n"
       "\ti - ignore\n"
       "\ts - substitute\n"
       "\tt - tokenize\n"
       "positions:\n"
       "\tf - file name (d0)\n"
       "\td - directory name (d1)\n"
       "\th - histogram name (default)\n"
       "\tn - composed name\n"
       "modifiers:\n"
       "\tg - all occurances\n")
    ;
    
    po::positional_options_description pos;
    pos.add("input",-1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(desc).positional(pos).run(), vm);
    if (argc == 1) {
      cout << desc << endl;
      return 0;
    }
    if (vm.count("conf")) {
      po::store( po::parse_config_file<char>(
        vm["conf"].as<string>().c_str(), desc), vm);
    }
    po::notify(vm);
  } catch (exception& e) {
    cerr << "program_options: " <<  e.what() << endl;
    return 1;
  }
  // end options ---------------------------------------------------

  vector<TFile*> files;
  for (const auto& f : ifname) {
    test(f)
    files.push_back(new TFile(f.c_str(),"read"));
    test(files.back()->GetName())
    get_hists(files.back());
  }

  if (re.size()==0) { // name by non-overlapping tokens
    if (ifname.size()==1) { // compare directories

    } else { // default to comparing files

    }
  } else { // use regular expression for grouping
    for (auto h=hists.begin(); h!=hists.end();) {
      bool erase = false;
      for (const auto& r : re)
        if ((erase = r(h))) break;
      if (!erase) ++h;
    }
  }

  for (const auto& h : hists) {
    cout << h.name << endl;
    for (const auto& a : h.tokens)
      for (const auto& b : a)
        cout << "    " << b << endl;
  }

  for (auto* f : files) delete f;

  return 0;
}
