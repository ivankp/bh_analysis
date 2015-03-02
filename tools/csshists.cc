#include "csshists.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>

#include <boost/regex.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>

#include <TH1.h>

using namespace std;

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

// Properties *******************************************************

enum prop_t {
  kClass,
  kBins,
  kFillColor,
  kFillColorAlpha,
  kFillStyle,
  kXLabelColor,
  kYLabelColor,
  kXLabelFont,
  kYLabelFont,
  kXLabelOffset,
  kYLabelOffset,
  kXLabelSize,
  kYLabelSize,
  kLineColor,
  kLineColorAlpha,
  kLineStyle,
  kLineWidth,
  kMarkerColor,
  kMarkerColorAlpha,
  kMarkerSize,
  kMarkerStyle,
  kMaximum,
  kMinimum,
  kOption,
  kStats,
  kTickLength,
  kTitle,
  kTitleFont,
  kTitleOffset,
  kTitleSize,
  kXTitle,
  kYTitle,
  kSumw2
};

#define DEFMAKER2(name,type) \
  pair<prop_t,const prop*> make_prop_##name(const string& str) { \
    return make_pair(type, new const prop_##name(str)); \
  }

#define DEFMAKER(name) DEFMAKER2(name,k##name)

struct prop {
  virtual void apply(TH1* h) const =0;
  virtual ~prop() { }
};

struct prop_int: public prop {
  int x;
  prop_int(const string& str): x( atoi(str.c_str()) ) { }
  virtual ~prop_int() { }
  virtual void apply(TH1* h) const =0;
};

struct prop_string: public prop {
  string str;
  prop_string(const string& str): str(str) { }
  virtual ~prop_string() { }
  virtual void apply(TH1* h) const =0;
};

struct prop_bool: public prop {
  bool b;
  prop_bool(const string& str) {
    if (str=="1") b = true;
    else if (boost::algorithm::to_lower_copy(str)=="true") b = true;
    else b = false;
  }
  virtual ~prop_bool() { }
  virtual void apply(TH1* h) const =0;
};

struct prop_Class: public prop {
  string name;
  prop_Class(const string& name): name(name) { }
  virtual ~prop_Class() { }
  // this property is used for construction of TH1
  virtual void apply(TH1* h) const { }
};
DEFMAKER(Class)

struct prop_Bins: public prop {
  bool isrange;
  prop_Bins(bool isrange): isrange(isrange) { }
  virtual ~prop_Bins() { }
  // this property is used for construction of TH1
  virtual void apply(TH1* h) const { }
};

struct prop_Bins_range: public prop_Bins {
  Double_t nbinsx, xlow, xup;
  prop_Bins_range(const string& str): prop_Bins(true) {
    stringstream ss(str);
    ss >> nbinsx >> xlow >> xup;
  }
};
DEFMAKER2(Bins_range,kBins)

struct prop_Bin_edges: public prop_Bins {
  vector<Double_t> xbins;
  prop_Bin_edges(const string& str): prop_Bins(false) {
    stringstream ss( boost::replace_all_copy(str, ",", " ") );
    Double_t bin;
    while(ss >> bin) xbins.push_back(bin);
  }
};
DEFMAKER2(Bin_edges,kBins)

struct prop_Title: public prop_string {
  prop_Title(const string& str): prop_string(str) { }
  virtual void apply(TH1* h) const {
    h->SetTitle(str.c_str());
  }
};
DEFMAKER(Title)

struct prop_Sumw2: public prop_bool {
  prop_Sumw2(const string& str): prop_bool(str) { }
  virtual void apply(TH1* h) const {
    h->Sumw2(b);
  }
};
DEFMAKER(Sumw2)

struct prop_LineColor: public prop_int {
  prop_LineColor(const string& str): prop_int(str) { }
  virtual void apply(TH1* h) const {
    h->SetLineColor(x);
  }
};
DEFMAKER(LineColor)

struct prop_LineWidth: public prop_int {
  prop_LineWidth(const string& str): prop_int(str) { }
  virtual void apply(TH1* h) const {
    h->SetLineWidth(x);
  }
};
DEFMAKER(LineWidth)

/*

virtual void	TAttFill::SetFillColor(Color_t fcolor)
virtual void	TAttFill::SetFillColorAlpha(Color_t fcolor, Float_t falpha)
virtual void	TAttFill::SetFillStyle(Style_t fstyle)
virtual void	TH1::SetLabelColor(Color_t color = 1, Option_t* axis = "X")
virtual void	TH1::SetLabelFont(Style_t font = 62, Option_t* axis = "X")
virtual void	TH1::SetLabelOffset(Float_t offset = 0.005, Option_t* axis = "X")
virtual void	TH1::SetLabelSize(Float_t size = 0.02, Option_t* axis = "X")
virtual void	TAttLine::SetLineColorAlpha(Color_t lcolor, Float_t lalpha)
virtual void	TAttLine::SetLineStyle(Style_t lstyle)
virtual void	TAttMarker::SetMarkerColor(Color_t mcolor = 1)
virtual void	TAttMarker::SetMarkerColorAlpha(Color_t mcolor, Float_t malpha)
virtual void	TAttMarker::SetMarkerSize(Size_t msize = 1)
virtual void	TAttMarker::SetMarkerStyle(Style_t mstyle = 1)
virtual void	TH1::SetMaximum(Double_t maximum = -1111)
virtual void	TH1::SetMinimum(Double_t minimum = -1111)
virtual void	TH1::SetOption(Option_t* option = " ")
virtual void	TH1::SetStats(Bool_t stats = kTRUE)
virtual void	TH1::SetTickLength(Float_t length = 0.02, Option_t* axis = "X")
virtual void	TH1::SetTitleFont(Style_t font = 62, Option_t* axis = "X")
virtual void	TH1::SetTitleOffset(Float_t offset = 1, Option_t* axis = "X")
virtual void	TH1::SetTitleSize(Float_t size = 0.02, Option_t* axis = "X")
virtual void	TH1::SetXTitle(const char* title)
virtual void	TH1::SetYTitle(const char* title)

*/

// Map strings to property makers ***********************************

#define add_prop_type2(str,name) \
  _map.insert(make_pair(str, &make_prop_##name));

#define add_prop_type(name) add_prop_type2(#name,name)

class prop_factory_map {
  typedef boost::unordered_map< string,
    pair<prop_t,const prop*> (*)(const string& str) > map_t;

  map_t _map;

public:
  prop_factory_map() {
    add_prop_type2("Bins",    Bins_range)
    add_prop_type2("BinEdges",Bin_edges)
    add_prop_type(Class)
    add_prop_type(Title)
    // add_prop_type(FillColor)
    // add_prop_type(FillColorAlpha)
    // add_prop_type(FillStyle)
    // add_prop_type(XLabelColor)
    // add_prop_type(YLabelColor)
    // add_prop_type(XLabelFont)
    // add_prop_type(YLabelFont)
    // add_prop_type(XLabelOffset)
    // add_prop_type(YLabelOffset)
    // add_prop_type(XLabelSize)
    // add_prop_type(YLabelSize)
    add_prop_type(LineColor)
    // add_prop_type(LineColorAlpha)
    // add_prop_type(LineStyle)
    add_prop_type(LineWidth)
    // add_prop_type(MarkerColor)
    // add_prop_type(MarkerColorAlpha)
    // add_prop_type(MarkerSize)
    // add_prop_type(MarkerStyle)
    // add_prop_type(Maximum)
    // add_prop_type(Minimum)
    // add_prop_type(Option)
    // add_prop_type(Stats)
    // add_prop_type(TickLength)
    // add_prop_type(TitleFont)
    // add_prop_type(TitleOffset)
    // add_prop_type(TitleSize)
    // add_prop_type(XTitle)
    // add_prop_type(YTitle)
    add_prop_type(Sumw2)
  }

  pair<prop_t,const prop*> operator()(const string& str) const {
    // separate property: value
    const size_t sep = str.find_first_of(":");

    // get property type and maker function
    string opt = str.substr(0,sep);
    boost::algorithm::trim(opt);

    map_t::const_iterator maker = _map.find(opt);
    if (maker==_map.end()) {
      throw std::runtime_error("Undefined histogram option: "+opt);
    }

    // create property object and return with type
    return maker->second(
      boost::algorithm::trim_copy( str.substr(sep+1) )
    );
  }

} const mk_prop;

// PIMPL ************************************************************

struct csshists::impl {

  typedef boost::unordered_map<prop_t, const prop*> prop_map;

  typedef vector< pair<const boost::regex*,prop_map*> > rules_t;
  typedef rules_t::iterator iter;
  rules_t rules;

  ~impl() {
    for (size_t i=0,n=rules.size();i<n;++i) {
      delete rules[i].first;
      for (prop_map::iterator it=rules[i].second->begin(),
           end=rules[i].second->end();it!=end;++it)
      {
        delete it->second;
      }
      delete rules[i].second;
    }
  }
};

// Constructor ******************************************************

csshists::csshists(const string& cssfilename): _impl( new impl )
{
  // Read CSS file
  ifstream css(cssfilename.c_str());
  char c;
  bool brak=false, q1=false, q2=false;
  vector< pair< string, vector<string> > > rule_str(1);
  while ( css.get(c) ) {
    if (brak) {
      if (!q2) if (c=='\'') q1 = !q1;
      if (!q1) if (c=='\"') q2 = !q2;

      if (!(q1||q2)) {
        if (c=='}') {
          brak=false;
          // remove whitespace from last property string
          boost::algorithm::trim(rule_str.back().second.back());
          // ref to last rule string
          string& last_rule = rule_str.back().first;
          // remove whitespace from last rule string
          boost::algorithm::trim(last_rule);
          // unquote last rule string
          if ( (last_rule[0]=='\''&&(*(last_rule.end()-1))=='\'') ||
               (last_rule[0]=='\"'&&(*(last_rule.end()-1))=='\"') )
            last_rule = last_rule.substr(1,last_rule.size()-2);
          // start new rule
          rule_str.push_back(pair< string, vector<string> >());
          continue;
        }
        if (c=='{') { throw runtime_error(
          "consecutive { in file \""+cssfilename+"\""
        ); }
      }

      if (c==';') {
        // remove whitespace from last property string
        boost::algorithm::trim(rule_str.back().second.back());
        // start new property string
        rule_str.back().second.push_back(string());
      } else {
        rule_str.back().second.back() += c;
      }

    } else {
      if (c=='{') {
        brak=true;
        // start first property string
        rule_str.back().second.push_back(string());
        continue;
      }
      if (c=='}') { throw runtime_error(
        "} before { in file \""+cssfilename+"\""
      ); }

      rule_str.back().first += c;
    }
  }

  // Convert rules to regex and props
  _impl->rules.reserve(rule_str.size());
  for (size_t i=0,n=rule_str.size();i<n;++i) {
    if (rule_str[i].second.size()==0) continue; // skip blank rule

    _impl->rules.push_back( make_pair(
      new const boost::regex( rule_str[i].first ),
      new impl::prop_map()
    ) );
    //_impl->rules.back().second.reserve( rule_str[i].second.size() );
    for (size_t j=0,m=rule_str[i].second.size();j<m;++j) {
      const string& prop_str = rule_str[i].second[j];
      if (prop_str.size()==0) continue; // skip blank property

      _impl->rules.back().second->insert( mk_prop( prop_str ) );
    }
  }


  // Test print
  //for (size_t i=0,n=_impl->rules.size();i<n;++i) {
  //  cout << _impl->rules[i].first->str() << endl;
  //  for (impl::prop_map::iterator j=_impl->rules[i].second.begin(),
  //       m=_impl->rules[i].second.end();j!=m;++j)
  //    cout << j->first <<' '<< j->second->str << endl;
  //    //cout <<j<<". " << _impl->rules[i].second[j]->str <<';'<< endl;
  //  cout << "-----------" << endl;
  //}
}

// Make Historgram **************************************************

TH1* csshists::mkhist(const string& name) const {
  impl::prop_map props;
  TH1* h;

  boost::smatch result;
  for (impl::iter it=_impl->rules.begin(),
       end=_impl->rules.end();it!=end;++it) {
    if (boost::regex_match(name, result, *it->first)) {
      for (impl::prop_map::iterator jt=it->second->begin(),
           end2=it->second->end();jt!=end2;++jt)
        props[jt->first] = jt->second;
    }
  }

  if (!props.size()) throw runtime_error(
    "no rules matched for histogram \""+name+"\""
  );

  if (!props.count(kClass)) throw runtime_error(
    "cannot find class for histogram \""+name+"\""
  );
  string class_name(static_cast<const prop_Class*>(props[kClass])->name);

  if (!props.count(kBins)) throw runtime_error(
    "cannot find binning for histogram \""+name+"\""
  );
  if (static_cast<const prop_Bins*>(props[kBins])->isrange) {
    const prop_Bins_range* b
      = static_cast<const prop_Bins_range*>(props[kBins]);

    if (!class_name.compare("TH1F")) {
      h = new TH1F(name.c_str(),"",b->nbinsx,b->xlow,b->xup);
    } else if (!class_name.compare("TH1D")) {
      h = new TH1D(name.c_str(),"",b->nbinsx,b->xlow,b->xup);
    } else if (!class_name.compare("TH1I")) {
      h = new TH1I(name.c_str(),"",b->nbinsx,b->xlow,b->xup);
    } else if (!class_name.compare("TH1C")) {
      h = new TH1C(name.c_str(),"",b->nbinsx,b->xlow,b->xup);
    } else if (!class_name.compare("TH1S")) {
      h = new TH1S(name.c_str(),"",b->nbinsx,b->xlow,b->xup);
    } else throw runtime_error(
      "undefined histogram class "+class_name
    );

  } else {
    const prop_Bin_edges* b
      = static_cast<const prop_Bin_edges*>(props[kBins]);

    if (!class_name.compare("TH1F")) {
      h = new TH1F(name.c_str(),"",b->xbins.size()-1,&b->xbins[0]);
    } else if (!class_name.compare("TH1D")) {
      h = new TH1D(name.c_str(),"",b->xbins.size()-1,&b->xbins[0]);
    } else if (!class_name.compare("TH1I")) {
      h = new TH1I(name.c_str(),"",b->xbins.size()-1,&b->xbins[0]);
    } else if (!class_name.compare("TH1C")) {
      h = new TH1C(name.c_str(),"",b->xbins.size()-1,&b->xbins[0]);
    } else if (!class_name.compare("TH1S")) {
      h = new TH1S(name.c_str(),"",b->xbins.size()-1,&b->xbins[0]);
    } else throw runtime_error(
      "undefined histogram class "+class_name
    );
  }
  props.erase(kClass);
  props.erase(kBins);

  for (impl::prop_map::iterator it=props.begin(),end=props.end();it!=end;++it)
    it->second->apply(h);

  return h;
}

// Destructor *******************************************************

csshists::~csshists() { delete _impl; }
