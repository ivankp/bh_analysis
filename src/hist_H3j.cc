#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <stdexcept>
#include <memory>
#include <algorithm>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "BHEvent.hh"
#include "SJClusterAlg.hh"
#include "weight.hh"
#include "timed_counter.hh"
#include "csshists.hh"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

template<typename T> inline T sq(const T x) { return x*x; }

// Histogram wrappers ***********************************************
struct hist {
  static unique_ptr<const csshists> css;
  virtual void Fill(Double_t x) noexcept =0;
};
unique_ptr<const csshists> hist::css;

class hist_wt: public hist {
  unordered_map<const weight*,TH1*> h;
public:
  hist_wt(const string& name) {
    TH1* hist = css->mkhist(name);
    // hist->Sumw2(false); // in ROOT6 true seems to be the default
    for (auto& wt : weight::all) {
      const weight *w = wt.get();
      dirs[w]->cd();
      h[w] = static_cast<TH1*>( hist->Clone() );
    }
    delete hist;
  }

  virtual void Fill(Double_t x) noexcept {
    for (auto& wt : weight::all)
      h[wt.get()]->Fill(x,wt->is_float ? wt->w.f : wt->w.d);
  }

  static unordered_map<const weight*,TDirectory*> dirs;
};
unordered_map<const weight*,TDirectory*> hist_wt::dirs;

// istream operators ************************************************
namespace std {
  template<class A, class B>
  istream& operator>> (istream& is, pair<A,B>& r) {
    string str;
    is >> str;
    const size_t sep = str.find(':');
    if (sep==string::npos) {
      r.first = 0;
      stringstream(str) >> r.second;
    } else {
      stringstream(str.substr(0,sep)) >> r.first;
      stringstream(str.substr(sep+1)) >> r.second;
    }
    return is;
  }
}

fastjet::JetDefinition* JetDef(string& str) {
  string::iterator it = --str.end();
  while (isdigit(*it)) --it;
  ++it;
  string name;
  transform(str.begin(), it, back_inserter(name), ::tolower);
  fastjet::JetAlgorithm alg;
  if (!name.compare("antikt")) alg = fastjet::antikt_algorithm;
  else if (!name.compare("kt")) alg = fastjet::kt_algorithm;
  else if (!name.compare("cambridge")) alg = fastjet::cambridge_algorithm;
  else throw runtime_error("Undefined jet clustering algorithm: "+name);
  return new fastjet::JetDefinition(
    alg,
    atof( string(it,str.end()).c_str() )/10.
  );
}

// ******************************************************************

struct Jet {
private:
  inline Double_t _tau(Double_t Eta) noexcept {
    return sqrt( pT*pT + mass*mass )/( 2.*cosh(eta - Eta) );
  }
public:
  Double_t mass, pT, eta, tau;
  Jet(const TLorentzVector& p, Double_t Eta) noexcept
  : mass(p.M()), pT(p.Pt()), eta(p.Eta()), tau(_tau(Eta))
  { }
  Jet(const fastjet::PseudoJet& p, Double_t Eta) noexcept
  : mass(p.m()), pT(sqrt(p.kt2())), eta(p.eta()), tau(_tau(Eta))
  { }
};

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, sj_files, wt_files, weights;
  string output_file, css_file, jet_alg;
  double pt_cut, eta_cut;
  pair<Long64_t,Long64_t> num_events;
  bool quiet;

  bool sj_given = false, wt_given = false;

  try {
    // General Options ------------------------------------
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "produce help message")
      ("bh", po::value< vector<string> >(&bh_files)->required(),
       "*add input BlackHat root file")
      ("sj", po::value< vector<string> >(&sj_files),
       "add input SpartyJet root file")
      ("wt", po::value< vector<string> >(&wt_files),
       "add input weights root file")
      ("output,o", po::value<string>(&output_file)->required(),
       "*output root file with histograms")
      ("cluster,c", po::value<string>(&jet_alg)->default_value("AntiKt4"),
       "jet clustering algorithm: e.g. antikt4, kt6\n"
       "without --sj: select FastJet algorithm\n"
       "with --sj: read jets from SpartyJet ntuple")
      ("weight,w", po::value<vector<string>>(&weights),
       "weight branchs; if skipped:\n"
       "  without --wt: ntuple weight is used\n"
       "  with --wt: all weights from wt files")
      ("jet-pt-cut", po::value<double>(&pt_cut)->default_value(30.),
       "jet pT cut in GeV")
      ("jet-eta-cut", po::value<double>(&eta_cut)->default_value(4.4,"4.4"),
       "jet eta cut in GeV")
      ("style,s", po::value<string>(&css_file)
       ->default_value(CONFDIR"/Hj.css","Hj.css"),
       "CSS style file for histogram binning and formating")
      ("num-events,n", po::value<pair<Long64_t,Long64_t>>(&num_events),
       "process only this many events,\nnum or first:num")
      ("quiet,q", po::bool_switch(&quiet),
       "Do not print exception messages")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (argc == 1 || vm.count("help")) {
      cout << desc << endl;
      return 0;
    }
    po::notify(vm);
    if (vm.count("sj")) sj_given = true;
    if (vm.count("wt")) wt_given = true;
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Setup input files **********************************************
  TChain*    tree = new TChain("t3");
  TChain* sj_tree = (sj_given ? new TChain("SpartyJet_Tree") : nullptr);
  TChain* wt_tree = (wt_given ? new TChain("weights") : nullptr);

  // Add trees from all the files to the TChains
  cout << "BH files:" << endl;
  for (auto& f : bh_files) {
    cout << "  " << f << endl;
    if (!tree->AddFile(f.c_str(),-1) ) exit(1);
  }
  if (sj_given) {
    cout << "SJ files:" << endl;
    for (auto& f : sj_files) {
      cout << "  " << f << endl;
      if (!sj_tree->AddFile(f.c_str(),-1) ) exit(1);
    }
  }
  if (wt_given) {
    cout << "Weight files:" << endl;
    for (auto& f : wt_files) {
      cout << "  " << f << endl;
      if (!wt_tree->AddFile(f.c_str(),-1) ) exit(1);
    }
  }
  cout << endl;

  // Find number of events to process
  if (num_events.second>0) {
    const Long64_t need_events = num_events.first + num_events.second;
    if (need_events>tree->GetEntries()) {
      cerr << "Fewer entries in BH chain (" << tree->GetEntries()
         << ") then requested (" << need_events << ')' << endl;
      exit(1);
    }
    if (sj_given) if (need_events>sj_tree->GetEntries()) {
      cerr << "Fewer entries in SJ chain (" << sj_tree->GetEntries()
         << ") then requested (" << need_events << ')' << endl;
      exit(1);
    }
    if (wt_given) if (need_events>wt_tree->GetEntries()) {
      cerr << "Fewer entries in weights chain (" << wt_tree->GetEntries()
         << ") then requested (" << need_events << ')' << endl;
      exit(1);
    }
  } else {
    num_events.second = tree->GetEntries();
    if (sj_given) if (num_events.second!=sj_tree->GetEntries()) {
      cerr << num_events.second << " entries in BH chain, but "
           << sj_tree->GetEntries() << " entries in SJ chain" << endl;
      exit(1);
    }
    if (wt_given) if (num_events.second!=wt_tree->GetEntries()) {
      cerr << num_events.second << " entries in BH chain, but "
           << wt_tree->GetEntries() << " entries in weights chain" << endl;
      exit(1);
    }
  }

  // Friend BlackHat tree with SpartyJet and Weight trees
  if (sj_given) tree->AddFriend(sj_tree,"SJ");
  if (wt_given) tree->AddFriend(wt_tree,"weights");

  // BlackHat tree branches
  BHEvent event;
  event.SetTree(tree, BHEvent::kinematics);

  // Jet Clustering Algorithm
  unique_ptr<fastjet::JetDefinition> jet_def;
  unique_ptr<SJClusterAlg> sj_alg;

  if (sj_given) {
    sj_alg.reset( new SJClusterAlg(tree,jet_alg) );
  } else {
    jet_def.reset( JetDef(jet_alg) );
    cout << "Clustering with " << jet_def->description() << endl << endl;
  }

  // Weights tree branches
  if (wt_given) {
    if (weights.size()) {
      cout << "Selected weights:" << endl;
      for (auto& w : weights) {
        cout << w << endl;
        weight::add(tree,w);
      }
    } else {
      cout << "Using all weights:" << endl;
      const TObjArray *br = wt_tree->GetListOfBranches();
      for (Int_t i=0,n=br->GetEntries();i<n;++i) {
        auto w = br->At(i)->GetName();
        cout << w << endl;
        weight::add(tree,w);
      }
    }
  } else weight::add(tree,"weight",false); // Use default ntuple weight
  cout << endl;

  // Read CSS file with histogram properties
  cout << "Histogram CSS file: " << css_file << endl;
  hist::css.reset( new csshists(css_file) );
  cout << endl;

  // Open output file with histograms *******************************
  TFile* fout = new TFile(output_file.c_str(),"recreate");
  if (fout->IsZombie()) exit(1);
  else cout << "Output file: " << fout->GetName() << endl << endl;

  // Make directories ***********************************************
  for (auto& w : weight::all) {
    hist_wt::dirs[w.get()] = fout->mkdir(w->name.c_str());
  }

  fout->cd();

  // Book histograms ************************************************
  TH1* h_N   = hist::css->mkhist("N");
  TH1* h_pid = hist::css->mkhist("pid");

  #define h_(name) h_##name(#name)

  /* NOTE:
   * excl = exactly the indicated number of jets, zero if no j in name
   * incl = that many or more jets
   * VBF = vector boson fusion cut
   */

  // Book Histograms
  hist_wt
    h_(H_mass),

    h_(jets_N_incl), h_(jets_N_excl), h_(jets_N_incl_pT50), h_(jets_N_excl_pT50),

    h_(H_pT_3j), h_(H_pT_3j_excl), h_(H_eta_3j), h_(H_eta_3j_excl),
    h_(jet3_mass), h_(jet3_pT), h_(jet3_eta), h_(jet3_tau),

    h_(H_pT_2j), h_(H_pT_2j_excl), h_(H_eta_2j), h_(H_eta_2j_excl),
    h_(jet2_mass), h_(jet2_pT), h_(jet2_eta), h_(jet2_tau),

    h_(H_pT_1j), h_(H_pT_1j_excl), h_(H_eta_1j), h_(H_eta_1j_excl),
    h_(jet1_mass), h_(jet1_pT), h_(jet1_eta), h_(jet1_tau),

    h_(H_pT_0j), h_(H_pT_0j_excl), h_(H_eta_0j), h_(H_eta_0j_excl),

    h_(jets_HT), h_(jets_tau_max), h_(jets_tau_sum)
  ;

  // Reading events from the input TChain ***************************
  Long64_t numOK = 0;
  Int_t prev_id = -1;
  cout << "Reading " << num_events.second << " entries";
  if (num_events.first>0) cout << " starting at " << num_events.first << endl;
  else cout << endl;
  num_events.second += num_events.first;
  timed_counter counter;

  for (Long64_t ent = num_events.first; ent < num_events.second; ++ent) {
    counter(ent);
    tree->GetEntry(ent);

    if (event.nparticle>BHMAXNP) {
      cerr << "More particles in the event then BHMAXNP" << endl
           << "Increase array length to " << event.nparticle << endl;
      exit(1);
    }

    // Find Higgs
    Int_t hi = 0; // Higgs index
    while (hi<event.nparticle) {
      if (event.kf[hi]==25) break;
      else ++hi;
    }
    if (hi==event.nparticle) {
      cerr << "No Higgs in event " << ent << endl;
      continue;
    }
    ++numOK;

    // Count number of events (not entries)
    if (prev_id!=event.eid) h_N->Fill(0.5);
    prev_id = event.eid;

    // Higgs 4-vector
    const TLorentzVector higgs(event.px[hi],event.py[hi],event.pz[hi],event.E[hi]);

    const Double_t H_mass = higgs.M();   // Higgs Mass
    const Double_t H_pT   = higgs.Pt();  // Higgs Pt
    const Double_t H_eta  = higgs.Eta(); // Higgs Pseudo-rapidity

    // Fill histograms ***********************************
    for (Int_t i=0;i<event.nparticle;i++) h_pid->Fill(event.kf[i]);

    h_H_mass  .Fill(H_mass);
    h_H_pT_0j .Fill(H_pT);
    h_H_eta_0j.Fill(H_eta);

    // Jet clustering *************************************
    vector<Jet> jets;
    if (sj_given) { // Read jets from SpartyJet ntuple
      for (auto& jet : sj_alg->jetsByPt(pt_cut,eta_cut))
        jets.emplace_back(jet,H_eta);

    } else { // Clusted with FastJet on the fly
      vector<fastjet::PseudoJet> particles;
      particles.reserve(event.nparticle-1);

      for (Int_t i=0; i<event.nparticle; ++i) {
        if (i==hi) continue;
        particles.emplace_back(
          event.px[i],event.py[i],event.pz[i],event.E[i]
        );
      }

      // Cluster, sort jets by pT, and apply pT cut
      const vector<fastjet::PseudoJet> fj_jets = sorted_by_pt(
        fastjet::ClusterSequence(particles, *jet_def).inclusive_jets(pt_cut)
      );

      // Apply eta cut
      jets.reserve(fj_jets.size());
      for (auto& j : fj_jets) {
        if (abs(j.eta()) < eta_cut) jets.emplace_back(j,H_eta);
      }
    }
    const size_t njets = jets.size(); // number of jets

    // ****************************************************

    int njets50 = 0;
    for (auto& j : jets) {
      if (j.pT>=50.) ++njets50;
      else break;
    }

    // Number of jets hists
    h_jets_N_excl.Fill(njets);
    h_jets_N_excl_pT50.Fill(njets50);
    for (unsigned char i=0;i<4;i++) {
      if (njets >= i) {
        h_jets_N_incl.Fill(i);
        if (njets50 >= i) h_jets_N_incl_pT50.Fill(i);
      }
    }

    if (njets==0) { // njets == 0;

      h_H_pT_0j_excl.Fill(H_pT);
      h_H_pT_0j_excl.Fill(H_eta);

    }
    else { // njets > 0;

      h_H_pT_1j  .Fill(H_pT);
      h_H_eta_1j .Fill(H_eta);

      h_jet1_mass.Fill(jets[0].mass);
      h_jet1_pT  .Fill(jets[0].pT);
      h_jet1_eta .Fill(jets[0].eta);
      h_jet1_tau .Fill(jets[0].tau);

      Double_t jets_HT = 0, jets_tau_max = 0, jets_tau_sum = 0;

      for (auto& jet : jets) {
        jets_HT += jet.pT;
        jets_tau_sum += jet.tau;
        if (jet.tau > jets_tau_max) jets_tau_max = jet.tau;
      }
      h_jets_HT.Fill(jets_HT);
      h_jets_tau_max.Fill(jets_tau_max);
      h_jets_tau_sum.Fill(jets_tau_sum);

      if (njets==1) { // njets == 1;

        h_H_pT_1j_excl .Fill(H_pT);
        h_H_eta_1j_excl.Fill(H_eta);

      }
      else { // njets > 1;

        h_H_pT_2j  .Fill(H_pT);
        h_H_eta_2j .Fill(H_eta);

        h_jet2_mass.Fill(jets[1].mass);
        h_jet2_pT  .Fill(jets[1].pT);
        h_jet2_eta .Fill(jets[1].eta);
        h_jet2_tau .Fill(jets[1].tau);

        if (njets==2) { // njets == 1;

          h_H_pT_2j_excl .Fill(H_pT);
          h_H_eta_2j_excl.Fill(H_eta);

        }
        else { // njets > 2;

          h_H_pT_3j  .Fill(H_pT);
          h_H_eta_3j .Fill(H_eta);

          h_jet3_mass.Fill(jets[2].mass);
          h_jet3_pT  .Fill(jets[2].pT);
          h_jet3_eta .Fill(jets[2].eta);
          h_jet3_tau .Fill(jets[2].tau);

          if (njets==3) { // njets == 1;

            h_H_pT_3j_excl .Fill(H_pT);
            h_H_eta_3j_excl.Fill(H_eta);

          }

        } // END njets > 2;

      } // END njets > 1;

    } // END njets > 0;

  } // END of event loop

  counter.prt(num_events.second);
  cout << endl;
  cout << "Successfully processed events: " << numOK << endl;

  // Close files
  fout->Write();
  fout->Close();
  delete fout;
  delete tree;
  delete sj_tree;
  delete wt_tree;

  return 0;
}
