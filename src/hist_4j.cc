#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
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
#include "hist_wt.hh"
#include "fj_jetdef.hh"
#include "int_range.hh"
#include "timed_counter.hh"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

template<typename T> inline T sq(const T x) { return x*x; }

// Constants ********************************************************
constexpr unsigned njets  = 4; // number of jets
constexpr unsigned n2jets = 6; // number of pairs

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, sj_files, wt_files, weights;
  string output_file, css_file, jet_alg;
  double pt_cut1, pt_cut4, eta_cut, dR_cut;
  int_range<Long64_t> ents {0,0};
  bool counter_newline, quiet;

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
      ("pt-cut1", po::value<double>(&pt_cut1)->default_value(100.,"100"),
       "first jet pT cut in GeV")
      ("pt-cut4", po::value<double>(&pt_cut4)->default_value(64.,"64"),
       "fourth jet pT cut in GeV")
      ("eta-cut", po::value<double>(&eta_cut)->default_value(2.8,"2.8"),
       "jet eta cut")
      ("dR-cut", po::value<double>(&dR_cut)->default_value(0.65,"0.65"),
       "jet minimum deltaR cut")
      ("style,s", po::value<string>(&css_file)
       ->default_value(CONFDIR"/4j.css","4j.css"),
       "CSS style file for histogram binning and formating")
      ("num-ent,n", po::value<int_range<Long64_t>>(&ents),
       "process only this many entries,\nnum or first:num")
      ("counter-newline", po::bool_switch(&counter_newline),
       "do not overwrite previous counter message")
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

  // Find number of entries to process
  if (ents.len>0) {
    const Long64_t need_ent = ents.first + ents.len;
    if (need_ent>tree->GetEntries()) {
      cerr << "Fewer entries in BH chain (" << tree->GetEntries()
         << ") then requested (" << need_ent << ')' << endl;
      exit(1);
    }
    if (sj_given) if (need_ent>sj_tree->GetEntries()) {
      cerr << "Fewer entries in SJ chain (" << sj_tree->GetEntries()
         << ") then requested (" << need_ent << ')' << endl;
      exit(1);
    }
    if (wt_given) if (need_ent>wt_tree->GetEntries()) {
      cerr << "Fewer entries in weights chain (" << wt_tree->GetEntries()
         << ") then requested (" << need_ent << ')' << endl;
      exit(1);
    }
  } else {
    ents.len = tree->GetEntries();
    if (sj_given) if (ents.len!=sj_tree->GetEntries()) {
      cerr << ents.len << " entries in BH chain, but "
           << sj_tree->GetEntries() << " entries in SJ chain" << endl;
      exit(1);
    }
    if (wt_given) if (ents.len!=wt_tree->GetEntries()) {
      cerr << ents.len << " entries in BH chain, but "
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
    jet_def.reset( fj_jetdef(jet_alg) );
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
  shared_ptr<csshists> hist_css( new csshists(css_file) );
  hist_wt::css.reset( hist_css.get() );
  cout << endl;

  // Open output file with histograms *******************************
  TFile* fout = new TFile(output_file.c_str(),"recreate");
  if (fout->IsZombie()) exit(1);
  else cout << "Output file: " << fout->GetName() << endl << endl;

  // Make directories ***********************************************
  for (auto& w : weight::all) {
    hist_wt::dirs[w.get()] = fout->mkdir((w->name+"_Jet"+jet_alg).c_str());
  }

  fout->cd();

  // Book histograms ************************************************
  TH1* h_N   = hist_css->mkhist("N");
  TH1* h_pid = hist_css->mkhist("pid");

  #define h_(name) h_##name(#name)

  /* NOTE:
   * excl = exactly the indicated number of jets, zero if no j in name
   * incl = that many or more jets
   *
   * y   = rapidity
   * eta = pseudo-rapidity
   */

  // Book Histograms
  hist_wt
    h_(jets_N_incl), h_(jets_N_excl),

    h_(jet1_pT), h_(jet2_pT), h_(jet3_pT), h_(jet4_pT),
    h_(4j_HT), h_(2j_HT),

    h_(jet1_y), h_(jet2_y), h_(jet3_y), h_(jet4_y),

    h_(4j_mass),

    h_(2j_mass_min),     h_(2j_mass_max),

    h_(2j_deltaphi_min), h_(2j_deltaphi_max),
    h_(2j_deltay_min),   h_(2j_deltay_max),

    h_(3j_deltaphi_min), h_(3j_deltaphi_max),
    h_(3j_deltay_min),   h_(3j_deltay_max)
  ;

  // Reading entries from the input TChain ***************************
  Long64_t num_selected = 0;
  Int_t prev_id = -1;
  cout << "Reading " << ents.len << " entries";
  if (ents.first>0) cout << " starting at " << ents.first;
  cout << endl;
  timed_counter counter(counter_newline);

  for (Long64_t ent = ents.first, ent_end = ents.end(); ent < ent_end; ++ent) {
    counter(ent);
    tree->GetEntry(ent);

    if (event.nparticle>BHMAXNP) {
      cerr << "More particles in the entry then BHMAXNP" << endl
           << "Increase array length to " << event.nparticle << endl;
      exit(1);
    }

    // Count number of events (not entries)
    if (prev_id!=event.eid) {
      h_N->Fill(0.5);
      ++num_selected;
    }
    prev_id = event.eid;

    // Fill histograms ***********************************
    for (Int_t i=0;i<event.nparticle;i++) h_pid->Fill(event.kf[i]);

    // Jet clustering *************************************
    vector<TLorentzVector> jets;
    jets.reserve(njets);
    if (sj_given) { // Read jets from SpartyJet ntuple
      jets = sj_alg->jetsByPt(pt_cut4,eta_cut);

    } else { // Clusted with FastJet on the fly
      vector<fastjet::PseudoJet> particles;
      particles.reserve(event.nparticle-1);

      for (Int_t i=0; i<event.nparticle; ++i) {
        particles.emplace_back(
          event.px[i],event.py[i],event.pz[i],event.E[i]
        );
      }

      // Cluster
      const vector<fastjet::PseudoJet> fj_jets =
        fastjet::ClusterSequence(particles, *jet_def).inclusive_jets(pt_cut4);

      // Apply pT cut & convert to TLorentzVector
      for (auto& j : fj_jets) {
        if (j.pt() < pt_cut4) continue;
        jets.emplace_back(j.px(),j.py(),j.pz(),j.E());
      }

      // Sort by pT in descending order
      std::sort( jets.begin(), jets.end(),
        [](const TLorentzVector& i, const TLorentzVector& j)
          { return i.Pt() > j.Pt(); }
      );
    }
    const size_t this_njets = jets.size(); // number of jets

    // ****************************************************

    // pT cut on the first jet
    if (this_njets) if (jets.front().Pt()<pt_cut1) continue;

    // Number of jets hists *******************************
    h_jets_N_excl.Fill(this_njets);
    for (unsigned i=0;i<this_njets;++i)
      if (this_njets >= i) h_jets_N_incl.Fill(i);

    if (this_njets < njets) continue;

    // Jets pT ********************************************
    static array<double,njets> pT, rap, phi;
    Double_t HT = 0.;
    for (size_t i=0;i<njets;++i) {
      HT += pT[i] = jets[i].Pt();
      rap[i] = jets[i].Rapidity();
      phi[i] = jets[i].Phi();
    }

    h_4j_HT.Fill(HT);
    h_jet1_pT.Fill(pT[0]);
    h_jet2_pT.Fill(pT[1]);
    h_jet3_pT.Fill(pT[2]);
    h_jet4_pT.Fill(pT[3]);

    h_jet1_y.Fill(rap[0]);
    h_jet2_y.Fill(rap[1]);
    h_jet3_y.Fill(rap[2]);
    h_jet4_y.Fill(rap[3]);

    // Sum of all jets ************************************
    const TLorentzVector all4 = jets[0] + jets[1] + jets[2] + jets[3];
    const Double_t m4 = all4.M();

    h_4j_mass.Fill(m4);

    // Jet pairs ******************************************
    static array<double,n2jets> dphi2_, dy2_;

    Double_t    m2_min =             (jets[0]+jets[1]).M();
    Double_t dphi2_min = dphi2_[0] = fabs(phi[0] - phi[1]);
    Double_t   dy2_min =   dy2_[0] = fabs(rap[0] - rap[1]);
    Double_t    m2_max =    m2_min;
    Double_t dphi2_max = dphi2_min;
    Double_t   dy2_max =   dy2_min;
    
    // To flatten traceless triangular matrix:
    // k = i*(i-1)/2 + j
    
    for (size_t i=2,k=1;i<njets;++i) {
      for (size_t j=0;j<i;++j,++k) {
        const Double_t    m2 =             (jets[i]+jets[j]).M();
        const Double_t dphi2 = dphi2_[k] = fabs(phi[i] - phi[j]);
        const Double_t   dy2 =   dy2_[k] = fabs(rap[i] - rap[j]);
        if (   m2 <    m2_min)    m2_min =    m2;
        if (   m2 >    m2_max)    m2_max =    m2;
        if (dphi2 < dphi2_min) dphi2_min = dphi2;
        if (dphi2 > dphi2_max) dphi2_max = dphi2;
        if (  dy2 <   dy2_min)   dy2_min =   dy2;
        if (  dy2 >   dy2_max)   dy2_max =   dy2;
      }
    }

    h_2j_mass_min    .Fill(   m2_min/m4 );
    h_2j_mass_max    .Fill(   m2_max/m4 );
    h_2j_deltaphi_min.Fill(dphi2_min);
    h_2j_deltaphi_max.Fill(dphi2_max);
    h_2j_deltay_min  .Fill(  dy2_min);
    h_2j_deltay_max  .Fill(  dy2_max);

    // Jet triplets ***************************************
    Double_t dphi3_min = dphi2_[0] + dphi2_[1];
    Double_t   dy3_min =   dy2_[0] +   dy2_[1];
    Double_t dphi3_max = dphi3_min;
    Double_t   dy3_max =   dy3_min;

    for (size_t i=2;i<n2jets;++i) {
      for (size_t j=0;j<i;++j) {
        if (i+j == n2jets-1) continue;

        const Double_t dphi3 = dphi2_[i] + dphi2_[j];
        const Double_t   dy3 =   dy2_[i] +   dy2_[j];

        if (dphi3 < dphi3_min) dphi3_min = dphi3;
        if (dphi3 > dphi3_max) dphi3_max = dphi3;
        if (  dy3 <   dy3_min)   dy3_min =   dy3;
        if (  dy3 >   dy3_max)   dy3_max =   dy3;
      }
    }

    h_3j_deltaphi_min.Fill(dphi3_min);
    h_3j_deltaphi_max.Fill(dphi3_max);
    h_3j_deltay_min  .Fill(  dy3_min);
    h_3j_deltay_max  .Fill(  dy3_max);

    // Sort by rapidity in ascending order ****************
    std::sort( jets.begin(), jets.end(),
      [](const TLorentzVector& i, const TLorentzVector& j)
        { return fabs(i.Rapidity()) < fabs(j.Rapidity()); }
    );

    h_2j_HT.Fill( jets[0].Pt() + jets[1].Pt() );

  } // END of event loop

  counter.prt(ents.end());
  cout << endl;
  cout << "Selected events: " << num_selected << endl;

  // Close files
  fout->Write();
  fout->Close();
  delete fout;
  delete tree;
  delete sj_tree;
  delete wt_tree;

  return 0;
}
