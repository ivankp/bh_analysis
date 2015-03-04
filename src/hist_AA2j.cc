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
#include "hist_wt.hh"
#include "fj_jetdef.hh"
#include "int_range.hh"
#include "timed_counter.hh"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

template<typename T> inline T sq(const T x) { return x*x; }

// ******************************************************************

struct Jet {
private:
  inline Double_t _tau(Double_t Y) noexcept {
    // need rapidity here
    return sqrt( pT*pT + mass*mass )/( 2.*cosh(y - Y) );
  }
public:
  TLorentzVector *p;
  Double_t mass, pT, y, tau;
  Jet(const TLorentzVector& p, Double_t Y, bool keep=false) noexcept
  : p(keep ? new TLorentzVector(p) : nullptr),
    mass(p.M()), pT(p.Pt()), y(p.Rapidity()), tau(_tau(Y))
  { }
  Jet(const fastjet::PseudoJet& p, Double_t Y, bool keep=false) noexcept
  : p(keep ? new TLorentzVector(p.px(),p.py(),p.pz(),p.E()) : nullptr),
    mass(p.m()), pT(sqrt(p.kt2())), y(p.rapidity()), tau(_tau(Y))
  { }
  ~Jet() { delete p; }
};

// ******************************************************************
int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, sj_files, wt_files, weights;
  string output_file, css_file, jet_alg;
  double jet_pt_cut, jet_eta_cut;
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
      ("jet-pt-cut", po::value<double>(&jet_pt_cut)->default_value(30.,"30"),
       "jet pT cut in GeV")
      ("jet-eta-cut", po::value<double>(&jet_eta_cut)->default_value(4.4,"4.4"),
       "jet eta cut")
      ("style,s", po::value<string>(&css_file)
       ->default_value(CONFDIR"/H2j.css","H2j.css"),
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
   * VBF = vector boson fusion cut
   *
   * y   = rapidity
   * eta = pseudo-rapidity
   */

  // Book Histograms
  hist_wt
    h_(H_mass),

    h_(jets_N_incl), h_(jets_N_excl), h_(jets_N_incl_pT50), h_(jets_N_excl_pT50),

    h_(H_pT_2j), h_(H_pT_2j_excl),
    h_(H_pT_1j), h_(H_pT_1j_excl),
    h_(H_pT_0j), h_(H_pT_0j_excl),


    h_(H_y_2j), h_(H_y_2j_excl),
    h_(H_y_1j), h_(H_y_1j_excl),
    h_(H_y_0j), h_(H_y_0j_excl),

    h_(H2j_pT), h_(H2j_pT_excl),
    h_(H1j_pT), h_(H1j_pT_excl),

    h_(jet1_mass), h_(jet2_mass),
    h_(jet1_pT),   h_(jet2_pT),
    h_(jet1_y),    h_(jet2_y),
    h_(jet1_tau),  h_(jet2_tau),

    h_(jets_HT), h_(jets_tau_max), h_(jets_tau_sum),

    h_(H2j_mass),

    h_(H_2j_deltaphi), h_(H_2j_deltaphi_excl),
    h_(H_2j_deltay), h_(H_2j_deltay_excl),

    h_(2j_mass),
    h_(j_j_deltaphi), h_(j_j_deltaphi_excl), h_(j_j_deltaphi_VBF),
    h_(j_j_deltay),

    h_(loose), h_(tight)
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

    // Find two photons
    Int_t Ai1 = -1, Ai2 = -1; // Higgs indices
    for (Int_t i=0; i<event.nparticle; ++i) {
      if (event.kf[i]==22) {
        if (Ai1==-1) Ai1 = event.kf[i];
        else {
          if (Ai2==-1) Ai2 = event.kf[i];
          break;
        }
      }
    }
    if (Ai2==-1) {
      cerr << "Event " << ent << " doesn't have 2 photons" << endl;
      continue;
    }

    // Count number of events (not entries)
    if (prev_id!=event.eid) {
      h_N->Fill(0.5);
      ++num_selected;
    }
    prev_id = event.eid;

    const TLorentzVector A1(event.px[Ai1], event.py[Ai1],
                            event.pz[Ai1], event.E [Ai1]);
    const TLorentzVector A2(event.px[Ai2], event.py[Ai2],
                            event.pz[Ai2], event.E [Ai2]);

    // Photon cuts
    if ( A1.Pt() > A2.Pt() ) {
      if ( A1.Et() < 43.75 ) continue;
      if ( A2.Et() < 31.35 ) continue;
    } else {
      if ( A2.Et() < 43.75 ) continue;
      if ( A1.Et() < 31.35 ) continue;
    }
    if ( A1.Eta() > 2.37 ) continue;
    if ( A2.Eta() > 2.37 ) continue;


    // Higgs 4-vector
    const TLorentzVector higgs = A1 + A2;

    const Double_t H_mass = higgs.M();        // Higgs Mass

    if ( H_mass < 115 || 135 < H_mass ) continue;

    const Double_t H_pT   = higgs.Pt();       // Higgs Pt
    const Double_t H_y    = higgs.Rapidity(); // Higgs Rapidity

    // Fill histograms ***********************************
    for (Int_t i=0;i<event.nparticle;i++) h_pid->Fill(event.kf[i]);

    h_H_mass .Fill(H_mass);
    h_H_pT_0j.Fill(H_pT);
    h_H_y_0j .Fill(H_y);

    // Jet clustering *************************************
    vector<Jet> jets;
    if (sj_given) { // Read jets from SpartyJet ntuple
      const vector<TLorentzVector> sj_jets = sj_alg->jetsByPt(jet_pt_cut,jet_eta_cut);
      jets.reserve(sj_jets.size());
      for (auto& jet : sj_jets) {
        jets.emplace_back(jet,H_y,jets.size()<2);
      }

    } else { // Clusted with FastJet on the fly
      vector<fastjet::PseudoJet> particles;
      particles.reserve(event.nparticle-1);

      for (Int_t i=0; i<event.nparticle; ++i) {
        if (i==Ai1) continue;
        if (i==Ai2) continue;
        particles.emplace_back(
          event.px[i],event.py[i],event.pz[i],event.E[i]
        );
      }

      // Cluster, sort jets by pT, and apply pT cut
      const vector<fastjet::PseudoJet> fj_jets = sorted_by_pt(
        fastjet::ClusterSequence(particles, *jet_def).inclusive_jets(jet_pt_cut)
      );

      // Apply eta cut
      jets.reserve(fj_jets.size());
      for (auto& jet : fj_jets) {
        if (abs(jet.eta()) < jet_eta_cut)
          jets.emplace_back(jet,H_y,jets.size()<2);
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

    if (njets==0) { // njets == 0; --------------------------------=0

      h_H_pT_0j_excl.Fill(H_pT);
      h_H_y_0j_excl .Fill(H_y);

    }
    else { // njets > 0; ------------------------------------------>0

      h_H_pT_1j  .Fill(H_pT);
      h_H_y_1j   .Fill(H_y);

      h_jet1_mass.Fill(jets[0].mass);
      h_jet1_pT  .Fill(jets[0].pT);
      h_jet1_y   .Fill(jets[0].y);
      h_jet1_tau .Fill(jets[0].tau);

      const Double_t H1j_pT = (higgs+(*jets[0].p)).Pt();

      h_H1j_pT   .Fill(H1j_pT);

      Double_t jets_HT = 0, jets_tau_max = 0, jets_tau_sum = 0;

      for (auto& jet : jets) {
        jets_HT += jet.pT;
        jets_tau_sum += jet.tau;
        if (jet.tau > jets_tau_max) jets_tau_max = jet.tau;
      }
      h_jets_HT     .Fill(jets_HT);
      h_jets_tau_max.Fill(jets_tau_max);
      h_jets_tau_sum.Fill(jets_tau_sum);

      if (njets==1) { // njets == 1; ------------------------------=1

        h_H_pT_1j_excl.Fill(H_pT);
        h_H_y_1j_excl .Fill(H_y);
        h_H1j_pT_excl .Fill(H1j_pT);

      }
      else { // njets > 1; ---------------------------------------->1

        h_H_pT_2j  .Fill(H_pT);
        h_H_y_2j   .Fill(H_y);

        h_jet2_mass.Fill(jets[1].mass);
        h_jet2_pT  .Fill(jets[1].pT);
        h_jet2_y   .Fill(jets[1].y);
        h_jet2_tau .Fill(jets[1].tau);

        const TLorentzVector jj = (*jets[0].p)+(*jets[1].p);
        const TLorentzVector H2j = higgs+jj;

        const Double_t H2j_mass      = H2j.M();
        const Double_t H2j_pT        = H2j.Pt();
        const Double_t H_2j_deltaphi = higgs.Phi() - jj.Phi();
        const Double_t H_2j_deltay   = H_y - jj.Rapidity();

        const Double_t jj_mass       = jj.M();
        const Double_t j_j_deltaphi  = jets[0].p->Phi() - jets[1].p->Phi();
        const Double_t j_j_deltay    = jets[0].y - jets[1].y;

        h_H2j_mass     .Fill(H2j_mass);
        h_H2j_pT       .Fill(H2j_pT);
        h_H_2j_deltaphi.Fill(H_2j_deltaphi);
        h_H_2j_deltay  .Fill(H_2j_deltay);
        h_2j_mass      .Fill(jj_mass);

        h_j_j_deltaphi .Fill(j_j_deltaphi);
        h_j_j_deltay   .Fill(j_j_deltay);

        if (j_j_deltay>2.8) { // VBF cuts
          if (jj_mass>400) {
            h_j_j_deltaphi_VBF.Fill(j_j_deltaphi);
            h_loose.Fill(0.5);
            if (H_2j_deltaphi>2.6) h_tight.Fill(0.5);
          }
        }

        if (njets==2) { // njets == 2; ----------------------------=2

          h_H_pT_2j_excl      .Fill(H_pT);
          h_H_y_2j_excl       .Fill(H_y);
          h_H2j_pT_excl       .Fill(H2j_pT);
          h_H_2j_deltaphi_excl.Fill(H_2j_deltaphi);
          h_H_2j_deltay_excl  .Fill(H_2j_deltay);
          h_j_j_deltaphi_excl .Fill(j_j_deltaphi);

        }

      } // END njets > 1;

    } // END njets > 0;

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
