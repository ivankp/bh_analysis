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
#include <TKey.h>
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
#include "catstr.hh"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

template<typename T> inline T sq(T x) { return x*x; }

template<typename T> inline bool between(T a, T x, T b) {
  if (b < a) swap(a,b);
  return ( (a<=x) && (x<=b) );
}

// ******************************************************************
struct Jet {
private:
  inline Double_t _tau(Double_t Y) noexcept {
    // need rapidity here
    return sqrt( pT*pT + mass*mass )/( 2.*cosh(y - Y) );
  }
public:
  TLorentzVector p;
  Double_t mass, pT, y, tau;
  Jet(const TLorentzVector& _p, Double_t Y) noexcept
  : p(_p), mass(p.M()), pT(p.Pt()), y(p.Rapidity()), tau(_tau(Y))
  { }
  Jet(const fastjet::PseudoJet& _p, Double_t Y) noexcept
  : p(_p.px(),_p.py(),_p.pz(),_p.E()),
    mass(p.M()), pT(p.Pt()), y(p.Rapidity()), tau(_tau(Y))
  { }
};
// ******************************************************************

int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, sj_files, wt_files, weights;
  string output_file, css_file, jet_alg;
  size_t njets;
  double jet_pt_cut, jet_eta_cut;
  int_range<Long64_t> ents;
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
      ("njets,j", po::value<size_t>(&njets)->required(),
       "*minimum number of jets per ntuple entry")
      ("cluster,c", po::value<string>(&jet_alg)->default_value("AntiKt4"),
       "jet clustering algorithm: e.g. antikt4, kt6\n"
       "without --sj: select FastJet algorithm\n"
       "with --sj: read jets from SpartyJet ntuple")
      ("weight,w", po::value<vector<string>>(&weights),
       "weight branchs; if skipped:\n"
       "  without --wt: ntuple weight2 is used\n"
       "  with --wt: all weights from wt files")
      ("jet-pt-cut", po::value<double>(&jet_pt_cut)->default_value(30.,"30"),
       "jet pT cut in GeV")
      ("jet-eta-cut", po::value<double>(&jet_eta_cut)->default_value(4.4,"4.4"),
       "jet eta cut")
      ("style,s", po::value<string>(&css_file)
       ->default_value(CONFDIR"/Hjets.css","Hjets.css"),
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
  if (ents.len > 0) {
    const Long64_t need_ent = ents.end();
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
    cout << "Clustering with " << jet_def->description() << endl;
  }

  // Weights tree branches
  if (wt_given) {
    cout << endl;
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
    cout << endl;
  } else weight::add(tree,"weight2",false); // Use default ntuple weight2

  // Read CSS file with histogram properties
  cout << "Histogram CSS file: " << css_file << endl;
  shared_ptr<csshists> hist_css( new csshists(css_file) );
  hist_wt::css = hist_css;

  // Open output file with histograms *******************************
  TFile* fout = new TFile(output_file.c_str(),"recreate");
  if (fout->IsZombie()) exit(1);
  else cout << "Output file: " << fout->GetName() << endl << endl;

  // Make directories ***********************************************
  for (auto& w : weight::all) {
    hist_wt::dirs[w.get()] = fout->mkdir((w->name+"_Jet"+jet_alg).c_str());
  }

  fout->cd();

  const size_t njetsR = njets + 1;

  // Book histograms ************************************************
  #define h_(name) hist_wt h_##name(#name);

  #define h_jj(name) hist_wt h_##name( njets>1 ? #name : string() );

  #define h_jet_var(var) \
    vector<hist_wt> h_jet_##var; \
    h_jet_##var.reserve(njetsR); \
    for (size_t j=1; j<=njetsR; ++j) { \
      h_jet_##var.emplace_back(cat("jet",j,'_',#var)); \
    }

  /* NOTE:
   * excl = exactly the indicated number of jets, zero if no j in name
   * incl = that many or more jets
   *
   * VBF = vector boson fusion cut
   *
   * y   = rapidity
   * eta = pseudo-rapidity
   */

  // Book Histograms ************************************************
  TH1* h_N = hist_css->mkhist("N");

  h_(AA_pT); h_(AA_y); h_(AA_mass);

  h_(AA_cos_theta_star); h_(AA_pTt); h_(AA_dy);

  h_jet_var(pT);

  h_(jets_HT);

  vector<hist_wt> h_AAnj_pT; h_AAnj_pT.reserve(njetsR);
  for (size_t j=1; j<=njetsR; ++j) {
    h_AAnj_pT.emplace_back(cat('H',j,"j_pT"));
  }

  h_jet_var(y); h_jet_var(mass); h_jet_var(tau);

  h_(jets_tau_max); h_(jets_tau_sum);

  h_(jets_N_incl); h_(jets_N_excl);

  h_jj(jjpT_dy);
  hist_wt h_jjpT_dy_nj_excl ( njets>1 ? cat("jjpT_dy_",njets, "j_excl") : string() ),
          h_jjpT_dy_nRj_excl( njets>1 ? cat("jjpT_dy_",njetsR,"j_excl") : string() );

  h_jj(jjdy_dy);
  hist_wt h_jjdy_dy_nj_excl ( njets>1 ? cat("jjdy_dy_",njets, "j_excl") : string() ),
          h_jjdy_dy_nRj_excl( njets>1 ? cat("jjdy_dy_",njetsR,"j_excl") : string() );

  h_jj(jjpT_dphi); h_jj(jjdy_dphi);
  h_jj(AA_jjpT_dphi_VBF); h_jj(AA_jjdy_dphi_VBF);

  h_jj(jjpT_mass);    h_jj(jjdy_mass);
  h_jj(AAjjpT_mass);  h_jj(AAjjdy_mass);
  h_jj(AA_jjpT_dy);   h_jj(AA_jjdy_dy);
  h_jj(AA_jjpT_dphi); h_jj(AA_jjdy_dphi);

  h_jj(jjpT_N_jhj_excl); h_jj(jjdy_N_jhj_excl);

  const size_t ndy = 6;

  vector<vector<hist_wt>> h_jet_pT_jjpT(njetsR);
  if (njets>1) for (size_t j=0; j<njetsR; ++j) {
    h_jet_pT_jjpT[j].reserve(ndy);
    for (size_t i=0; i<ndy; ++i)
      h_jet_pT_jjpT[j].emplace_back(cat("jet",j+1,"_pT_jjpT_mindy",i+1));
  }

  vector<vector<hist_wt>> h_jet_pT_jjdy(njetsR);
  if (njets>1) for (size_t j=0; j<njetsR; ++j) {
    h_jet_pT_jjdy[j].reserve(ndy);
    for (size_t i=0; i<ndy; ++i)
      h_jet_pT_jjdy[j].emplace_back(cat("jet",j+1,"_pT_jjdy_mindy",i+1));
  }

  h_jj(jjpT_loose_VBF); h_jj(jjdy_loose_VBF);
  h_jj(jjpT_tight_VBF); h_jj(jjdy_tight_VBF);

  // Reading entries from the input TChain **************************
  Long64_t num_selected = 0, num_events = 0;
  Int_t prev_id = -1;
  cout << "Reading " << ents.len << " entries";
  if (ents.first>0) cout << " starting at " << ents.first;
  cout << endl;
  timed_counter counter(counter_newline);

  // variables
  Double_t jjpT_dy   = 0, jjdy_dy   = 0;
  Double_t jjpT_dphi = 0, jjdy_dphi = 0;

  // LOOP
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
        if (Ai1==-1) Ai1 = i;
        else {
          if (Ai2==-1) Ai2 = i;
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
      prev_id = event.eid;
      ++num_events;
      
      for (auto h : hist_wt::all) h->FillSumw2();
    }

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

    // diphoton 4-vector
    const TLorentzVector AA = A1 + A2;

    const Double_t AA_mass = AA.M();        // diphoton Mass

    if ( AA_mass < 115 || 135 < AA_mass ) continue;

    const Double_t AA_pT   = AA.Pt();       // diphoton Pt
    const Double_t AA_y    = AA.Rapidity(); // diphoton Rapidity
    const Double_t AA_phi  = AA.Phi();      // diphoton Phi

    // |costheta*| from 1307.1432
    const Double_t cts = abs(sinh(A1.Eta()-A2.Eta())) /
            						 sqrt(1.+sq(AA_pT/AA_mass)) *
            						 2.*A1.Pt()*A2.Pt()/sq(AA_mass);

    const Double_t AA_pTt =
      fabs(A1.Px()*A2.Py()-A2.Px()*A1.Py()) / ((A1-A2).Pt()*2);

    const Double_t AA_dy = fabs(A1.Rapidity() - A2.Rapidity());

    // Jet clustering *************************************
    vector<Jet> jets;
    jets.reserve(njetsR);
    size_t nj;

    if (sj_given) { // Read jets from SpartyJet ntuple
      const vector<TLorentzVector> sj_jets = sj_alg->jetsByPt(jet_pt_cut,jet_eta_cut);
      nj = sj_jets.size();

      // skip entry if not enough jets
      if (nj >= njets) {

        // add to jets container
        for (const auto& jet : sj_jets) jets.emplace_back(jet,AA_y);

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

      // Cluster and apply pT cut
      vector<fastjet::PseudoJet> fj_jets =
        fastjet::ClusterSequence(particles, *jet_def).inclusive_jets(jet_pt_cut);

      // Apply eta cut
      for (auto it=fj_jets.begin(); it!=fj_jets.end(); ) {
        if (abs(it->eta()) > jet_eta_cut) fj_jets.erase(it);
        else ++it;
      }

      nj = fj_jets.size();

      // skip entry if not enough jets
      if (nj >= njets) {

        // add to jets container
        for (const auto& jet : fj_jets) jets.emplace_back(jet,AA_y);

        // sort by pt
        sort(jets.begin(), jets.end(), [](const Jet& i, const Jet& j){
          return ( i.pT > j.pT );
        });

      }
    }
    // ****************************************************

    // Number of jets hists
    h_jets_N_excl.Fill(nj);
    for (size_t i=0; i<njetsR; i++) {
      if (nj >= i) h_jets_N_incl.Fill(i);
    }

    // Do no process entries with fewer then njets
    if (nj < njets) continue;

    // Increment selected entries
    ++num_selected;

    size_t jjdy_1 = 1, jjdy_2 = 0;

    // find Δy
    if (njets>1) {
      // jjpT
      jjpT_dy = abs(jets[0].y - jets[1].y);
      jjpT_dphi = fmod( abs(jets[0].p.Phi() - jets[1].p.Phi()), M_PI);

      // jjdy
      jjdy_dy = jjpT_dy;
      for (size_t i=2; i<nj; ++i) {
        for (size_t j=0; j<i; ++j) {
          const Double_t dy = abs(jets[i].y - jets[j].y);
          if (dy > jjdy_dy) {
            jjdy_dy = dy;
            jjdy_1 = i;
            jjdy_2 = j;
          }
        }
      }

      jjdy_dphi = fmod( abs(jets[jjdy_1].p.Phi() - jets[jjdy_2].p.Phi()), M_PI);
    }

    // Fill histograms ************************************
    h_AA_mass.Fill(AA_mass);
    h_AA_pT  .Fill(AA_pT);
    h_AA_y   .Fill(AA_y);

    h_AA_pTt .Fill(AA_pTt);
    h_AA_dy  .Fill(AA_dy);
    h_AA_cos_theta_star.Fill(cts);

    TLorentzVector AAnj = AA;
    for (size_t j=0, _nj=min(nj,njetsR); j<_nj; ++j) {
      h_jet_mass[j].Fill(jets[j].mass);
      h_jet_pT  [j].Fill(jets[j].pT  );
      h_jet_y   [j].Fill(jets[j].y   );
      h_jet_tau [j].Fill(jets[j].tau );

      h_AAnj_pT  [j].Fill( (AAnj += jets[j].p).Pt() );

      if (njets>1) {
        // dy by pT
        for (size_t i=0; i<ndy; ++i)
          if ( jjpT_dy > (i+1) )
            h_jet_pT_jjpT[j][i].Fill(jets[j].pT);

        // dy by dy
        for (size_t i=0; i<ndy; ++i)
          if ( jjdy_dy > (i+1) )
            h_jet_pT_jjdy[j][i].Fill(jets[j].pT);
      }
    }

    if (njets>1) {

      h_jjpT_dy.Fill(jjpT_dy);
      h_jjdy_dy.Fill(jjdy_dy);

      h_jjpT_dphi.Fill(jjpT_dphi);
      h_jjdy_dphi.Fill(jjdy_dphi);

      if (nj==njets) {
        h_jjpT_dy_nj_excl.Fill(jjpT_dy);
        h_jjdy_dy_nj_excl.Fill(jjdy_dy);
      } else if (nj==njetsR) {
        h_jjpT_dy_nRj_excl.Fill(jjpT_dy);
        h_jjdy_dy_nRj_excl.Fill(jjdy_dy);
      }

      TLorentzVector jj = jets[0].p + jets[1].p;
      Double_t jj_mass = jj.M();
      Double_t AA_jj_dphi = fmod( abs(AA_phi - jj.Phi()), M_PI);
      h_jjpT_mass   .Fill( jj_mass );
      h_AAjjpT_mass .Fill( (AA + jj).M() );
      h_AA_jjpT_dy  .Fill( abs(AA_y - jj.Rapidity()) );
      h_AA_jjpT_dphi.Fill( AA_jj_dphi );

      if (jjpT_dy>2.8 && jj_mass>400) {
        h_AA_jjpT_dphi_VBF.Fill( AA_jj_dphi );
        h_jjpT_loose_VBF .Fill( 0.5 );
        if (AA_jj_dphi>2.6) h_jjpT_tight_VBF.Fill( 0.5 );
      }

      if (jjdy_1 != 1) {
        jj = jets[jjdy_1].p + jets[jjdy_2].p;
        jj_mass   = jj.M();
        AA_jj_dphi = fmod( abs(AA_phi - jj.Phi()), M_PI);
      }
      h_jjdy_mass  .Fill( jj_mass );
      h_AAjjdy_mass .Fill( (AA + jj).M() );
      h_AA_jjdy_dy  .Fill( abs(AA_y - jj.Rapidity()) );
      h_AA_jjdy_dphi.Fill( AA_jj_dphi );

      if (jjdy_dy>2.8 && jj_mass>400) {
        h_AA_jjdy_dphi_VBF.Fill( AA_jj_dphi );
        h_jjdy_loose_VBF .Fill( 0.5 );
        if (AA_jj_dphi>2.6) h_jjdy_tight_VBF.Fill( 0.5 );
      }

      h_jjpT_N_jhj_excl.Fill( between(jets[0     ].y,AA_y,jets[1     ].y) ? nj : 0 );
      h_jjdy_N_jhj_excl.Fill( between(jets[jjdy_1].y,AA_y,jets[jjdy_1].y) ? nj : 0 );
    }

    Double_t jets_HT = 0, jets_tau_max = 0, jets_tau_sum = 0;
    for (const auto& jet : jets) {
      jets_HT += jet.pT;
      jets_tau_sum += jet.tau;
      if (jet.tau > jets_tau_max) jets_tau_max = jet.tau;
    }

    h_jets_HT.Fill(jets_HT);
    h_jets_tau_max.Fill(jets_tau_max);
    h_jets_tau_sum.Fill(jets_tau_sum);

  } // END of event loop ********************************************
  
  // finish correcting Sumw2
  for (auto h : hist_wt::all) {
    h->FillSumw2();
    h->AdoptSumw2();
  }

  counter.prt(ents.end());
  cout << endl;
  cout << "Selected entries: " << num_selected << endl;
  cout << "Processed events: " << num_events << endl;

  // Close files
  fout->Write();
  fout->Close();
  delete fout;
  delete tree;
  delete sj_tree;
  delete wt_tree;

  return 0;
}
