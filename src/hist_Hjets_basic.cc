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
#include <TKey.h>
#include <TH1.h>
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "BHEvent.hh"
#include "weight.hh"
#include "hist_wt.hh"
#include "fj_jetdef.hh"
#include "int_range.hh"
#include "real_range.hh"
#include "timed_counter2.hh"
#include "catstr.hh"
#include "senum.hh"
#include "nest.hh"

#include "dphi.hh"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

template <typename T> [[ gnu::const ]]
inline T sq(T x) noexcept { return x*x; }
template <typename T, typename... TT> [[ gnu::const ]]
inline T sq(T x, TT... xx) noexcept { return sq(x)+sq(xx...); }


template <typename T> [[ gnu::const ]]
inline bool between(T a, T x, T b) noexcept {
  if (__builtin_expect(b<a,0)) swap(a,b);
  return ( (a<=x) && (x<=b) );
}

// ******************************************************************
struct Jet {
  TLorentzVector p;
  Double_t mass, pT, y, eta, phi;

  template <typename P>
  Jet(const P& _p) noexcept
  : p(_p[0],_p[1],_p[2],_p[3]),
    mass(p.M()), pT(p.Pt()), y(p.Rapidity()), eta(p.Eta()), phi(p.Phi())
  { }
};
// ******************************************************************

senum(VBF,(none)(hardest)(any))

int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, wt_files, weights;
  string output_file, css_file, jet_alg, tree_name;
  size_t njets;
  double jet_pt_cut, jet_eta_cut;
  real_range<Double_t> AA_mass_cut;
  int_range<Long64_t> ents;
  VBF::type VBFcut;
  bool AAntuple, quiet, strict;
  Long64_t cache_size;

  bool wt_given = false;
  bool apply_AA_mass_cut = false;

  try {
    // General Options ------------------------------------
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "produce help message")

      ("bh", po::value< vector<string> >(&bh_files)->required(),
       "*add input BlackHat root file")
      ("wt", po::value< vector<string> >(&wt_files),
       "add input weights root file")
      ("output,o", po::value<string>(&output_file)->required(),
       "*output root file with histograms")

      ("weight,w", po::value<vector<string>>(&weights),
       "weight branchs; if skipped:\n"
       "  without --wt: ntuple weight2 is used\n"
       "  with --wt: all weights from wt files")
      ("tree-name", po::value(&tree_name)->default_value("t3"),
       "change ntuple TTree name")

      ("njets,j", po::value<size_t>(&njets)->required(),
       "*minimum number of jets per ntuple entry")
      ("jet-alg,c", po::value<string>(&jet_alg)->default_value("AntiKt4"),
       "jet clustering algorithm: e.g. antikt4, kt6\n"
       "without --sj: select FastJet algorithm\n"
       "with --sj: read jets from SpartyJet ntuple")
      ("jet-pt-cut", po::value<double>(&jet_pt_cut)->default_value(30.,"30"),
       "jet pT cut in GeV")
      ("jet-eta-cut", po::value<double>(&jet_eta_cut)->default_value(4.4,"4.4"),
       "jet eta cut")

      ("VBF", po::value<VBF::type>(&VBFcut)->default_value(VBF::none),
       "apply Vector Bososn Fusion cuts\n(jj_dy>2.8 && jj_mass>400)\n"
       "to hardest or any jets")

      ("AA", po::bool_switch(&AAntuple),
       "make Higgs from diphoton\nand produce AA histograms")
      ("AA-mass-cut", po::value<real_range<Double_t>>(&AA_mass_cut),
       "apply a mass cut to the diphoton,\ne.g. 115:135")

      ("style,s", po::value<string>(&css_file)
       ->default_value(CONFDIR"/Hjets_basic.css","Hjets_basic.css"),
       "CSS style file for histogram binning and formating")
      ("cache", po::value<Long64_t>(&cache_size)->default_value(50),
       "cache size in Mb")

      ("num-ent,n", po::value<int_range<Long64_t>>(&ents),
       "process only this many entries,\nnum or first:num")
      ("quiet,q", po::bool_switch(&quiet),
       "Do not print exception messages")
      ("strict,z", po::bool_switch(&strict),
       "Exit on no particle of interest")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (argc == 1 || vm.count("help")) {
      cout << desc << endl;
      return 0;
    }
    po::notify(vm);
    if (vm.count("wt")) wt_given = true;
    if (AAntuple && vm.count("AA-mass-cut")) apply_AA_mass_cut = true;
  } catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    return 1;
  }
  // END OPTIONS ****************************************************

  // Setup input files **********************************************
  TChain* const bh_tree = new TChain(tree_name.c_str());
  TChain* const wt_tree = (wt_given ? new TChain("weights") : nullptr);

  // Add trees from all the files to the TChains
  cout << "BH files:" << endl;
  for (auto& f : bh_files) {
    cout << "  " << f << endl;
    if (!bh_tree->AddFile(f.c_str(),-1) ) return 1;
  }
  if (wt_given) {
    cout << "Weight files:" << endl;
    for (auto& f : wt_files) {
      cout << "  " << f << endl;
      if (!wt_tree->AddFile(f.c_str(),-1) ) return 1;
    }
  }
  cout << endl;

  // Find number of entries to process
  if (ents.len > 0) {
    const Long64_t need_ent = ents.end();
    if (need_ent>bh_tree->GetEntries()) {
      cerr << "Fewer entries in BH chain (" << bh_tree->GetEntries()
         << ") then requested (" << need_ent << ')' << endl;
      return 1;
    }
    if (wt_given) if (need_ent>wt_tree->GetEntries()) {
      cerr << "Fewer entries in weights chain (" << wt_tree->GetEntries()
         << ") then requested (" << need_ent << ')' << endl;
      return 1;
    }
  } else {
    ents.len = bh_tree->GetEntries();
    if (wt_given) if (ents.len!=wt_tree->GetEntries()) {
      cerr << ents.len << " entries in BH chain, but "
           << wt_tree->GetEntries() << " entries in weights chain" << endl;
      return 1;
    }
  }

  // BlackHat tree branches
  BHEvent event;
  event.SetTree(bh_tree, BHEvent::kinematics);

  // Jet Clustering Algorithm
  unique_ptr<fastjet::JetDefinition> jet_def;

  jet_def.reset( fj_jetdef(jet_alg) );
  cout << "Clustering with " << jet_def->description() << endl;

  // Weights tree branches ******************************************
  if (wt_given) {
    cout << endl;
    if (weights.size()) {
      wt_tree->SetBranchStatus("*",0);
      cout << "Selected weights:" << endl;
      for (auto& w : weights) {
        cout << w << endl;
        weight::add(wt_tree,w);
      }
    } else {
      cout << "Using all weights:" << endl;
      const TObjArray *br = wt_tree->GetListOfBranches();
      for (Int_t i=0,n=br->GetEntries();i<n;++i) {
        auto w = br->At(i)->GetName();
        cout << w << endl;
        weight::add(wt_tree,w);
      }
    }
    cout << endl;
  } else { // Use default ntuple weight2
    weight::add(bh_tree,"weight2",false);
  }

  // Use TTreeCache *************************************************
  if (cache_size > 0) {
    cache_size *= 1024*1024;
    if (wt_given) cache_size /= 2;

    bh_tree->SetCacheSize(cache_size);
    auto *brs = bh_tree->GetListOfBranches();
    for ( Int_t i=0, n=brs->GetEntries(); i<n; ++i ) {
      TBranch *br = static_cast<TBranch*>( brs->At(i) );
      if (!br->TestBit(kDoNotProcess)) {
        #if ROOT_VERSION_CODE >= ROOT_VERSION(6,04,00)
          if (bh_tree->AddBranchToCache(br,kTRUE) < 0) {
            cerr << "In BH Tree: Could not cache branch "
                 << br->GetName() << endl;
            return 1;
          }
        #else
          bh_tree->AddBranchToCache(br,kTRUE);
        #endif
      }
    }
    bh_tree->StopCacheLearningPhase();

    if (wt_given) {
      wt_tree->SetCacheSize(cache_size);
      brs = wt_tree->GetListOfBranches();
      for ( Int_t i=0, n=brs->GetEntries(); i<n; ++i ) {
        TBranch *br = static_cast<TBranch*>( brs->At(i) );
        if (!br->TestBit(kDoNotProcess)) {
          #if ROOT_VERSION_CODE >= ROOT_VERSION(6,04,00)
            if (wt_tree->AddBranchToCache(br,kTRUE) < 0) {
              cerr << "In weight Tree: Could not cache branch "
                   << br->GetName() << endl;
              return 1;
            }
          #else
            wt_tree->AddBranchToCache(br,kTRUE);
          #endif
        }
      }
      wt_tree->StopCacheLearningPhase();
    }
  }

  // Read CSS file with histogram properties ************************
  cout << "Histogram CSS file: " << css_file << endl;
  shared_ptr<csshists> hist_css( new csshists(css_file) );
  hist_wt::css = hist_css;

  // Open output file with histograms *******************************
  TFile* fout = new TFile(output_file.c_str(),"recreate");
  if (fout->IsZombie()) return 1;
  else cout << "Output file: " << fout->GetName() << endl << endl;

  // Make directories ***********************************************
  for (auto& w : weight::all) {
    hist_wt::dirs[w.get()] = fout->mkdir((w->name+"_Jet"+jet_alg).c_str());
  }

  fout->cd();

  const size_t njetsR = njets + 1; // number of jest in Real part

  // Book histograms ************************************************
  #define h_(name) hist_wt h_##name(#name);

  #define h_opt(name,opt) hist_wt h_##name( opt ? #name : string() );

  #define h_jj(name) h_opt(name, njets>1)

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

  h_(H_pT)
  h_(H_y)
  h_(H_eta)
  h_(H_phi)
  h_(H_mass)

  h_opt(AA_cos_theta_star, AAntuple)
  h_opt(AA_dy, AAntuple)

  h_(jets_N_incl) h_(jets_N_excl)
  h_jet_var(pT)
  h_jet_var(y)
  h_jet_var(eta)
  h_jet_var(phi)
  h_jet_var(mass)

  h_jj(jjpT_dpT)  h_jj(jjfb_dpT)
  h_jj(jjpT_dy)   h_jj(jjfb_dy)
  h_jj(jjpT_deta) h_jj(jjfb_deta)
  h_jj(jjpT_dphi) h_jj(jjfb_dphi)
  h_jj(jjpT_mass) h_jj(jjfb_mass)

  // Reading entries from the input TChain **************************
  Long64_t num_selected = 0, num_events = 0;
  Int_t prev_id = -1;
  cout << "Reading " << ents.len << " entries";
  if (ents.first>0) cout << " starting at " << ents.first;
  cout << endl;

  // Variables ******************************************************
  size_t jb = 0, jf = 0;

  TLorentzVector A1, A2, higgs, jjpT, jjfb;

  Double_t H_pT, H_y, H_eta, H_phi, H_mass;

  Double_t jjpT_dpT=0, jjpT_dy=0, jjpT_deta=0, jjpT_dphi=0, jjpT_mass=0,
           jjfb_dpT=0, jjfb_dy=0, jjfb_deta=0, jjfb_dphi=0, jjfb_mass=0;

  // LOOP ***********************************************************
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(ents.first,ents.end()); ent.ok(); ++ent) {
    bh_tree->GetEntry(ent);

    if (event.nparticle>BHMAXNP) {
      cerr << "\033[31mMore particles in entry than BHMAXNP = "
           << BHMAXNP << "\033[0m\n"
           << "Increase array length to " << event.nparticle << endl;
      break;
    }

    // Find Higgs or AA
    Int_t Hi = 0, Ai1 = -1, Ai2 = -1; // Higgs or phton indices

    if (AAntuple) { // AA

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
        if (!quiet)
          cerr << "\033[31mEvent " << ent
               << " doesn't have 2 photons\033[0m" << endl;
        if (strict) break; else continue;
      }

    } else { // Higgs

      while (Hi<event.nparticle) {
        if (event.kf[Hi]==25) break;
        else ++Hi;
      }
      if (Hi==event.nparticle) {
        if (!quiet)
          cerr << "\033[31mNo Higgs in event " << ent <<"\033[0m"<< endl;
        if (strict) return 1; else continue;
      }

    }

    // Count number of events (not entries)
    if (prev_id!=event.eid) {
      h_N->Fill(0.5);
      prev_id = event.eid;
      ++num_events;

      for (auto h : hist_wt::all) h->FillSumw2();
    }

    if (AAntuple) {
      A1 = TLorentzVector(event.px[Ai1], event.py[Ai1],
                          event.pz[Ai1], event.E [Ai1]);
      A2 = TLorentzVector(event.px[Ai2], event.py[Ai2],
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

      higgs = A1 + A2;

    } else {
      higgs = TLorentzVector(event.px[Hi], event.py[Hi],
                             event.pz[Hi], event.E [Hi]);
    }

    H_mass = higgs.M(); // Higgs Mass

    if (apply_AA_mass_cut)
      if ( H_mass < AA_mass_cut.min || AA_mass_cut.max < H_mass ) continue;


    // Jet clustering *************************************
    vector<Jet> jets;
    jets.reserve(njetsR);
    size_t nj;

    { // Clustering with FastJet on the fly
      vector<fastjet::PseudoJet> particles;
      particles.reserve(event.nparticle-1);

      for (Int_t i=0; i<event.nparticle; ++i) {
        if (AAntuple) {
          if (i==Ai1) continue;
          if (i==Ai2) continue;
        } else if (i==Hi) continue;
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
        for (const auto& jet : fj_jets) jets.emplace_back(jet);

        // sort by pt
        sort(jets.begin(), jets.end(), [](const Jet& i, const Jet& j){
          return ( i.pT > j.pT );
        });

      }
    }
    // ****************************************************

    if (njets>1) {
      // Leading (pT) dijet
      // These variables are computed here
      // because they are needed for the VBF cuts
      // The rest of the lead dijet variables are computed after the cuts
      jjpT      = jets[0].p + jets[1].p; // Four-momentum
      jjpT_dy   = abs(jets[0].y - jets[1].y);
      jjpT_dphi = dphi(jets[0].phi, jets[1].phi);
      jjpT_mass = jjpT.M();

      // Apply VBF cuts
      // Always refering to the two hardest jets
      bool passedVBF = true;
      if (VBFcut!=VBF::none) {
        passedVBF = (jjpT_dy>2.8 && jjpT_mass>400);

        if (nj>2 && VBFcut==VBF::any && !passedVBF) {
          for (size_t i=2; i<nj && !passedVBF; ++i)
            for (size_t j=0; j<i && !passedVBF; ++j)
              if ( abs(jets[i].y - jets[j].y)>2.8 &&
                      (jets[i].p + jets[j].p).M()>400 )
                passedVBF = true;
        }
      }
      if (!passedVBF) continue;
    }

    if (wt_given) wt_tree->GetEntry(ent);

    // Number of jets hists
    h_jets_N_excl.Fill(nj);
    h_jets_N_incl.FillIncl(nj);

    // Do not process entries with fewer then njets
    if (nj < njets) continue;

    // Increment selected entries
    ++num_selected;

    if (njets>1) {
      // Leading (pT) dijet
      jjpT_dpT  = jets[0].pT - jets[1].pT;
      jjpT_deta = abs(jets[0].eta - jets[1].eta);
      jjpT_dphi = dphi(jets[0].phi, jets[1].phi);

      // Forward-backward (fb) dijet
      jb = jf = 0;
      for (size_t j=1; j<nj; ++j) {
        if (jets[j].y < jets[jb].y) jb = j;
        if (jets[j].y > jets[jf].y) jf = j;
      }
      jjfb = jets[jb].p + jets[jf].p;
      jjfb_dpT  = jets[jf].pT - jets[jb].pT;
      jjfb_dy   = abs(jets[jf].y - jets[jb].y);
      jjfb_deta = abs(jets[jf].eta - jets[jb].eta);
      jjfb_dphi = dphi(jets[jb].phi, jets[jf].phi);
      jjfb_mass = jjfb.M();
    }

    // Higgs variables
    H_pT  = higgs.Pt();
    H_y   = higgs.Rapidity();
    H_eta = higgs.Eta();
    H_phi = higgs.Phi();

    // ****************************************************
    // Fill histograms ************************************

    // Higgs
    h_H_pT  .Fill(H_pT);
    h_H_y   .Fill(H_y);
    h_H_eta .Fill(H_eta);
    h_H_phi .Fill(H_phi);
    h_H_mass.Fill(H_mass);

    if (AAntuple) { // Diphoton
      // |cos(theta*)| from 1307.1432
      const Double_t cts =
        abs(sinh(A1.Eta()-A2.Eta())) / sqrt(1.+sq(H_pT/H_mass)) *
        A1.Pt()*A2.Pt()*2 / sq(H_mass);

      h_AA_cos_theta_star.Fill(cts);
      h_AA_dy.Fill( abs(A1.Rapidity() - A2.Rapidity()) );
    }

    // Individual jets
    for (size_t j=0, _nj=min(nj,njetsR); j<_nj; ++j) {
      h_jet_pT  [j].Fill(jets[j].pT  );
      h_jet_y   [j].Fill(jets[j].y   );
      h_jet_eta [j].Fill(jets[j].eta );
      h_jet_phi [j].Fill(jets[j].phi );
      h_jet_mass[j].Fill(jets[j].mass);
    }

    // Dijet
    if (njets>1) {
      // Leading pT jets
      h_jjpT_dpT .Fill(jjpT_dpT);
      h_jjpT_dy  .Fill(jjpT_dy);
      h_jjpT_deta.Fill(jjpT_deta);
      h_jjpT_dphi.Fill(jjpT_dphi);
      h_jjpT_mass.Fill(jjpT_mass);

      // Forward-backward jets
      h_jjfb_dpT .Fill(jjfb_dpT);
      h_jjfb_dy  .Fill(jjfb_dy);
      h_jjfb_deta.Fill(jjfb_deta);
      h_jjfb_dphi.Fill(jjfb_dphi);
      h_jjfb_mass.Fill(jjfb_mass);
    }

  } // END of event loop ********************************************

  // finish correcting Sumw2
  for (auto h : hist_wt::all) {
    h->FillSumw2();
    h->AdoptSumw2();
  }

  cout << "Selected entries: " << num_selected << endl;
  cout << "Processed events: " << num_events << endl;

  // Close files
  fout->Write();
  fout->Close();
  delete fout;
  delete bh_tree;
  delete wt_tree;

  return 0;
}
