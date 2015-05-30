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
#include "real_range.hh"
#include "timed_counter.hh"
#include "catstr.hh"
#include "senum.hh"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

template<typename T> inline T sq(T x) noexcept { return x*x; }

template<typename T> inline bool between(T a, T x, T b) noexcept {
  if (b < a) swap(a,b);
  return ( (a<=x) && (x<=b) );
}

inline Double_t dphi(Double_t phi1, Double_t phi2) noexcept {
  return fmod( abs(phi1 - phi2), M_PI);
}

inline Double_t fphi2(const TLorentzVector& b, const TLorentzVector& f) noexcept {
  const Double_t bx=b.Px(), by=b.Py(), fx=f.Px(), fy=f.Py();

  const Double_t phi2 = acos( ( bx*fx+by*fy ) /
    ( sqrt(bx*bx+by*by) * sqrt(fx*fx+fy*fy) )
  );
  if ( ( bx*fy - by*fx ) < 0.) return -phi2;
  else return phi2;
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
  Double_t mass, pT, y, phi, tau;
  Jet(const TLorentzVector& _p, Double_t Y) noexcept
  : p(_p), mass(p.M()), pT(p.Pt()), y(p.Rapidity()), phi(p.Phi()), tau(_tau(Y))
  { }
  Jet(const fastjet::PseudoJet& _p, Double_t Y) noexcept
  : p(_p.px(),_p.py(),_p.pz(),_p.E()),
    mass(p.M()), pT(p.Pt()), y(p.Rapidity()), phi(p.Phi()), tau(_tau(Y))
  { }
};
// ******************************************************************

senum(VBF,(none)(hardest)(any))

int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, sj_files, wt_files, weights;
  string output_file, css_file, jet_alg;
  size_t njets;
  double jet_pt_cut, jet_eta_cut;
  real_range<Double_t> AA_mass_cut;
  int_range<Long64_t> ents;
  VBF::type VBFcut;
  bool AAntuple, counter_newline, quiet, strict;
  Long64_t cache_size;

  bool sj_given = false, wt_given = false;
  bool apply_AA_mass_cut = false;

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

      ("weight,w", po::value<vector<string>>(&weights),
       "weight branchs; if skipped:\n"
       "  without --wt: ntuple weight2 is used\n"
       "  with --wt: all weights from wt files")

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
       ->default_value(CONFDIR"/Hjets.css","Hjets.css"),
       "CSS style file for histogram binning and formating")
      ("cache", po::value<Long64_t>(&cache_size)->default_value(50),
       "cache size in Mb")

      ("num-ent,n", po::value<int_range<Long64_t>>(&ents),
       "process only this many entries,\nnum or first:num")
      ("counter-newline", po::bool_switch(&counter_newline),
       "do not overwrite previous counter message")
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
    if (vm.count("sj")) sj_given = true;
    if (vm.count("wt")) wt_given = true;
    if (AAntuple && vm.count("AA-mass-cut")) apply_AA_mass_cut = true;
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Setup input files **********************************************
  TChain* const bh_tree = new TChain("t3");
  TChain* const sj_tree = (sj_given ? new TChain("SpartyJet_Tree") : nullptr);
  TChain* const wt_tree = (wt_given ? new TChain("weights") : nullptr);

  // Add trees from all the files to the TChains
  cout << "BH files:" << endl;
  for (auto& f : bh_files) {
    cout << "  " << f << endl;
    if (!bh_tree->AddFile(f.c_str(),-1) ) exit(1);
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
    if (need_ent>bh_tree->GetEntries()) {
      cerr << "Fewer entries in BH chain (" << bh_tree->GetEntries()
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
    ents.len = bh_tree->GetEntries();
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
  //if (sj_given) tree->AddFriend(sj_tree,"SJ");
  //if (wt_given) tree->AddFriend(wt_tree,"weights");

  // BlackHat tree branches
  BHEvent event;
  event.SetTree(bh_tree, BHEvent::kinematics);

  // Jet Clustering Algorithm
  unique_ptr<fastjet::JetDefinition> jet_def;
  unique_ptr<SJClusterAlg> sj_alg;

  if (sj_given) {
    sj_alg.reset( new SJClusterAlg(sj_tree,jet_alg) );
  } else {
    jet_def.reset( fj_jetdef(jet_alg) );
    cout << "Clustering with " << jet_def->description() << endl;
  }

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
    if (wt_given && sj_given) cache_size /= 3;
    else if (wt_given || sj_given) cache_size /= 2;
  
    bh_tree->SetCacheSize(cache_size);
    auto *brs = bh_tree->GetListOfBranches();
    for ( Int_t i=0, n=brs->GetEntries(); i<n; ++i ) {
      TBranch *br = static_cast<TBranch*>( brs->At(i) );
      if (!br->TestBit(kDoNotProcess)) {
        // cout << br->GetName() << " added to cache" << endl;
        bh_tree->AddBranchToCache(br,kTRUE);
      }
      bh_tree->StopCacheLearningPhase();
    }

    if (wt_given) {
      wt_tree->SetCacheSize(cache_size);
      brs = wt_tree->GetListOfBranches();
      for ( Int_t i=0, n=brs->GetEntries(); i<n; ++i ) {
        TBranch *br = static_cast<TBranch*>( brs->At(i) );
        if (!br->TestBit(kDoNotProcess)) {
          // cout << br->GetName() << " added to cache" << endl;
          wt_tree->AddBranchToCache(br,kTRUE);
        }
      }
      wt_tree->StopCacheLearningPhase();
    }

    if (sj_given) {
      sj_tree->SetCacheSize(cache_size);
      brs = sj_tree->GetListOfBranches();
      for ( Int_t i=0, n=brs->GetEntries(); i<n; ++i ) {
        TBranch *br = static_cast<TBranch*>( brs->At(i) );
        if (!br->TestBit(kDoNotProcess)) {
          // cout << br->GetName() << " added to cache" << endl;
          sj_tree->AddBranchToCache(br,kTRUE);
        }
      }
      sj_tree->StopCacheLearningPhase();
    }
  }

  // Read CSS file with histogram properties ************************
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

  #define h_opt(name,opt) hist_wt h_##name( opt ? #name : string() );

  #define h_jj(name) h_opt(name, njets>1)
  
  #define h_jj_excl(name) \
    hist_wt h_##name##_nj_excl ( njets>1 ? cat(#name "_",njets, "j_excl") : string() ), \
            h_##name##_nRj_excl( njets>1 ? cat(#name "_",njetsR,"j_excl") : string() );
    

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

  h_(H_pT) h_(H_y) h_(H_mass)

  h_opt(AA_cos_theta_star,           AAntuple)
  h_opt(AA_cos_theta_star_AApT80,    AAntuple)
  h_opt(AA_cos_theta_star_80AApT200, AAntuple)
  h_opt(AA_cos_theta_star_200AApT,   AAntuple)

  h_opt(AA_pTt, AAntuple)
  h_opt(AA_dy,  AAntuple)

  h_jet_var(pT)

  h_(jets_HT)

  vector<hist_wt> h_Hnj_pT; h_Hnj_pT.reserve(njetsR);
  for (size_t j=1; j<=njetsR; ++j) {
    h_Hnj_pT.emplace_back(cat('H',j,"j_pT"));
  }

  h_jet_var(y) h_jet_var(mass) h_jet_var(tau)

  h_(jets_tau_max) h_(jets_tau_sum)

  h_(jets_N_incl) h_(jets_N_excl)

  h_jj(jjpT_dy)   h_jj_excl(jjpT_dy)
  h_jj(jjfb_dy)   h_jj_excl(jjfb_dy)
  h_jj(jjpT_dphi) h_jj_excl(jjpT_dphi)
  h_jj(jjfb_dphi) h_jj_excl(jjfb_dphi)

  h_jj(jjpT_mass) h_jj(jjfb_mass);

  h_jj(H_jjpT_dy)   h_jj_excl(H_jjpT_dy)   h_jj(H_jjpT_dy_avgyjj)
  h_jj(H_jjfb_dy)   h_jj_excl(H_jjfb_dy)   h_jj(H_jjfb_dy_avgyjj)
  h_jj(H_jjpT_dphi) h_jj_excl(H_jjpT_dphi)
  h_jj(H_jjfb_dphi) h_jj_excl(H_jjfb_dphi)

  h_jj(H_jjfb_phi2)
  h_jj(HjjpT_mass) h_jj(Hjjfb_mass);

  h_jj(jjpT_N_jhj_incl) h_jj(jjfb_N_jhj_incl)
  h_jj(jjpT_N_jhj_excl) h_jj(jjfb_N_jhj_excl)

  const size_t ndy = 6;

  vector<vector<hist_wt>> h_jet_pT_jjpT(njetsR);
  if (njets>1) for (size_t j=0; j<njetsR; ++j) {
    h_jet_pT_jjpT[j].reserve(ndy);
    for (size_t i=0; i<ndy; ++i)
      h_jet_pT_jjpT[j].emplace_back(cat("jet",j+1,"_pT_jjpT_mindy",i+1));
  }

  vector<vector<hist_wt>> h_jet_pT_jjfb(njetsR);
  if (njets>1) for (size_t j=0; j<njetsR; ++j) {
    h_jet_pT_jjfb[j].reserve(ndy);
    for (size_t i=0; i<ndy; ++i)
      h_jet_pT_jjfb[j].emplace_back(cat("jet",j+1,"_pT_jjfb_mindy",i+1));
  }

  h_jj(jj_loose) h_jj(jjpT_tight) h_jj(jjfb_tight)

  h_opt(jjpT_j_dy_veto, njets>2)
  h_opt(jjfb_j_dy_veto, njets>2)

  // Reading entries from the input TChain **************************
  Long64_t num_selected = 0, num_events = 0;
  Int_t prev_id = -1;
  cout << "Reading " << ents.len << " entries";
  if (ents.first>0) cout << " starting at " << ents.first;
  cout << endl;
  timed_counter counter(counter_newline);
  
  // Variables ******************************************************
  size_t jb = 0, jf = 0;

  TLorentzVector A1, A2, higgs, jjpT, jjfb;
  
  Double_t H_mass, H_pT, H_y, H_phi;

  Double_t jjpT_dy=0, jjpT_dphi=0, jjpT_ycenter=0, jjpT_mass=0,
           jjfb_dy=0, jjfb_dphi=0, jjfb_ycenter=0, jjfb_mass=0;
  
  // LOOP ***********************************************************
  for (Long64_t ent = ents.first, ent_end = ents.end(); ent < ent_end; ++ent) {
    counter(ent);
    bh_tree->GetEntry(ent);

    if (event.nparticle>BHMAXNP) {
      cerr << "\033[31mMore particles in entry then BHMAXNP = "
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
        if (strict) exit(1); else continue;
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

    H_y = higgs.Rapidity(); // Higgs Rapidity

    // Jet clustering *************************************
    vector<Jet> jets;
    jets.reserve(njetsR);
    size_t nj;

    if (sj_given) { // Read jets from SpartyJet ntuple
      sj_tree->GetEntry(ent);
      const vector<TLorentzVector> sj_jets = sj_alg->jetsByPt(jet_pt_cut,jet_eta_cut);
      nj = sj_jets.size();

      // skip entry if not enough jets
      if (nj >= njets) {

        // add to jets container
        for (const auto& jet : sj_jets) jets.emplace_back(jet,H_y);

      }

    } else { // Clusted with FastJet on the fly
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
        for (const auto& jet : fj_jets) jets.emplace_back(jet,H_y);

        // sort by pt
        sort(jets.begin(), jets.end(), [](const Jet& i, const Jet& j){
          return ( i.pT > j.pT );
        });

      }
    }
    // ****************************************************

    H_pT  = higgs.Pt();  // Higgs Pt
    H_phi = higgs.Phi(); // Higgs Phi
    
    if (njets>1) {
      jjpT = jets[0].p + jets[1].p;
      jjpT_dy = abs(jets[0].y - jets[1].y);
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

    // Do no process entries with fewer then njets
    if (nj < njets) continue;

    // Increment selected entries
    ++num_selected;
    
    if (njets>1) {
      jjpT_dphi = dphi(jets[0].phi, jets[1].phi);
      jjpT_ycenter = (jets[0].y + jets[1].y)/2;

      // find most forward and backward jets
      jb = jf = 0;
      for (size_t j=1; j<nj; ++j) {
        if (jets[j].y < jets[jb].y) jb = j;
        if (jets[j].y > jets[jf].y) jf = j;
      }
      jjfb_dy   = jets[jf].y - jets[jb].y;
      jjfb_dphi = dphi(jets[jb].phi, jets[jf].phi);
      jjfb_ycenter = (jets[jf].y + jets[jb].y)/2;
    
      jjfb = jets[jb].p + jets[jf].p;
      jjfb_mass = jjfb.M();
    }

    // ****************************************************
    // Fill histograms ************************************

    h_H_mass.Fill(H_mass);
    h_H_pT  .Fill(H_pT);
    h_H_y   .Fill(H_y);

    TLorentzVector Hnj = higgs;
    for (size_t j=0, _nj=min(nj,njetsR); j<_nj; ++j) {
      h_jet_mass[j].Fill(jets[j].mass);
      h_jet_pT  [j].Fill(jets[j].pT  );
      h_jet_y   [j].Fill(jets[j].y   );
      h_jet_tau [j].Fill(jets[j].tau );

      h_Hnj_pT  [j].Fill( (Hnj += jets[j].p).Pt() );

      if (njets>1) {
        // dy by pT
        for (size_t i=0; i<ndy; ++i)
          if ( jjpT_dy > (i+1) )
            h_jet_pT_jjpT[j][i].Fill(jets[j].pT);

        // dy by dy
        for (size_t i=0; i<ndy; ++i)
          if ( jjfb_dy > (i+1) )
            h_jet_pT_jjfb[j][i].Fill(jets[j].pT);
      }
    }

    if (njets>1) {

      // Leading pT jets are tagging
      // ----------------------------------------
      Double_t H_jj_dy   = abs(H_y - jjpT.Rapidity());
      Double_t H_jj_dphi = dphi(jjpT.Phi(),H_phi);
      
      h_jjpT_dy    .Fill(jjpT_dy);
      h_jjpT_dphi  .Fill(jjpT_dphi);
      h_H_jjpT_dy  .Fill(H_jj_dy);
      h_H_jjpT_dphi.Fill(H_jj_dphi);
      
      if (nj==njets) {
        h_jjpT_dy_nj_excl     .Fill(jjpT_dy);
        h_jjpT_dphi_nj_excl   .Fill(jjpT_dphi);
        
        h_H_jjpT_dy_nj_excl   .Fill(H_jj_dy);
        h_H_jjpT_dphi_nj_excl .Fill(H_jj_dphi);
      } else {
        h_jjpT_dy_nRj_excl    .Fill(jjpT_dy);
        h_jjpT_dphi_nRj_excl  .Fill(jjpT_dphi);
        
        h_H_jjpT_dy_nRj_excl  .Fill(H_jj_dy);
        h_H_jjpT_dphi_nRj_excl.Fill(H_jj_dphi);
      }
      
      h_jjpT_mass .Fill( jjpT_mass );
      h_HjjpT_mass.Fill( (higgs + jjpT).M() );

      h_H_jjpT_dy_avgyjj.Fill( abs(H_y - jjpT_ycenter) );
      
      h_jj_loose.Fill( 0.5 );
      if (H_jj_dphi > 2.6) h_jjpT_tight.Fill( 0.5 );

      size_t NjHj = ( between(jets[0].y,H_y,jets[1].y) ? nj : 0 );
      h_jjpT_N_jhj_excl.Fill(NjHj);
      h_jjpT_N_jhj_incl.FillIncl(NjHj);
      
      if (njets>2) {
        // Third photon veto:
        // ydists is the smallest distance between the centre of the
        // tagging jets and any possible further jet
        Double_t y_dists = 100.;
        for (size_t j=2; j<nj; ++j) {
          Double_t y_distt = abs(jets[j].y - jjpT_ycenter);
          if (y_distt < y_dists) y_dists = y_distt;
        }
        h_jjpT_j_dy_veto.FillIncl(y_dists);
      }
      // ----------------------------------------

      // Forward-backward jets are tagging
      // ----------------------------------------
      H_jj_dy   = abs(H_y - jjfb.Rapidity());
      H_jj_dphi = dphi(jjfb.Phi(),H_phi);
      
      h_jjfb_dy    .Fill(jjfb_dy);
      h_jjfb_dphi  .Fill(jjfb_dphi);
      h_H_jjfb_dy  .Fill(H_jj_dy);
      h_H_jjfb_dphi.Fill(H_jj_dphi);
      
      if (nj==njets) {
        h_jjfb_dy_nj_excl     .Fill(jjfb_dy);
        h_jjfb_dphi_nj_excl   .Fill(jjfb_dphi);
        
        h_H_jjfb_dy_nj_excl   .Fill(H_jj_dy);
        h_H_jjfb_dphi_nj_excl .Fill(H_jj_dphi);
      } else {
        h_jjfb_dy_nRj_excl    .Fill(jjfb_dy);
        h_jjfb_dphi_nRj_excl  .Fill(jjfb_dphi);
        
        h_H_jjfb_dy_nRj_excl  .Fill(H_jj_dy);
        h_H_jjfb_dphi_nRj_excl.Fill(H_jj_dphi);
      }
      
      h_jjfb_mass .Fill( jjfb_mass );
      h_Hjjfb_mass.Fill( (higgs + jjfb).M() );

      h_H_jjfb_dy_avgyjj.Fill( abs(H_y - jjfb_ycenter) );
      
      if (H_jj_dphi > 2.6) h_jjfb_tight.Fill( 0.5 );

      NjHj = ( between(jets[jb].y,H_y,jets[jb].y) ? nj : 0 );
      h_jjfb_N_jhj_excl.Fill(NjHj);
      h_jjfb_N_jhj_incl.FillIncl(NjHj);
      
      if (njets>2) {
        // Third photon veto:
        // ydists is the smallest distance between the centre of the
        // tagging jets and any possible further jet
        Double_t y_dists = 100.;
        for (size_t j=0; j<nj; ++j) {
          if (j==jb || j==jf) continue;
          Double_t y_distt = abs(jets[j].y - jjfb_ycenter);
          if (y_distt < y_dists) y_dists = y_distt;
        }
        h_jjfb_j_dy_veto.FillIncl(y_dists);
      }
      
      // Calculation of phi_2 from arXiv:1001.3822  
      // phi_2 = azimuthal angle between the vector sum of jets 
      // forward and jets backward of the Higgs boson

      TLorentzVector vsumf(0.,0.,0.,0.);
      TLorentzVector vsumb(0.,0.,0.,0.);
  
      Double_t f_nonzero = false, b_nonzero = false;
      for (const auto& j : jets) {
        if (j.y > H_y) { vsumf += j.p; f_nonzero = true; }
        else           { vsumb += j.p; b_nonzero = true; }
      }

      // Calculate phi_2
      Double_t phi2;
      if (f_nonzero && b_nonzero) {
        phi2 = fphi2(vsumb,vsumf);
      } else if (!f_nonzero) {
        vsumb -= jets[jf].p;
        phi2 = fphi2(vsumb,jets[jf].p);
      } else { 
        vsumf -= jets[jb].p;
        phi2 = fphi2(jets[jb].p,vsumf);
      }
  
      h_H_jjfb_phi2.Fill(phi2);

    } // END if (njets>1)

    // Diphoton histograms
    if (AAntuple) {
    
      // |costheta*| from 1307.1432
      const Double_t cts =
        abs(sinh(A1.Eta()-A2.Eta())) / sqrt(1.+sq(H_pT/H_mass)) *
        A1.Pt()*A2.Pt()*2 / sq(H_mass);

      h_AA_cos_theta_star.Fill(cts);
      if (H_pT<80) h_AA_cos_theta_star_AApT80.Fill(cts);
      else if (H_pT<200) h_AA_cos_theta_star_80AApT200.Fill(cts);
      else h_AA_cos_theta_star_200AApT.Fill(cts);

      h_AA_pTt.Fill(
        abs(A1.Px()*A2.Py()-A2.Px()*A1.Py()) / ((A1-A2).Pt()*2)
      );

      h_AA_dy.Fill( abs(A1.Rapidity() - A2.Rapidity()) );
      
    } // END AA

    // Jets tau
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
  delete bh_tree;
  delete sj_tree;
  delete wt_tree;

  return 0;
}
