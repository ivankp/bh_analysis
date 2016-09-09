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
#include "weight.hh"
#include "hist_wt.hh"
#include "fj_jetdef.hh"
#include "int_range.hh"
#include "real_range.hh"
#include "timed_counter.hh"
#include "catstr.hh"

#include "dphi.hh"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

template<typename T> inline T sq(T x) noexcept { return x*x; }

template<typename T> inline bool between(T a, T x, T b) noexcept {
  if (b < a) swap(a,b);
  return ( (a<=x) && (x<=b) );
}

Double_t fphi2(const TLorentzVector& b, const TLorentzVector& f) noexcept {
  const Double_t bx=b.Px(), by=b.Py(), fx=f.Px(), fy=f.Py();

  const Double_t phi2 = acos( ( bx*fx+by*fy ) /
    ( sqrt(bx*bx+by*by) * sqrt(fx*fx+fy*fy) )
  );
  if ( ( bx*fy - by*fx ) < 0.) return -phi2;
  else return phi2;
}

// ******************************************************************
struct Particle {
  TLorentzVector p;
  Double_t mass, pT, y, eta, phi, tau;
  Particle(const TLorentzVector& _p) noexcept
  : p(_p), mass(p.M()), pT(p.Pt()), y(p.Rapidity()), eta(p.Eta()), phi(p.Phi())
  { }
  Particle(const fastjet::PseudoJet& _p) noexcept
  : p(_p.px(),_p.py(),_p.pz(),_p.E()),
    mass(p.M()), pT(p.Pt()), y(p.Rapidity()), eta(p.Eta()), phi(p.Phi())
  { }
};
struct Jet: public Particle {
private:
  inline Double_t _tau(Double_t Y) noexcept {
    // need rapidity here
    return sqrt( pT*pT + mass*mass )/( 2.*cosh(y - Y) );
  }
public:
  Double_t tau;
  template<typename T>
  Jet(T&& p, Double_t Y) noexcept: Particle(p), tau(_tau(Y)) { }
};
// ******************************************************************

int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, wt_files, weights;
  string output_file, css_file, jet_alg, tree_name;
  size_t njets;
  double jet_pt_cut, jet_eta_cut;
  int_range<Long64_t> ents;
  bool counter_newline, quiet, strict;
  Long64_t cache_size;

  bool wt_given = false;

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
       "jet clustering algorithm: e.g. antikt4, kt6")
      ("jet-pt-cut", po::value<double>(&jet_pt_cut)->default_value(30.,"30"),
       "jet pT cut in GeV")
      ("jet-eta-cut", po::value<double>(&jet_eta_cut)->default_value(4.4,"4.4"),
       "jet eta cut")

      ("style,s", po::value<string>(&css_file)
       ->default_value(CONFDIR"/Ajets.css","Ajets.css"),
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
    if (vm.count("wt")) wt_given = true;
  }
  catch(exception& e) {
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
  unique_ptr<fastjet::JetDefinition> jet_def( fj_jetdef(jet_alg) );
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

  h_(A_pT) h_(A_y) h_(A_eta) h_(A_phi)

  h_jet_var(pT)

  h_(jets_HT)

  vector<hist_wt> h_Anj_pT; h_Anj_pT.reserve(njetsR);
  for (size_t j=1; j<=njetsR; ++j) {
    h_Anj_pT.emplace_back(cat('A',j,"j_pT"));
  }

  h_jet_var(y) h_jet_var(eta) h_jet_var(mass) h_jet_var(tau)

  h_(jets_tau_max) h_(jets_tau_sum)

  h_(jets_N_incl) h_(jets_N_excl)

  h_jj(jjpT_dy)   h_jj_excl(jjpT_dy)
  h_jj(jjfb_dy)   h_jj_excl(jjfb_dy)
  h_jj(jjpT_dphi) h_jj_excl(jjpT_dphi)
  h_jj(jjfb_dphi) h_jj_excl(jjfb_dphi)

  h_jj(jjpT_mass) h_jj(jjfb_mass);

  h_jj(A_jjpT_dy)   h_jj_excl(A_jjpT_dy)   h_jj(A_jjpT_dy_avgyjj)
  h_jj(A_jjfb_dy)   h_jj_excl(A_jjfb_dy)   h_jj(A_jjfb_dy_avgyjj)
  h_jj(A_jjpT_dphi) h_jj_excl(A_jjpT_dphi)
  h_jj(A_jjfb_dphi) h_jj_excl(A_jjfb_dphi)

  h_jj(A_jjfb_phi2)
  h_jj(AjjpT_mass) h_jj(Ajjfb_mass);

  h_jj(jjpT_N_jAj_incl) h_jj(jjfb_N_jAj_incl)
  h_jj(jjpT_N_jAj_excl) h_jj(jjfb_N_jAj_excl)

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

  TLorentzVector jjpT, jjfb;

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

    // Find photon (A)
    Int_t Ai = 0; // photon index
    while (Ai<event.nparticle) {
      if (event.kf[Ai]==22) break;
      else ++Ai;
    }
    if (Ai==event.nparticle) {
      if (!quiet)
        cerr << "\033[31mNo Higgs in event " << ent <<"\033[0m"<< endl;
      if (strict) return 1; else continue;
    }

    // Count number of events (not entries)
    if (prev_id!=event.eid) {
      h_N->Fill(0.5);
      prev_id = event.eid;
      ++num_events;

      for (auto h : hist_wt::all) h->FillSumw2();
    }

    const Particle A(
      TLorentzVector(event.px[Ai], event.py[Ai], event.pz[Ai], event.E[Ai]));

    // Jet clustering *************************************
    vector<Jet> jets;
    jets.reserve(njetsR);
    size_t nj;

    { // Clusted with FastJet on the fly
      vector<fastjet::PseudoJet> particles;
      particles.reserve(event.nparticle-1);

      for (Int_t i=0; i<event.nparticle; ++i) {
        if (i==Ai) continue;
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
        for (const auto& jet : fj_jets) jets.emplace_back(jet,A.y);

        // sort by pt
        sort(jets.begin(), jets.end(), [](const Jet& i, const Jet& j){
          return ( i.pT > j.pT );
        });

      }
    }
    // ****************************************************

    if (njets>1) {
      jjpT = jets[0].p + jets[1].p;
      jjpT_dy = abs(jets[0].y - jets[1].y);
      jjpT_mass = jjpT.M();
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

    h_A_pT .Fill(A.pT);
    h_A_y  .Fill(A.y);
    h_A_eta.Fill(A.eta);
    h_A_phi.Fill(A.phi);

    TLorentzVector Anj = A.p;
    for (size_t j=0, _nj=min(nj,njetsR); j<_nj; ++j) {
      h_jet_mass[j].Fill(jets[j].mass);
      h_jet_pT  [j].Fill(jets[j].pT  );
      h_jet_y   [j].Fill(jets[j].y   );
      h_jet_eta [j].Fill(jets[j].eta );
      h_jet_tau [j].Fill(jets[j].tau );

      h_Anj_pT  [j].Fill( (Anj += jets[j].p).Pt() );

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
      Double_t A_jj_dy   = abs(A.y - jjpT.Rapidity());
      Double_t A_jj_dphi = dphi(jjpT.Phi(),A.phi);

      h_jjpT_dy    .Fill(jjpT_dy);
      h_jjpT_dphi  .Fill(jjpT_dphi);
      h_A_jjpT_dy  .Fill(A_jj_dy);
      h_A_jjpT_dphi.Fill(A_jj_dphi);

      if (nj==njets) {
        h_jjpT_dy_nj_excl     .Fill(jjpT_dy);
        h_jjpT_dphi_nj_excl   .Fill(jjpT_dphi);

        h_A_jjpT_dy_nj_excl   .Fill(A_jj_dy);
        h_A_jjpT_dphi_nj_excl .Fill(A_jj_dphi);
      } else {
        h_jjpT_dy_nRj_excl    .Fill(jjpT_dy);
        h_jjpT_dphi_nRj_excl  .Fill(jjpT_dphi);

        h_A_jjpT_dy_nRj_excl  .Fill(A_jj_dy);
        h_A_jjpT_dphi_nRj_excl.Fill(A_jj_dphi);
      }

      h_jjpT_mass .Fill( jjpT_mass );
      h_AjjpT_mass.Fill( (A.p + jjpT).M() );

      h_A_jjpT_dy_avgyjj.Fill( abs(A.y - jjpT_ycenter) );

      h_jj_loose.Fill( 0.5 );
      if (A_jj_dphi > 2.6) h_jjpT_tight.Fill( 0.5 );

      size_t NjAj = ( between(jets[0].y,A.y,jets[1].y) ? nj : 0 );
      h_jjpT_N_jAj_excl.Fill(NjAj);
      h_jjpT_N_jAj_incl.FillIncl(NjAj);

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
      A_jj_dy   = abs(A.y - jjfb.Rapidity());
      A_jj_dphi = dphi(jjfb.Phi(),A.phi);

      h_jjfb_dy    .Fill(jjfb_dy);
      h_jjfb_dphi  .Fill(jjfb_dphi);
      h_A_jjfb_dy  .Fill(A_jj_dy);
      h_A_jjfb_dphi.Fill(A_jj_dphi);

      if (nj==njets) {
        h_jjfb_dy_nj_excl     .Fill(jjfb_dy);
        h_jjfb_dphi_nj_excl   .Fill(jjfb_dphi);

        h_A_jjfb_dy_nj_excl   .Fill(A_jj_dy);
        h_A_jjfb_dphi_nj_excl .Fill(A_jj_dphi);
      } else {
        h_jjfb_dy_nRj_excl    .Fill(jjfb_dy);
        h_jjfb_dphi_nRj_excl  .Fill(jjfb_dphi);

        h_A_jjfb_dy_nRj_excl  .Fill(A_jj_dy);
        h_A_jjfb_dphi_nRj_excl.Fill(A_jj_dphi);
      }

      h_jjfb_mass .Fill( jjfb_mass );
      h_Ajjfb_mass.Fill( (A.p + jjfb).M() );

      h_A_jjfb_dy_avgyjj.Fill( abs(A.y - jjfb_ycenter) );

      if (A_jj_dphi > 2.6) h_jjfb_tight.Fill( 0.5 );

      NjAj = ( between(jets[jb].y,A.y,jets[jb].y) ? nj : 0 );
      h_jjfb_N_jAj_excl.Fill(NjAj);
      h_jjfb_N_jAj_incl.FillIncl(NjAj);

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
        if (j.y > A.y) { vsumf += j.p; f_nonzero = true; }
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

      h_A_jjfb_phi2.Fill(phi2);

    } // END if (njets>1)

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
  delete wt_tree;

  return 0;
}
