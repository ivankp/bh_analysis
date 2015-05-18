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

inline Double_t dphi(const TLorentzVector& v1, Double_t phi2) noexcept {
  return dphi(v1.Phi(),phi2);
}

inline Double_t fphi2(const TLorentzVector& b, const TLorentzVector& f) noexcept {
  const Double_t bx=b.Px(), by=b.Py(), fx=f.Px(), fy=f.Py();

  const Double_t phi2 = acos( ( bx*fx+by*fy ) /
    sqrt(bx*bx+by*by) * sqrt(fx*fx+fy*fy)
  );
  if ( ( bx*fy - by*fx ) < 0.) return -phi2;
  else return phi2;
}

senum(jjtype,(pT)(fb))

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

int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, sj_files, wt_files, weights;
  string output_file, css_file, jet_alg;
  size_t njets;
  double jet_pt_cut, jet_eta_cut;
  real_range<Double_t> AA_mass_cut;
  int_range<Long64_t> ents;
  bool AAntuple, VBF, counter_newline, quiet;
  jjtype::type jj_type = jjtype::pT;

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

      ("jj", po::value<jjtype::type>(&jj_type),
       "tagging jets type: pT or fb")

      ("VBF", po::bool_switch(&VBF),
       "apply Vector Bososn Fusion cuts\njj_dy>2.8 && jj_mass>400")

      ("AA", po::bool_switch(&AAntuple),
       "make Higgs from diphoton\nand produce AA histograms")
      ("AA-mass-cut", po::value<real_range<Double_t>>(&AA_mass_cut),
       "apply a mass cut to the diphoton,\ne.g. 115:135")

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
    if (vm.count("AA_mass_cut")) apply_AA_mass_cut = true;
    if (njets>2 && !vm.count("jj")) {
      throw runtime_error(
        "Type of two tagging jets has to be specified for njets > 2\n"
        "Pass either --jj=pT or --jj=fb flag");
    }
    if (njets<2) VBF = false;
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

  h_(H_pT); h_(H_y); h_(H_mass);

  h_opt(AA_cos_theta_star, AAntuple);
  h_opt(AA_pTt, AAntuple);
  h_opt(AA_dy,  AAntuple);

  h_jet_var(pT);

  h_(jets_HT);

  vector<hist_wt> h_Hnj_pT; h_Hnj_pT.reserve(njetsR);
  for (size_t j=1; j<=njetsR; ++j) {
    h_Hnj_pT.emplace_back(cat('H',j,"j_pT"));
  }

  h_jet_var(y); h_jet_var(mass); h_jet_var(tau);

  h_(jets_tau_max); h_(jets_tau_sum);

  h_(jets_N_incl_nojetcuts); h_(jets_N_excl_nojetcuts);

  h_jj(jj_dy);
  hist_wt h_jj_dy_nj_excl ( njets>1 ? cat("jj_dy_",njets, "j_excl") : string() ),
          h_jj_dy_nRj_excl( njets>1 ? cat("jj_dy_",njetsR,"j_excl") : string() );

  h_jj(jj_dphi);
  h_jj(jj_mass);

  h_jj(Hjj_mass);
  h_jj(H_jj_dy); h_jj(H_jj_dy_avgyjj); h_jj(H_jj_dphi);
  h_opt(H_jj_phi2, njets>1 && jj_type==jjtype::fb);

  h_jj(jj_N_jHj_excl);

  const size_t ndy = 6;

  vector<vector<hist_wt>> h_jet_pT_jj(njetsR);
  if (njets>1) for (size_t j=0; j<njetsR; ++j) {
    h_jet_pT_jj[j].reserve(ndy);
    for (size_t i=0; i<ndy; ++i)
      h_jet_pT_jj[j].emplace_back(cat("jet",j+1,"_pT_minjjdy",i+1));
  }

  h_jj(jj_loose); h_jj(jj_tight);

  h_opt(jj_j_dy_veto, njets>2);

  // Reading entries from the input TChain **************************
  Long64_t num_selected = 0, num_events = 0;
  Int_t prev_id = -1;
  cout << "Reading " << ents.len << " entries";
  if (ents.first>0) cout << " starting at " << ents.first;
  cout << endl;
  timed_counter counter(counter_newline);

  // Variables ******************************************************
  size_t j1 = 0, j2 = 0;

  TLorentzVector A1, A2, higgs, jj;

  Double_t H_mass, H_pT, H_y, H_phi;

  Double_t jj_dy=0, jj_dphi=0, jj_ycenter=0, jj_mass=0;

  // LOOP ***********************************************************
  for (Long64_t ent = ents.first, ent_end = ents.end(); ent < ent_end; ++ent) {
    counter(ent);
    tree->GetEntry(ent);

    if (event.nparticle>BHMAXNP) {
      cerr << "\033[31mMore particles in the entry then BHMAXNP\033[0m" << endl
           << "Increase array length to " << event.nparticle << endl;
      exit(1);
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
        cerr << "\033[31mEvent " << ent
             << " doesn't have 2 photons\033[0m" << endl;
        continue;
      }

    } else { // Higgs

      while (Hi<event.nparticle) {
        if (event.kf[Hi]==25) break;
        else ++Hi;
      }
      if (Hi==event.nparticle) {
        cerr << "\033[31mNo Higgs in event " << ent <<"\033[0m"<< endl;
        continue;
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

    // Number of jets hists
    h_jets_N_excl_nojetcuts.Fill(nj);
    for (size_t i=0; i<=nj; ++i)
      h_jets_N_incl_nojetcuts.Fill(i);

    // Do no process entries with fewer then njets
    if (nj < njets) continue;

    if (njets==2) {
      j1 = 0; j2 = 1;
    } else if (njets>2) {
      switch (jj_type) {
        case jjtype::pT:
          j1 = 0; j2 = 1;
          break;
        case jjtype::fb:
          for (size_t j=1; j<nj; ++j) {
            if (jets[j].y < jets[j1].y) j1 = j;
            if (jets[j].y > jets[j2].y) j2 = j;
          }
          break;
      }
    }

    if (njets > 1) {
      jj = jets[j1].p + jets[j2].p;
      jj_dy = abs(jets[j1].y - jets[j2].y);
      jj_mass = jj.M();

      if (VBF) // Apply VBF cuts
        if (jj_dy>2.8 && jj_mass>400) continue;

      jj_dphi = dphi(jets[j1].phi, jets[j2].phi);
      jj_ycenter = (jets[j1].y + jets[j2].y)/2;
    }

    // Increment selected entries
    ++num_selected;

    // ****************************************************

    H_pT  = higgs.Pt();  // Higgs Pt
    H_phi = higgs.Phi(); // Higgs Phi

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
        for (size_t i=0; i<ndy; ++i)
          if ( jj_dy > (i+1) )
            h_jet_pT_jj[j][i].Fill(jets[j].pT);
      }
    }

    if (njets>1) {

      h_jj_dy  .Fill(jj_dy);
      h_jj_dphi.Fill(jj_dphi);

      if (nj==njets) h_jj_dy_nj_excl .Fill(jj_dy);
      else           h_jj_dy_nRj_excl.Fill(jj_dy);

      const Double_t H_jj_dphi = dphi(jj,H_phi);

      h_jj_mass  .Fill( jj_mass );
      h_Hjj_mass .Fill( (higgs + jj).M() );
      h_H_jj_dy  .Fill( abs(H_y - jj.Rapidity()) );
      h_H_jj_dphi.Fill( H_jj_dphi );
      h_H_jj_dy_avgyjj.Fill( abs(H_y - jj_ycenter) );

      h_jj_loose .Fill( 0.5 );
      if (H_jj_dphi > 2.6) h_jj_tight.Fill( 0.5 );

      h_jj_N_jHj_excl.Fill( between(jets[j1].y,H_y,jets[j2].y) ? nj : 0 );

	    if (njets>2) {
        // Third photon veto:
        // ydists is the smallest distance between the centre of the
        // tagging jets and any possible further jet
        Double_t y_dists = 100.;
        for (size_t j=0; j<nj; ++j) {
          if (j==j1 || j==j2) continue;
          Double_t y_distt = fabs(jets[j].y - jj_ycenter);
          if (y_distt < y_dists) y_dists = y_distt;
        }
        h_jj_j_dy_veto.FillIncl(y_dists);
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
      if (jj_type == jjtype::fb) {
        // j1 is backward, j2 is forward
  	    Double_t phi2;
  	    if (f_nonzero && b_nonzero) {
  	      phi2 = fphi2(vsumb,vsumf);
  	    } else if (!f_nonzero) {
  	      vsumb -= jets[j2].p;
  	      phi2 = fphi2(vsumb,jets[j2].p);
  	    } else {
  	      vsumf -= jets[j1].p;
  	      phi2 = fphi2(jets[j1].p,vsumf);
  	    }

  	    h_H_jj_phi2.Fill(phi2);
      }

		  // Diphoton histograms
		  if (AAntuple) {

	      // |costheta*| from 1307.1432
	      h_AA_cos_theta_star.Fill(
          abs(sinh(A1.Eta()-A2.Eta())) / sqrt(1.+sq(H_pT/H_mass)) *
          A1.Pt()*A2.Pt()*2 / sq(H_mass)
        );

        h_AA_pTt.Fill(
          fabs(A1.Px()*A2.Py()-A2.Px()*A1.Py()) / ((A1-A2).Pt()*2)
        );

        h_AA_dy.Fill( fabs(A1.Rapidity() - A2.Rapidity()) );

		  } // END AA

    } // END if (njets>1)

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
