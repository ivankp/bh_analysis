#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TChain.h>
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "BHEvent.hh"
#include "weight.hh"
#include "fj_jetdef.hh"
#include "int_range.hh"
#include "real_range.hh"
#include "timed_counter.hh"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
namespace po = boost::program_options;

int main(int argc, char** argv)
{
  // START OPTIONS **************************************************
  vector<string> bh_files, wt_files, weights;
  string output_file, jet_alg;
  size_t njets;
  double jet_pt_cut, jet_eta_cut;
  real_range<Double_t> AA_mass_cut;
  int_range<Long64_t> ents;
  bool AAntuple, counter_newline, quiet, strict;
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
      // ("output,o", po::value<string>(&output_file)->required(),
      //  "*output file with graph")

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

      ("AA", po::bool_switch(&AAntuple),
       "make Higgs from diphoton\nand produce AA histograms")
      ("AA-mass-cut", po::value<real_range<Double_t>>(&AA_mass_cut),
       "apply a mass cut to the diphoton,\ne.g. 115:135")

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
    if (AAntuple && vm.count("AA-mass-cut")) apply_AA_mass_cut = true;
  }
  catch(exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Setup input files **********************************************
  TChain* const bh_tree = new TChain("t3");
  TChain* const wt_tree = (wt_given ? new TChain("weights") : nullptr);

  // Add trees from all the files to the TChains
  cout << "BH files:" << endl;
  for (auto& f : bh_files) {
    cout << "  " << f << endl;
    if (!bh_tree->AddFile(f.c_str(),-1) ) exit(1);
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
    if (wt_given) if (need_ent>wt_tree->GetEntries()) {
      cerr << "Fewer entries in weights chain (" << wt_tree->GetEntries()
         << ") then requested (" << need_ent << ')' << endl;
      exit(1);
    }
  } else {
    ents.len = bh_tree->GetEntries();
    if (wt_given) if (ents.len!=wt_tree->GetEntries()) {
      cerr << ents.len << " entries in BH chain, but "
           << wt_tree->GetEntries() << " entries in weights chain" << endl;
      exit(1);
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
  }

  // Reading entries from the input TChain **************************
  Long64_t num_entries = 0, num_events = 0;
  Int_t prev_id = -1;
  cout << "Reading " << ents.len << " entries";
  if (ents.first>0) cout << " starting at " << ents.first;
  cout << endl;
  timed_counter counter(counter_newline);

  const size_t nw = weight::all.size();
  vector<pair<Double_t,Double_t>> cross_section(nw,{0.,0.});

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

    if (AAntuple) {
      TLorentzVector A1(event.px[Ai1], event.py[Ai1],
                        event.pz[Ai1], event.E [Ai1]);
      TLorentzVector A2(event.px[Ai2], event.py[Ai2],
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

      if (apply_AA_mass_cut) {
        Double_t AA_mass = (A1 + A2).M();
        if ( AA_mass < AA_mass_cut.min || AA_mass_cut.max < AA_mass ) continue;
      }
    }

    // Jet clustering *************************************
    { // Clusted with FastJet on the fly
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

      size_t nj = 0;

      // Apply eta cut
      for (auto& j : fj_jets)
        if (abs(j.eta()) <= jet_eta_cut) ++nj;

      if (nj < njets) continue; // cut on the number of jets
    }
    // ****************************************************

    if (wt_given) wt_tree->GetEntry(ent);

    // Count number of events (not entries)
    ++num_entries;
    if (prev_id!=event.eid) {
      prev_id = event.eid;
      ++num_events;

      if (prev_id>0) for (size_t i=0; i<nw; ++i) {
        cross_section[i].second += cross_section[i].first;
        cross_section[i].first = 0.;
      }
    }

    for (size_t i=0; i<nw; ++i) {
      cross_section[i].first += weight::all[i]->get();
    }

  } // END of event loop ********************************************
  ++num_events;

  for (size_t i=0; i<nw; ++i) {
    // add last event
    cross_section[i].second += cross_section[i].first;

    // convert to cross section
    cross_section[i].second /= num_events;

    cout << weight::all[i]->name << ": "
         << cross_section[i].second << " pb" << endl;
  }

  counter.prt(ents.end());
  cout << endl;
  cout << "Selected entries: " << num_entries << endl;
  cout << "Selected events : " << num_events  << endl;

  // Close files
  delete bh_tree;
  delete wt_tree;

  return 0;
}
