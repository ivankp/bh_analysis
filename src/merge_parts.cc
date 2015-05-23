#include <iostream>
#include <string>
#include <unordered_map>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>

using namespace std;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

// Safely get an object from TDirectory
template<class T>
inline T* get(TDirectory* dir, const char* name) {
  T *obj = nullptr;
  dir->GetObject(name,obj);
  if (obj) return obj;
  else {
    cerr << "No object " << name << " in " << dir->GetName() << endl;
    exit(1);
  }
}

int main(int argc, char** argv)
{
  if (argc<3) {
    cout << "Usage: " << argv[0]
         << " output_fout.root input_born.root input_real.root ..." << endl;
    exit(0);
  }

  // Track histograms
  unordered_map<string,int> hist_names;

  // Output file
  TFile *fout = new TFile(argv[1],"recreate");
  if (fout->IsZombie()) exit(1);
  cout << "Output file: " << fout->GetName() << endl;

  for (int fi=2; fi<argc; ++fi) {
    static bool first_file = true;

    // Current input file
    TFile *f = new TFile(argv[fi],"read");
    if (f->IsZombie()) exit(1);
    cout << "Input file: " << f->GetName() << endl;

    // Number of events
    const Double_t N = get<TH1D>(f,"N")->GetAt(1);
    const Double_t N_scale = 1./N;
    cout << "Events: " << N << endl;

    static TKey *key1;
    TIter nextkey1(f->GetListOfKeys());
    while ((key1 = (TKey*)nextkey1())) { // loop over dirs
      static TObject *obj1;
      obj1 = key1->ReadObj();
      if (obj1->InheritsFrom(TDirectory::Class())) {

        static bool first_dir = true;

        TDirectory *dir = static_cast<TDirectory*>(obj1);
        TDirectory *out_dir =
          ( first_file ? fout->mkdir(dir->GetName())
                       : get<TDirectory>(fout,dir->GetName()) );
        if (first_file) out_dir->cd();

        static TKey *key2;
        TIter nextkey2(dir->GetListOfKeys());
        while ((key2 = (TKey*)nextkey2())) { // loop over hists
          static TObject *obj2;
          obj2 = key2->ReadObj();
          if (obj2->InheritsFrom(TH1::Class())) {

            TH1 *h = static_cast<TH1*>(obj2);

            if (first_dir) {
              hist_names[h->GetName()] = 1;
            } else {
              const auto it = hist_names.find(h->GetName());
              if (it==hist_names.end()) {
                cerr << "\033[33mWarning\033[0m:"
                        " Unexpected histogram " << h->GetName()
                     << " in file " << f->GetName()
                     << " in dir " << dir->GetName() << endl;
                // exit(1);
                continue;
              } else {
                ++it->second;
              }
            }

            // Scale and add histograms
            if (first_file) {
              h->Scale(N_scale);
              h->Clone();
            } else {
              TH1 *new_h = get<TH1>(out_dir,h->GetName());
              Double_t ent = new_h->GetEntries();
              new_h->Add(h,N_scale);
              new_h->SetEntries(ent+h->GetEntries());
            }

          }
        }

        // Check if all histograms from the first directory were present
        auto it = hist_names.begin();
        const int first_counts = it->second;
        for (auto end=hist_names.end(); it!=end; ++it) {
          if (it->second != first_counts) {
            cerr << "\n\033[33mWarning\033[0m: Expected hist "
                 << (first_counts > it->second ? it->first
                                               : hist_names.begin()->first)
                 << " is not in file " << f->GetName()
                 << " in dir " << dir->GetName();
            // exit(1);
            continue;
          }
        }

        // no longer first directory
        if (first_dir) first_dir = false;

      } else { // not a TDirectory
        if (obj1->InheritsFrom(TH1::Class())) {

            TH1 *h = static_cast<TH1*>(obj1);

            // Scale and add histograms
            if (first_file) {
              fout->cd();
              h->Scale(N_scale);
              h->Clone();
            } else {
              TH1 *new_h = get<TH1>(fout,h->GetName());
              Double_t ent = new_h->GetEntries();
              new_h->Add(h,N_scale);
              new_h->SetEntries(ent+h->GetEntries());
            }

        }
      }
    }

    delete f;
    // no longer first file
    if (first_file) first_file = false;
  }

  // Set N to 1
  get<TH1D>(fout,"N")->SetAt(1,1);

  // Write and close output root file
  fout->Write();
  cout << "\n\033[32mWrote file\033[0m: " << fout->GetName() << endl;
  delete fout;

  return 0;
}
