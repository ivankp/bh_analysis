#include <iostream>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>

using namespace std;

#define test(var) \
cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

void scale(TDirectory* dir, double factor) noexcept {
  TIter nextkey(dir->GetListOfKeys());
  TKey *key;
  TObject *obj;
  while ((key = static_cast<TKey*>(nextkey()))) {
    obj = key->ReadObj();
    if (obj->InheritsFrom(TDirectory::Class())) {
      scale(static_cast<TDirectory*>(obj),factor);
    } else if (obj->InheritsFrom(TH1::Class())) {
      static_cast<TH1*>(obj)->Scale(factor);
    }
  }
}

int main(int argc, char** argv)
{
  if (argc!=3) {
    cout << "Usage: " << argv[0] << " factor file.root" << endl;
    return 1;
  }

  TFile *f = new TFile(argv[2],"update");
  if (f->IsZombie()) exit(1);
  cout << "Input file: " << f->GetName() << endl;

  scale(f,atof(argv[1]));

  f->Write(0,TObject::kOverwrite);

  delete f;
  return 0;
}
