#ifndef __fj_jetdef_hh__
#define __fj_jetdef_hh__

#include <string>
#include <algorithm>
#include <stdexcept>

#include <fastjet/ClusterSequence.hh>

fastjet::JetDefinition* fj_jetdef(const std::string& str) {
  std::string::iterator it = --str.end();
  while (isdigit(*it)) --it;
  ++it;

  std::string name;
  std::transform(str.begin(), it, std::back_inserter(name), ::tolower);
  fastjet::JetAlgorithm alg;
  if (!name.compare("antikt")) alg = fastjet::antikt_algorithm;
  else if (!name.compare("kt")) alg = fastjet::kt_algorithm;
  else if (!name.compare("cambridge")) alg = fastjet::cambridge_algorithm;
  else throw runtime_error("Undefined jet clustering algorithm: "+name);

  return new fastjet::JetDefinition(
    alg,
    std::stod( std::string(it,str.end()) )/10.
  );
}

#endif
