#include <iostream>

#include <boost/foreach.hpp>
#include <boost/regex.hpp>

#include <spinlib/orcaparser.h>
#include <spinlib/nucleus.h>

#include <vector>
#include <algorithm>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv) {
  if (argc < 2) {
    cerr << "orca FILE missing" << endl;
    return 1;
  }

  vector<boost::regex> allowedNuclei;
  for(int i = 2; i < argc; ++i) {
    allowedNuclei.push_back(boost::regex(argv[i]));
  }

  OrcaParser parser(argv[1]);

  cout << "-------------------" << endl;
  cout << "ELECTRONIC G-MATRIX" << endl;
  cout << "-------------------" << endl;
  cout << endl;
  cout << "The g-matrix: " << endl;
  cout << parser.electronGMatrix() << endl;
  cout << endl;
  cout << "-----------------------------------------" << endl;
  cout << "ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE" << endl;
  cout << "-----------------------------------------" << endl;
  cout << endl;

  BOOST_FOREACH(const Nucleus& nuc, parser.nuclei()) {
    if (!allowedNuclei.empty()) {
      bool found = false;
      BOOST_FOREACH(const boost::regex& pattern, allowedNuclei) {
        if (boost::regex_match(nuc.name, pattern)) {
          found = true;
          break;
        }
      }
      if (!found) {
        continue;
      }
    }
    cout << "-----------------------------------------------------------" << endl;
    cout << " Nucleus  " << nuc.name << "  : A:ISTP=    " << nuc.isotope << " I=  " << float(nuc.twoJ)/2.0 << " P=0.0 MHz/au**3" << endl;
    cout << "                    Q:ISTP=    0 I=  0.0 Q=  0.0 barn" << endl;
    cout << "-----------------------------------------------------------" << endl;
    cout << endl;
    cout << " Raw HFC matrix (all values in MHz): " << endl;
    cout << nuc.A << endl;
    cout << endl;
  }

}