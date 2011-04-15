#include <iostream>

#include <boost/foreach.hpp>

#include "spinlib/orcaparser.h"
#include "spinlib/nucleus.h"
#include "spinlib/experiment.h"
#include "spinlib/helpers.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv) {
  if (argc < 2) {
    cerr << "orca FILE missing" << endl;
    return 1;
  }

  OrcaParser parser(argv[1]);

  cout << "electron G Matrix:" << endl << parser.electronGMatrix() << endl << endl;

  BOOST_FOREACH(const Nucleus& nuc, parser.nuclei()) {
    cout << "Nucleus: " << nuc.name << ", J = " << float(nuc.twoJ)/2 << ", g = " << nuc.g << endl;
    cout << "A = " << endl << nuc.A << endl;
    cout << "euler angles: " << parser.eulerRotation()[nuc.name].first.transpose() << endl;
    cout << "[Ax, Ay, Az]: " << parser.eulerRotation()[nuc.name].second.transpose() << endl;
    if (!nuc.Q.isZero()) {
      cout << "Q = " << nuc.Q.transpose() << endl;
    }
    if (!nuc.EFG.isZero()) {
      cout << "EFG = " << endl << nuc.EFG << endl;
      cout << "Q_full = " << endl << nuc.quadrupoleCouplingMatrix() << endl;
    }
    cout << endl;
  }

  Experiment exp(parser.nuclei());
  cout << "peak mem consumption per node for this orca file approx " << guessPeakMemConsumption(exp) << endl;

  return 0;
}
