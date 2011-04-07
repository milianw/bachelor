#include "spinlib/spinhamiltonian.h"
#include "spinlib/experiment.h"
#include "spinlib/resonancefield.h"
#include "spinlib/spins.h"
#include "spinlib/helpers.h"

#include <string>
#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace boost::filesystem;

int main(int argc, char** argv) {
  int spinHalf = 1;
  int spinOne = 0;
  if (argc > 1) {
    spinHalf = atoi(argv[1]);
  }
  if (argc > 2) {
    spinOne = atoi(argv[2]);
  }

  path orcaFilePath(argv[1]);
  string orcaFile;
  if (exists(orcaFilePath)) {
    orcaFile = argv[1];
  }

  Experiment exp = getExperiment(orcaFile, spinHalf, spinOne);
  exp.mwFreqGHz = 9.5;

  ResonanceField f(exp);
  fp from;
  fp to;
  guessBRange(exp, from, to);
  vector< fp > field = f.calculate(from, to);
  for(unsigned int i = 0; i < field.size(); ++i) {
    cout << field.at(i);
    if (i + 1 < field.size())
      cout << ", ";
  }
  cout << endl;
  return 0;
}
