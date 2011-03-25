#include "spinlib/spinhamiltonian.h"
#include "spinlib/experiment.h"
#include "spinlib/resonancefield.h"
#include "spinlib/spins.h"

#include <string>
#include <iostream>
#include <sstream>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "spinlib/operatorsum_p.h"

int main(int argc, char** argv) {
  int spinHalf = 1;
  int spinOne = 0;
  float B = 0.3;
  if (argc > 1) {
    spinHalf = atoi(argv[1]);
  }
  if (argc > 2) {
    spinOne = atoi(argv[2]);
  }

  cout << "number of J = 0.5: " << spinHalf << " + 1 electron" << endl;
  cout << "number of J = 1: " << spinOne << endl;

  Experiment exp = Experiment::generateDummy(spinHalf, spinOne);
  exp.mwFreqGHz = 9.5;

  ResonanceField f(exp);
  vector< fp > field = f.calculate(0, 1);
  for(int i = 0; i < field.size(); ++i) {
    cout << field.at(i);
    if (i + 1 < field.size())
      cout << ", ";
  }
  cout << endl;
  return 0;
}
