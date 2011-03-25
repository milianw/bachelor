#include "spinlib/spinhamiltonian.h"
#include "spinlib/experiment.h"
#include "spinlib/resonancefield.h"

#include <string>
#include <iostream>
#include <sstream>

using namespace std;

string formatSize(long long unsigned int size) {
  const char prefix[5] = {' ', 'K', 'M', 'G', 'T'};
  int i = 0;
  for(i; i < 4; ++i) {
    if (size > 1024) {
      size /= 1024;
    } else {
      break;
    }
  }
  stringstream stream;
  stream << size << prefix[i] << 'B';
  return stream.str();
}

int main(int argc, char** argv) {
  int spinHalf = 1;
  int spinOne = 0;
  if (argc > 1) {
    spinHalf = atoi(argv[1]);
  }
  if (argc > 2) {
    spinOne = atoi(argv[2]);
  }

  Experiment exp = Experiment::generateDummy(spinHalf, spinOne);
  cout << "peak mem consumption at least:" << formatSize(2 * (sizeof(complex<fp>) * exp.dimension * exp.dimension)) << endl;

  exp.mwFreqGHz = 9.5;
  {
  SpinHamiltonian H(0.3, exp);
  cout << H.calculateIntensity() << endl;
  }
  {
    const ResonanceField f(exp);
    const BisectNode a = f.diagonalizeNode(0.3);
    const BisectNode b = f.diagonalizeNode(0.6);
    f.checkSegment(a, b);
    f.findRootsInSegment(a, b);
  }
  return 0;
}
