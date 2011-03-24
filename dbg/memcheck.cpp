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

int main() {
  Experiment exp(4, 1);
  exp.mwFreqGHz = 9.5;
  {
  SpinHamiltonian H(0.3, exp);
  cout << "peak mem consumption at least:" << formatSize(2 * (sizeof(complex<fp>) * exp.dimension * exp.dimension)) << endl;
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
