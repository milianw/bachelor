#include <iostream>

#include <Eigen/Dense>

#include "bench_shared.h"

using namespace std;
using namespace Eigen;

int main() {
  MatrixXcd m(DIMENSION, DIMENSION);
  for(int i = 0; i < DIMENSION; ++i) {
    m(i, i) = i;
    for(int j = i+1; j < DIMENSION; ++j) {
      m(i, j) = j;
      m(j, i) = j;
    }
  }
  DEBUG(cout << m << endl;)
  SelfAdjointEigenSolver<MatrixXcd> solver(m);
  DEBUG(cout << solver.eigenvectors() << endl;)
  cout << solver.eigenvalues().transpose() << endl;
  return 0;
}
