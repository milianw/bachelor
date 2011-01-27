#include <Eigen/Dense>

#include <cmath>

using namespace Eigen;

void EIGEN_DONT_INLINE doStuff()
{
  const int n = 8;
  const int dim = pow(2, n);

  MatrixXcd m(dim, dim);
  SelfAdjointEigenSolver<MatrixXcd> eigenSolver(m);
  const VectorXd eigenValues = eigenSolver.eigenvalues();
  const MatrixXcd eigenVectors = eigenSolver.eigenvectors();
}

int main()
{
  doStuff();

  return 0;
}