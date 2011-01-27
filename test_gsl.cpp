#include <gsl/gsl_blas.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>

#include <Eigen/Dense>

#include <cmath>

void EIGEN_DONT_INLINE doStuff()
{
  const int n = 8;
  const int dim = pow(2, n);

  gsl_matrix_complex* m = gsl_matrix_complex_calloc(dim, dim);
  gsl_eigen_hermv_workspace* w = gsl_eigen_hermv_alloc(dim);
  gsl_vector* eigenvalues = gsl_vector_alloc(dim);
  gsl_matrix_complex* eigenvectors = gsl_matrix_complex_alloc(dim, dim);
  gsl_eigen_hermv(m, eigenvalues, eigenvectors,  w);
}

int main()
{
  doStuff();

  return 0;
}