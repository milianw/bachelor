#include <iostream>

#include <cstdio>

#include <mkl.h>
#include <mkl_lapack.h>

#include "bench_shared.h"

using namespace std;

void PrintArrayZ(MKL_INT m, MKL_INT n, MKL_Complex16 *a, MKL_INT *lda, char *name)
{
  MKL_INT   i, j;
  MKL_Complex16 *addr;

  printf("\n       ARRAY %s", name);

  for (i = 0; i < m; i++) {
    printf("\n         ");
    addr = a + i*(*lda);
    for (j = 0; j < n; j++)
      printf("(%6.2f,%6.2f)   ", (addr+j)->real,
            (addr+j)->imag);
  } /* for */
  printf("\n");
  return;
} /* PrintArrayZ */

int main() {
  int rmaxa = DIMENSION + 1;
  int cmaxa = DIMENSION;
  // row major:
  int lda = cmaxa;

  MKL_Complex16* m = (MKL_Complex16 *)calloc( rmaxa * cmaxa, sizeof(MKL_Complex16) );
  MKL_Complex16 c;
  for(int i = 0; i < DIMENSION; ++i) {
    c.imag = 0;
    c.real = i;
    // (i, i):
    *(m + i*lda + i) = c;
    for(int j = i+1; j < DIMENSION; ++j) {
      c.imag = 0;
      c.real = j;
      // (i, j):
      *(m + i*lda + j) = c;
      // (j, i):
      *(m + j*lda + i) = c;
    }
  }

  DEBUG(PrintArrayZ(DIMENSION, DIMENSION, m, &lda, "m");)

  int n = DIMENSION;
  int lwork = n;
  double *d = new double[n];
  double *e = new double[n];
  int info = 0;
  {
  MKL_Complex16* tau = (MKL_Complex16 *)calloc( (rmaxa-1) * (cmaxa-1), sizeof(MKL_Complex16) );
  MKL_Complex16* work = (MKL_Complex16 *)calloc( lwork * (lwork+1), sizeof(MKL_Complex16) );
  char uplo = 'L';
  zhetrd(&uplo, &n, m, &lda, d, e, tau, work, &lwork, &info);
  DEBUG(cout << "info after zhetrd: " << info << endl;)
  zungtr(&uplo, &n, m, &lda, tau, work, &lwork, &info);
  DEBUG(cout << "info after zungtr: " << info << endl;)
  DEBUG(PrintArrayZ(DIMENSION, DIMENSION, m, &lda, "m_zhetrd");)
  cfree(tau);
  cfree(work);
  }

  {
  char compz = 'V';
  double* work = new double[2*n-2];
  zsteqr(&compz, &n, d, e, m, &lda, work, &info);
  DEBUG(PrintArrayZ(DIMENSION, DIMENSION, m, &lda, "m_vec");)
  DEBUG(cout << "info after zsteqr: " << info << endl;)
  cout << "eigen values: ";
  for(int i = 0; i < DIMENSION; ++i) {
    cout << d[i] << "  ";
  }
  cout << endl;
  delete[] work;
  }

  cfree(m);
  delete[] d;
  delete[] e;
  return 0;
}
