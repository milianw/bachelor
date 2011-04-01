#include "spinlib/spinhamiltonian.h"
#include "spinlib/experiment.h"
#include "spinlib/resonancefield.h"
#include "spinlib/spins.h"
#include "spinlib/constants.h"
#include "spinlib/helpers.h"

#include <string>
#include <iostream>
#include <sstream>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "spinlib/operatorsum_p.h"

void dbg(const MatrixXc &m)
{
  SelfAdjointEigenSolver<MatrixXc> eigenSolver(m, EigenvaluesOnly);
  cout << m << endl;
  cout << "determinant:" << m.determinant() << endl;
  cout << "Eigen values:" << eigenSolver.eigenvalues().transpose() << endl;
  cout << endl;
}

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
  if (argc > 3) {
    B = atof(argv[3]);
  }

  cout << "number of J = 0.5: " << spinHalf << " + 1 electron" << endl;
  cout << "number of J = 1: " << spinOne << endl;
  cout << "B = " << B << "T" << endl;

  Experiment exp = getExperiment(getenv("ORCA_FILE") ? getenv("ORCA_FILE") : "", spinHalf, spinOne);
  printExperiment(cout, exp);
  SpinHamiltonian H(B, exp);
  {
    MatrixXc m = H.hamiltonian();
    m /= (Constants::h * 1E09);
    cout << "total hamiltonian:" << endl;
    dbg(m);
  }
  {
    MatrixXc m(exp.dimension, exp.dimension);
    m.setZero();
    H.addNuclearZeeman(m);
    m /= (Constants::h * 1E09);
    cout << "nuclear zeeman:" << endl;
    dbg(m);
  }
  {
    MatrixXc m(exp.dimension, exp.dimension);
    m.setZero();
    H.addElectronZeeman(m);
    m /= (Constants::h * 1E09);
    cout << "electron zeeman:" << endl;
    dbg(m);
  }
  {
    MatrixXc m(exp.dimension, exp.dimension);
    m.setZero();
    H.addNuclearZeeman(m);
    H.addElectronZeeman(m);
    m /= (Constants::h * 1E09);
    cout << "zeeman total:" << endl;
    dbg(m);
  }
  {
    MatrixXc m(exp.dimension, exp.dimension);
    m.setZero();
    H.addHyperFine(m);
    m /= (Constants::h * 1E09);
    cout << "hyperfine:" << endl;
    dbg(m);
  }
  /*
  {
    Spins s(exp.nProtons + 1, exp.nNitrogens);
    cout << "spin halfs:" << s.spinHalfs << endl;
    cout << "pauli matrices:" << endl << PauliMatrix_J_one::X << endl << PauliMatrix_J_one::Y << endl  << PauliMatrix_J_one::Z << endl;
    for(int bra = 0; bra < s.states; ++bra) {
      for(int ket = 0; ket < s.states; ++ket) {
        cout << bra << '\t' << ket << "\t:\t" << s.spinVector(bra, ket, 1).transpose() << endl;
      }
    }
  }
  */
  /*
  {
    SpinOperator<Matrix2c> S(PauliMatrix_J_half::X, PauliMatrix_J_half::Y, PauliMatrix_J_half::Z);
    SpinOperator<Matrix2c> I(PauliMatrix_J_half::X, PauliMatrix_J_half::Y, PauliMatrix_J_half::Z);
    SpinOperator<OperatorSum> I_total;

    I_total += S;
    I_total += S;

    cout << (I_total.X.toMatrix() + I_total.Y.toMatrix() + I_total.Z.toMatrix()) << endl;
  }
  {
    Spins spins(2, 0);
    for(int bra = 0; bra < spins.states; ++bra) {
      for(int ket = 0; ket < spins.states; ++ket) {
        cout << bra << '\t' << ket << ":" << endl;
        const Vector3c S = spins.spinVector(bra, ket, 0);
        const Vector3c I = spins.spinVector(bra, ket, 1);
        cout << '\t'  << spins.spinInState(0, bra) << '\t' << spins.spinInState(0, ket) << '\t' << S.transpose() << " (Electron)" << endl;
        cout << '\t'  << spins.spinInState(1, bra) << '\t' << spins.spinInState(0, ket) << '\t'  << I.transpose() << " (Nucleus)" << endl;
        cout << "\t \t \t" << (S(0) * I(0)) << ' ' << (S(1) * I(1)) << (S(2) * I(2)) <<  endl;
        cout << "\t \t \t" << (conj(S(0)) * I(0)) << ' ' << (conj(S(1)) * I(1)) << (conj(S(2)) * I(2)) <<  endl;
        cout << '\t'  << S.dot(I) << '\t' << I.dot(S) << (I.transpose() * S) << endl;
        cout << '\t'  << (S(0) * I(0) + S(1) * I(1) + S(2) * I(2)) << endl;
        cout << '\t'  << (conj(S(0)) * I(0) + conj(S(1)) * I(1) + conj(S(2)) * I(2)) << endl;
        cout << endl;
      }
    }
  }
  */
  return 0;
}
