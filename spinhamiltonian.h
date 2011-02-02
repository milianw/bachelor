/*
 * This file is part of my bachelor thesis.
 *
 * Copyright 2011 Milian Wolff <mail@milianw.de>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef MW_BACHELOR_SPINHAMILTONIAN_H
#define MW_BACHELOR_SPINHAMILTONIAN_H

#include <iostream>
#include <QtCore/QTextStream>
#include <string>
#include <cmath>
#include <bitset>

#include <complex>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "constants.h"
#include "experiment.h"

using namespace Constants;

typedef complex<double> c_double;

// Pauli Matrices
namespace PauliMatrix {
const Matrix2cd X = (Matrix2cd() << 0, 0.5, 0.5, 0).finished();
const Matrix2cd Y = (Matrix2cd() << 0, c_double(0, -0.5), c_double(0, 0.5), 0).finished();
const Matrix2cd Z = (Matrix2cd() << 0.5, 0, 0, -0.5).finished();
}

/**
 * TODO: check whether inlining some parts is noticeable
 *
 * notes on porting:
 *
 * ~~~~~~~~~~~
 * gsl_blas_zdotu (x, y, dotu)
 * => eigen: dotu = x.conjugate().dot(y)
 *
 * from eigen docs about .dot():
 * Note: If the scalar type is complex numbers, then this function returns the hermitian (sesquilinear) dot product,
 *       conjugate-linear in the first variable and linear in the second variable.
 *       but gsl_blas_zdotu seems to differ from this, hence use a.conjugate().dot(b) instead of a.dot(b)
 * NOTE: <orzel> in this very specific case, it might be that (a.transpose()*b)(0,0) is faster (taking the only element of the 1x1 matrix x^T.y
 * ~~~~~~~~~~~
 * gsl_blas_zgemv(CblasTrans, gsl_complex_rect(1,0), atensor, I, gsl_complex_rect(0,0), atensor_I);
 * => eigen: atensor_I = atensor.transpose() * I
 * ~~~~~~~~~~~
 * gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1,0),
 *                eigenvectors, moments, gsl_complex_rect(1,0), intermediate);
 * docs: C = \alpha op(A) op(B) + \beta C
 * => eigen: intermediate = eigenvectors.adjoint() * moments + intermediate
 */
class SpinHamiltonian {
  public:
    /// @p B magnetic field in Tesla
    SpinHamiltonian(const double B, const Experiment& experiment);
    ~SpinHamiltonian();

    /// list all possible transitions with their frequency for the given B field
    void calculateTransitions() const;
    /// calculate total intensity of all transitions that are valid for the
    /// incoming microwave frequency
    void calculateIntensity(const double mwFreq, QTextStream* out) const;

  private:
    /// the complete hamiltonian
    MatrixXcd hamiltonian() const;
    /// nuclear Zeeman Hamiltonian component
    MatrixXcd nuclearZeeman() const;
    /// hyper fine Hamiltonian component
    MatrixXcd hyperFine() const;
    /// electron Zeeman Hamiltonian component
    MatrixXcd electronZeeman() const;

    /// interprets @p i as binary number and returns the k-th bit of it
    inline bool spinState(int i, int k) const;

    /// return spin vector from pauli matrices
    inline Vector3cd spinVector(int i, int j, int k) const;

    /// all bits for states @p i, @p j, must match for nucleus @p k
    /// otherwise the integral will be zero anyways and we can skip
    /// <  i  |H|  j  >
    /// <10101|H|10001> = <1|H|0><1|1><0|0><0|0><1|1> => must be calculated
    ///    ^=k     ^=k
    /// <10100|H|10001> = <1|H|0><1|1><0|0><0|0><0|1>
    ///    ^=k     ^=k                            ^=0 => can be skipped
    enum {
      IncludeElectron = 0,
      IgnoreElectron = 1
    };
    inline bool stateContributes(int i, int j, int k, bool ignoreElectron = IgnoreElectron) const;

    /// moments
    /// TODO: better document
    MatrixXcd magneticMoments() const;

    inline c_double magneticMoment(const int i, const int j) const;

    /// intensitiy matrix with coefficients (i, j) = |< psi_j | M | psi_i>|^2
    /// psi_i being the i-th eigen vector
    /// M being the magnetic moment matrix
    MatrixXd intensityMatrix(const MatrixXcd& eigenVectors) const;

    const double m_B;
    const Experiment m_exp;
    Vector3cd m_staticBField;
};

SpinHamiltonian::SpinHamiltonian(const double B, const Experiment& experiment)
: m_B(B)
, m_exp(experiment)
, m_staticBField(m_exp.staticBField(B))
{
}

SpinHamiltonian::~SpinHamiltonian()
{
}

inline Vector3cd SpinHamiltonian::spinVector(int i, int j, int k) const
{
  const int a = spinState(i, k); //spin state of state k in row i
  const int b = spinState(j, k); //spin state of state k in column j
  return (Vector3cd() << PauliMatrix::X(a, b), PauliMatrix::Y(a, b), PauliMatrix::Z(a, b)).finished();
}

bool SpinHamiltonian::spinState(int i, int k) const
{
  // k-bit == 2^k = 0001000
  //                   ^k = 4
  const int kPow = (1 << k);
  return i & kPow;
}

inline bool SpinHamiltonian::stateContributes(int i, int j, int k, bool ignoreElectron) const
{
  // states are equal if: all bits except for k-bit are equal
  // k-bit == 2^k = 0001000
  //                   ^k = 4
  int kBit = (1 << k);
  if (ignoreElectron) {
    // ignore electron bit
    // electronBit == 2^nProtons == 100000...
    int electronBit = (1 << m_exp.nProtons);
    // essentially a fast variant of:
    // i % electronBit
    i &= electronBit - 1;
    j &= electronBit - 1;
  }
  return (i | kBit) == (j | kBit);
}

MatrixXcd SpinHamiltonian::hamiltonian() const
{
  return nuclearZeeman() + hyperFine() + electronZeeman();
}

MatrixXcd SpinHamiltonian::nuclearZeeman() const
{
  //Compute nZeeman============================================================  
  MatrixXcd nZeeman(m_exp.dimension, m_exp.dimension);
  nZeeman.setZero();
  //to turn off: return nZeeman;

  for (int i = 0; i < m_exp.dimension; ++i) {
    for (int j = 0; j < m_exp.dimension; ++j) {
      //m_exp.nProtons is always the index of the electronic spin state

      if (spinState(i, m_exp.nProtons) != spinState(i, m_exp.nProtons)) {
        continue;  //matrix elements between different e states are zero
      }

      for (int k = 0; k < m_exp.nProtons; ++k) {
        if (!stateContributes(i, j, k)) {
          continue;
        }

        // set cell to dot product of H_B and I
        nZeeman(i, j) += m_staticBField.dot(spinVector(i, j, k));
      }
    }
  }

  nZeeman *= (-1.0*g_1H*NUC_MAGNETON);

  // DEBUG:
  // cout << nZeeman << endl;

  return nZeeman;
}

MatrixXcd SpinHamiltonian::hyperFine() const
{
  //Compute Hyperfine couplings matrix=========================================
  MatrixXcd hyperfine(m_exp.dimension, m_exp.dimension);
  hyperfine.setZero();

  for (int i = 0; i < m_exp.dimension; ++i) {
    for (int j = 0; j < m_exp.dimension; ++j) {
      //compute elements of s vector
      const Vector3cd s = spinVector(i, j, m_exp.nProtons);

      for (int k = 0; k < m_exp.nProtons; ++k) {    //loop over nuclei
        if (!stateContributes(i, j, k)) {
          continue;
        }

        //multiply atensor by I
        //multiply s by atensor_I
        hyperfine(i, j) += s.dot(m_exp.aTensor * spinVector(i, j, k));
        ///TODO: before: it was this, BUT only conjugate because of zdotu instead of zdotc in original - OR?
        ///hyperfine(i, j) += (m_exp.aTensor.transpose() * spinVector(i, j, k)).conjugate().dot(s);
      }
    }
  }
  hyperfine *= h * 1.0E6;

  // DEBUG:
  // cout << hyperfine << endl;

  return hyperfine;
}

MatrixXcd SpinHamiltonian::electronZeeman() const
{
  //Compute eZeeman============================================================  
  MatrixXcd eZeeman(m_exp.dimension, m_exp.dimension);
  eZeeman.setZero();
  //first multiply the g tensor with the static magnetic field hamiltonian
  const Vector3cd gDotH_B = m_exp.gTensor * m_staticBField;

  //depending on the convention, i might have to tranpose the gtensor here
  for (int i = 0; i < m_exp.dimension; ++i) {
    for (int j = 0; j < m_exp.dimension; ++j) {
      //m_exp.nProtons is always the index of the electron spin
      if (!stateContributes(i, j, m_exp.nProtons)) {
        continue;
      }

      eZeeman(i,j) = gDotH_B.dot(spinVector(i, j, m_exp.nProtons));
    }
  }
  eZeeman *= Bohrm;

  // DEBUG
  // cout << eZeeman << endl;

  return eZeeman;
}

MatrixXcd SpinHamiltonian::magneticMoments() const
{
  MatrixXcd moments(m_exp.dimension, m_exp.dimension);

  for (int i = 0; i < m_exp.dimension; ++i) {
    for (int j = 0; j < m_exp.dimension; ++j) {
      //m_exp.nProtons is always the index of the electronic spin state
      moments(i, j) = magneticMoment(i, j);

//       cout << i << '\t' << j << '\t' << "FINAL:" << '\t' << moments(i, j) << endl;
    }
  }

  //cout << moments << endl;
  return moments;
}

inline c_double SpinHamiltonian::magneticMoment(const int i, const int j) const
{
  c_double ret = 0;
  for (int k = 0; k < m_exp.nProtons+1; ++k) {
    if (!stateContributes(i, j, k, IncludeElectron)) {
      continue;
    }

    const int a = spinState(i, k);  //spin state of state k in row i
    const int b = spinState(j, k);  //spin state of state k in column j
    c_double xMoment = PauliMatrix::X(a, b);

    if (k != m_exp.nProtons) {
      xMoment *= -1.0 * g_1H * NUC_MAGNETON;
    } else {
      xMoment *= 2.023 * Bohrm;
    }

//         cout << i << '\t' << j << '\t' << k << '\t' << xMoment << endl;
    ret += xMoment;
  }
  return ret;
}

///TODO: optimize!
MatrixXd SpinHamiltonian::intensityMatrix(const MatrixXcd& eigenVectors) const {
  const MatrixXcd moments = magneticMoments();
  ///TODO: take direction of B0 and B1 into account, integrate over plane
  MatrixXd intensities(m_exp.dimension, m_exp.dimension);
  for(int i = 0; i < m_exp.dimension; ++i) {
    /// right part: M | Psi_i >
    const MatrixXcd mTimesPsiI = moments * eigenVectors.col(i);
    for(int j = 0; j < m_exp.dimension; ++j) {
      if (i != j) {
        /// left part with < Psi_j | and abs squared: |< Psi_j | M | Psi_i >|^2
        intensities(i, j) = (eigenVectors.col(j).adjoint() * mTimesPsiI).norm();
        ///TODO: eq 3-24, p 52 says: |< j|M|i > dot H_1|^2
        ///meaning: what about H_1?
      } else {
        intensities(i, j) = 0;
      }
    }
  }
  return intensities;
}

void SpinHamiltonian::calculateIntensity(const double mwFreq, QTextStream* out) const
{
  //Diagonalize the total Hamiltonian matrix===================================
  SelfAdjointEigenSolver<MatrixXcd> eigenSolver(hamiltonian());
  const VectorXd eigenValues = eigenSolver.eigenvalues();
  const MatrixXcd eigenVectors = eigenSolver.eigenvectors();

  const MatrixXcd moments = magneticMoments();

  double intensity = 0;
  for (int i = 0;i < m_exp.dimension; ++i) {
    for (int j = i + 1; j < m_exp.dimension; ++j) {
        // transition frequency:
      const double freq = (1.0/h/1.0E9 * abs(eigenValues(i) - eigenValues(j)));
      // assume it's only seen when energy is below frequency threshold
      if (abs(mwFreq/freq - 1.0) > 5.0E-4) {
        continue;
      }
      intensity += (eigenVectors.col(j).adjoint() * moments * eigenVectors.col(i)).norm();
    }
  }
  (*out) << scientific << m_B << '\t' << (intensity * 2.0 * M_PI * (Bohrm / hbar) * (Bohrm / hbar)) << endl;
}

void SpinHamiltonian::calculateTransitions() const
{
  //Diagonalize the total Hamiltonian matrix===================================
  SelfAdjointEigenSolver<MatrixXcd> eigenSolver(hamiltonian());
  const VectorXd eigenValues = eigenSolver.eigenvalues();
  const MatrixXcd eigenVectors = eigenSolver.eigenvectors();

  MatrixXd probabilities = intensityMatrix(eigenVectors);
  probabilities /= probabilities.maxCoeff();
  
  cout << "\n\nAllowed Transitions: \t\t\tB= " << fixed << m_B << endl;
  cout << "---------------------------------------------------\n";
  cout << "Transition\tFrequency (GHz)\t\tProbability" << endl;
  cout << "---------------------------------------------------\n";
  int transitions = 1;
  for (int i = 0;i < m_exp.dimension; ++i) {
    for (int j = i + 1; j < m_exp.dimension; ++j) {
      const double probability = probabilities(i, j);
      if (probability > 1.0E-6) {
        cout.width(10);
        cout << right << transitions++ << '\t';
        cout.width(14);
        cout.precision(5);
        // transition frequency:
        cout << right << (1.0/h/1.0E9 * abs(eigenValues(i) - eigenValues(j))) << "\t\t";
        cout.precision(8);
        cout << probability << endl;
      }
    }
  }
  cout << "---------------------------------------------------\n";
  cout << transitions << " transitions in total" << endl;
}

#endif // MW_BACHELOR_SPINHAMILTONIAN_H
