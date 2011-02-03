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

#include "types.h"
#include "constants.h"
#include "experiment.h"

using namespace Constants;

// Pauli Matrices
namespace PauliMatrix {
const Matrix2c X = (Matrix2c() << 0, 0.5, 0.5, 0).finished();
const Matrix2c Y = (Matrix2c() << 0, c_fp(0, -0.5), c_fp(0, 0.5), 0).finished();
const Matrix2c Z = (Matrix2c() << 0.5, 0, 0, -0.5).finished();
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
    SpinHamiltonian(const fp B, const Experiment& experiment);
    ~SpinHamiltonian();

    /// list all possible transitions with their frequency for the given B field
    void calculateTransitions() const;
    /// calculate total intensity of all transitions that are valid for the
    /// incoming microwave frequency
    void calculateIntensity(const fp mwFreq, QTextStream* out) const;

  private:
    /// the complete hamiltonian
    MatrixXc hamiltonian() const;
    /// nuclear Zeeman Hamiltonian component
    MatrixXc nuclearZeeman() const;
    /// hyper fine Hamiltonian component
    MatrixXc hyperFine() const;
    /// electron Zeeman Hamiltonian component
    MatrixXc electronZeeman() const;

    /// interprets @p i as binary number and returns the k-th bit of it
    inline bool spinState(int i, int k) const;

    /// return spin vector from pauli matrices
    inline Vector3c spinVector(int i, int j, int k) const;

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
    MatrixXc magneticMoments() const;

    inline c_fp magneticMoment(const int i, const int j) const;

    /// intensitiy matrix with coefficients (i, j) = |< psi_j | M | psi_i>|^2
    /// psi_i being the i-th eigen vector
    /// M being the magnetic moment matrix
    MatrixX intensityMatrix(const MatrixXc& eigenVectors) const;

    const fp m_B;
    const Experiment m_exp;
    Vector3c m_staticBField;
};

SpinHamiltonian::SpinHamiltonian(const fp B, const Experiment& experiment)
: m_B(B)
, m_exp(experiment)
, m_staticBField(m_exp.staticBField(B))
{
}

SpinHamiltonian::~SpinHamiltonian()
{
}

Vector3c SpinHamiltonian::spinVector(int i, int j, int k) const
{
  const int a = spinState(i, k); //spin state of state k in row i
  const int b = spinState(j, k); //spin state of state k in column j
  return (Vector3c() << PauliMatrix::X(a, b), PauliMatrix::Y(a, b), PauliMatrix::Z(a, b)).finished();
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

MatrixXc SpinHamiltonian::hamiltonian() const
{
  return nuclearZeeman() + hyperFine() + electronZeeman();
}

MatrixXc SpinHamiltonian::nuclearZeeman() const
{
  //Compute nZeeman============================================================  
  MatrixXc nZeeman(m_exp.dimension, m_exp.dimension);
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

MatrixXc SpinHamiltonian::hyperFine() const
{
  //Compute Hyperfine couplings matrix=========================================
  MatrixXc hyperfine(m_exp.dimension, m_exp.dimension);
  hyperfine.setZero();

  for (int i = 0; i < m_exp.dimension; ++i) {
    for (int j = 0; j < m_exp.dimension; ++j) {
      //compute elements of s vector
      const Vector3c s = spinVector(i, j, m_exp.nProtons);

      for (int k = 0; k < m_exp.nProtons; ++k) {    //loop over nuclei
        if (!stateContributes(i, j, k)) {
          continue;
        }

        //multiply atensor by I
        //multiply s by atensor_I
        hyperfine(i, j) += s.dot(m_exp.aTensor * spinVector(i, j, k));
      }
    }
  }
  hyperfine *= h * 1.0E6;

  // DEBUG:
  // cout << hyperfine << endl;

  return hyperfine;
}

MatrixXc SpinHamiltonian::electronZeeman() const
{
  //Compute eZeeman============================================================  
  MatrixXc eZeeman(m_exp.dimension, m_exp.dimension);
  eZeeman.setZero();
  //first multiply the g tensor with the static magnetic field hamiltonian
  const Vector3c gDotH_B = m_exp.gTensor * m_staticBField;

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

MatrixXc SpinHamiltonian::magneticMoments() const
{
  MatrixXc moments(m_exp.dimension, m_exp.dimension);

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

c_fp SpinHamiltonian::magneticMoment(const int i, const int j) const
{
  c_fp ret = 0;
  for (int k = 0; k < m_exp.nProtons+1; ++k) {
    if (!stateContributes(i, j, k, IncludeElectron)) {
      continue;
    }

    const int a = spinState(i, k);  //spin state of state k in row i
    const int b = spinState(j, k);  //spin state of state k in column j
    c_fp xMoment = PauliMatrix::X(a, b);

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
MatrixX SpinHamiltonian::intensityMatrix(const MatrixXc& eigenVectors) const {
  const MatrixXc moments = magneticMoments();
  ///TODO: take direction of B0 and B1 into account, integrate over plane
  MatrixX intensities(m_exp.dimension, m_exp.dimension);
  for(int i = 0; i < m_exp.dimension; ++i) {
    /// right part: M | Psi_i >
    const MatrixXc mTimesPsiI = moments * eigenVectors.col(i);
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

void SpinHamiltonian::calculateIntensity(const fp mwFreq, QTextStream* out) const
{
  //Diagonalize the total Hamiltonian matrix===================================
  SelfAdjointEigenSolver<MatrixXc> eigenSolver(hamiltonian());
  const VectorX eigenValues = eigenSolver.eigenvalues();
  const MatrixXc eigenVectors = eigenSolver.eigenvectors();

  const MatrixXc moments = magneticMoments();

  fp intensity = 0;
  for (int i = 0;i < m_exp.dimension; ++i) {
    for (int j = i + 1; j < m_exp.dimension; ++j) {
        // transition frequency:
      const fp freq = (1.0/h/1.0E9 * abs(eigenValues(i) - eigenValues(j)));
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
  SelfAdjointEigenSolver<MatrixXc> eigenSolver(hamiltonian());
  const VectorX eigenValues = eigenSolver.eigenvalues();
  const MatrixXc eigenVectors = eigenSolver.eigenvectors();

  MatrixX probabilities = intensityMatrix(eigenVectors);
  probabilities /= probabilities.maxCoeff();
  
  cout << "\n\nAllowed Transitions: \t\t\tB= " << fixed << m_B << endl;
  cout << "---------------------------------------------------\n";
  cout << "Transition\tFrequency (GHz)\t\tProbability" << endl;
  cout << "---------------------------------------------------\n";
  int transitions = 1;
  for (int i = 0;i < m_exp.dimension; ++i) {
    for (int j = i + 1; j < m_exp.dimension; ++j) {
      const fp probability = probabilities(i, j);
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
