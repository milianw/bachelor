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

#include "spinhamiltonian.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>

#include "constants.h"
#include "experiment.h"
#include "spins.h"
#include "nucleus.h"

using namespace Constants;
using namespace Eigen;
using namespace std;

SpinHamiltonian::SpinHamiltonian(const fp B, const Experiment& experiment)
: m_B(B)
, m_exp(experiment)
, m_spins(m_exp.spinSystem())
, m_staticBField(m_exp.staticBField(B))
{
}

SpinHamiltonian::~SpinHamiltonian()
{
}

inline Vector3c SpinHamiltonian::spinVector(int bra, int ket, int k) const
{
  return m_spins.spinVector(bra, ket, k);
}

inline int SpinHamiltonian::spinState(int state, int k) const
{
  return m_spins.spinInState(k, state);
}

inline bool SpinHamiltonian::stateContributes(int bra, int ket, int k, bool ignoreElectron) const
{
  for(int n = 0; n < m_spins.elements; ++n) {
    if (n == k || (ignoreElectron && n == 0)) {
      continue;
    } else if (m_spins.spinInState(n, bra) != m_spins.spinInState(n, ket)) {
      return false;
    }
  }
  return true;
}

MatrixXc SpinHamiltonian::hamiltonian() const
{
  MatrixXc H(m_spins.states, m_spins.states);
  H.setZero();
  addNuclearZeeman(H);
  addHyperFine(H);
  addElectronZeeman(H);
  addQuadrupole(H);
  return H;
}

void SpinHamiltonian::addNuclearZeeman(MatrixXc& H) const
{
  for (int bra = 0; bra < m_spins.states; ++bra) {
    for (int ket = 0; ket < m_spins.states; ++ket) {
      if (spinState(bra, 0 /* = Electron */) != spinState(ket, 0 /* = Electron */)) {
        continue;  //matrix elements between different electron states are zero
      }

      // k = 1 to skip electron
      c_fp colVal = 0;
      for (int k = 1; k < m_spins.elements; ++k) {
        if (!stateContributes(bra, ket, k)) {
          continue;
        }

        // set cell to dot product of H_B and I
        ///TODO: .dot i.e. with adjoint or cwiseProduct (should not be a difference since B field is never complex?)
        c_fp val = m_staticBField.dot(spinVector(bra, ket, k));
        if (k < m_spins.spinHalfs) {
          // J = 1/2
          val *= g_1H;
        } else {
          // J = 1
          val *= g_14N;
        }
        colVal += val;
      }
      H(bra, ket) += colVal * -1.0 * NUC_MAGNETON;
    }
  }
}

void SpinHamiltonian::addHyperFine(MatrixXc& H) const
{
  for (int bra = 0; bra < m_spins.states; ++bra) {
    for (int ket = 0; ket < m_spins.states; ++ket) {
      //compute elements of s vector
      const Vector3c s = spinVector(bra, ket, 0 /* = Electron */);

      c_fp colVal = 0;
      // k = 1 to skip electron
      for (int k = 1; k < m_spins.elements; ++k) {    //loop over nuclei
        if (!stateContributes(bra, ket, k)) {
          continue;
        }

        //multiply atensor by I
        //multiply s by atensor_I
        ///TODO: proper aTensor for different nuclei
        ///NOTE: don't use .dot() as that would do s.adjoint() * ...
        ///      but for operators this is wrong
        colVal += s.cwiseProduct(m_exp.nuclei.at(k - 1).A * spinVector(bra, ket, k)).sum();
      }
      H(bra, ket) += colVal * h * 1.0E6;
    }
  }
}

// equation: \beta S * g * H
// \beta: borh magneton
// S: electron Spin Operator
// g: g Tensor
// H: static B Field hamiltonian
void SpinHamiltonian::addElectronZeeman(MatrixXc& H) const
{
  //first multiply the g tensor with the static magnetic field hamiltonian
  const Vector3c gDotH_B = m_exp.gTensor() * m_staticBField;
  //depending on the convention, i might have to tranpose the gtensor here
  for (int bra = 0; bra < m_spins.states; ++bra) {
    for (int ket = 0; ket < m_spins.states; ++ket) {
      // 0 is always the index of the electron spin
      if (!stateContributes(bra, ket, 0 /* = Electron */)) {
        continue;
      }

      ///NOTE: don't use .dot() as that would do s.adjoint() * ...
      ///      but for operators this is wrong
      H(bra, ket) += spinVector(bra, ket, 0 /* = Electron */).cwiseProduct(gDotH_B).sum() * Bohrm;
    }
  }
}

void SpinHamiltonian::addQuadrupole(MatrixXc& H) const
{

  // k = 1 to skip electron
  for (int k = 1; k < m_spins.elements; ++k) {    //loop over nuclei
    const Nucleus& nucleus = m_exp.nuclei.at(k - 1);
    if (nucleus.twoJ < 1) {
      // only spin one and higher has quadrupole interaction
      continue;
    } else if (nucleus.Q.isZero()) {
      continue;
    }

    ///NOTE: must be extended for J > 1 if ever supported
    const Matrix3c Isquared[3] = {PauliMatrix_J_one::X * PauliMatrix_J_one::X,
                                  PauliMatrix_J_one::Y * PauliMatrix_J_one::Y,
                                  PauliMatrix_J_one::Z * PauliMatrix_J_one::Z};

    for (int bra = 0; bra < m_spins.states; ++bra) {
      for (int ket = 0; ket < m_spins.states; ++ket) {
        if (!stateContributes(bra, ket, k, IncludeElectron)) {
          continue;
        }
        c_fp colVal = 0;

        int spinBra = spinState(bra, k);
        int spinKet = spinState(ket, k);

        // iterate over x,y,z
        ///FIXME: must be rewritten once Q is rotated
        for(int i = 0; i < 3; ++i) {
          colVal += nucleus.Q(i) * Isquared[i](spinBra, spinKet);
        }

        // Q is in MHz, convert to MKS
        H(bra, ket) += colVal * h * 1.0E6;
      }
    }
  }
}

MatrixXc SpinHamiltonian::magneticMoments() const
{
  MatrixXc moments(m_spins.states, m_spins.states);

  const Vector3c g_x = m_exp.gTensorEigenVectors().col(0) / m_exp.gTensorEigenVectors().col(0).norm();
  const fp g = m_exp.gTensorEigenValues()(0);

  for (int bra = 0; bra < m_spins.states; ++bra) {
    for (int ket = 0; ket < m_spins.states; ++ket) {
      c_fp moment = 0;

      // sum x-moments of nuclei and electron in system
      for (int k = 0; k < m_spins.elements; ++k) {
        if (!stateContributes(bra, ket, k, IncludeElectron)) {
          continue;
        }

        c_fp xMoment = g_x.dot(spinVector(bra, ket, k));

        if (k == 0) {
          // electron
          xMoment *= g * Bohrm;
        } else if (k < m_spins.spinHalfs) {
          // J = 1/2 nucleus
          xMoment *= -1.0 * g_1H * NUC_MAGNETON;
        } else {
          // J = 1 nucleus
          xMoment *= -1.0 * g_14N * NUC_MAGNETON;
        }

        moment += xMoment;
      }

      moments(bra, ket) = moment;
    }
  }

  return moments;
}

MatrixX SpinHamiltonian::intensityMatrix(const MatrixXc& eigenVectors) const {
  ///TODO: take direction of B0 and B1 into account, integrate over plane
  return (eigenVectors.adjoint() * magneticMoments() * eigenVectors).cwiseAbs2();
}

fp SpinHamiltonian::calculateIntensity() const
{
  //Diagonalize the total Hamiltonian matrix===================================
  SelfAdjointEigenSolver<MatrixXc> eigenSolver(hamiltonian());
  const VectorX& eigenValues = eigenSolver.eigenvalues();
  const MatrixXc& eigenVectors = eigenSolver.eigenvectors();

  ///NOTE: intensityMatrix has even higher peak mem consumption due to
  ///      matrix product allocations and temporary
  const MatrixXc moments = magneticMoments();

  const char* thresholdStr = getenv("FREQUENCY_THRESHOLD");
  float threshold = 1.0E-3;
  if (thresholdStr) {
    threshold = atof(thresholdStr);
  }

  fp intensity = 0;
  ///TODO: take direction of B0 and B1 into account, integrate over plane
  for (int i = 0;i < m_exp.dimension; ++i) {
    ///TODO: lower half - results x2?
    for (int j = i + 1; j < m_exp.dimension; ++j) {
      // transition frequency:
      const fp freq = (1.0/h/1.0E9 * abs(eigenValues(i) - eigenValues(j)));
      // assume it's only seen when energy is below frequency threshold
      if (abs(freq - m_exp.mwFreqGHz) > threshold * m_exp.mwFreqGHz) {
        continue;
      }
      /// < Psi_j | and abs squared: |< Psi_j | M | Psi_i >|^2 M | Psi_i >
      intensity += (eigenVectors.col(j).adjoint() * moments * eigenVectors.col(i)).squaredNorm();
      ///TODO: eq 3-24, p 52 says: |< j|M|i > dot H_1|^2
      ///meaning: what about H_1?
    }
  }
  return (intensity * 2.0 * M_PI * (Bohrm / hbar) * (Bohrm / hbar));
}

void SpinHamiltonian::calculateTransitions() const
{
  //Diagonalize the total Hamiltonian matrix===================================
  SelfAdjointEigenSolver<MatrixXc> eigenSolver(hamiltonian());
  const VectorX eigenValues = eigenSolver.eigenvalues();
  const MatrixXc eigenVectors = eigenSolver.eigenvectors();
  MatrixX probabilities = intensityMatrix(eigenVectors);
  probabilities /= probabilities.maxCoeff();

  float threshold = 1.0E-6;
  if (const char* thresholdStr = getenv("PROBABILITY_THRESHOLD")) {
    threshold = atof(thresholdStr);
  }

  multimap<fp, fp> intensities;
  for (int i = 0; i < m_exp.dimension; ++i) {
    for (int j = i + 1; j < m_exp.dimension; ++j) {
      const fp probability = probabilities(i, j);
      if (probability > threshold) {
        intensities.insert(pair<fp, fp>((1.0/h/1.0E9 * abs(eigenValues(i) - eigenValues(j))), probability));
      }
    }
  }

  cout << "\n\nAllowed Transitions: \t\t\tB= " << fixed << m_B << endl;
  cout << "---------------------------------------------------\n";
  cout << "Transition\tFrequency (GHz)\t\tProbability" << endl;
  cout << "---------------------------------------------------\n";
  multimap< fp, fp >::const_iterator it = intensities.begin();
  multimap< fp, fp >::const_iterator end = intensities.end();
  int transitions = 0;
  while(it != end) {
    cout.width(10);
    cout << right << ++transitions << '\t';
    cout.width(14);
    cout.precision(5);
    cout << right << it->first << "\t\t";
    cout.precision(8);
    cout << it->second << endl;
    ++it;
  }
  cout << "---------------------------------------------------\n";
  cout << transitions << " transitions in total" << endl;
}

VectorX SpinHamiltonian::calculateEigenValues() const
{
  //Diagonalize the total Hamiltonian matrix===================================
  SelfAdjointEigenSolver<MatrixXc> eigenSolver(hamiltonian(), EigenvaluesOnly);
  return eigenSolver.eigenvalues();
}
