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

#include "experiment.h"

#include "constants.h"
#include "nucleus.h"
#include "spins.h"

#include <algorithm>

#include <boost/foreach.hpp>

using namespace std;

int dimensionForNuclei(const vector<Nucleus>& nuclei)
{
  int dimension = 1;
  dimension *= 2; // electron
  BOOST_FOREACH(const Nucleus& nucleus, nuclei) {
    dimension *= (nucleus.twoJ + 1);
  }
  return dimension;
}

bool sortNuclei_cmp(const Nucleus& l, const Nucleus& r)
{
  return l.twoJ <= r.twoJ;
}

vector<Nucleus> sortNuclei(vector<Nucleus> nuclei)
{
  sort(nuclei.begin(), nuclei.end(), sortNuclei_cmp);
  return nuclei;
}

Experiment::Experiment(const vector<Nucleus>& nuclei_)
: nuclei(sortNuclei(nuclei_))
, dimension(dimensionForNuclei(nuclei))
, mwFreqGHz(0)
{
  setGTensor(Matrix3::Identity() * Constants::g_E);
}

Experiment Experiment::generateDummy(int protons, int nitrogens)
{
  vector<Nucleus> nuclei;
  const Matrix3c A = Matrix3c::Identity() * 1420;
  for(int i = 0; i < protons; ++i) {
    nuclei.push_back(Nucleus("1H", 1, 1, A, Constants::g_1H));
  }
  for(int i = 0; i < nitrogens; ++i) {
    nuclei.push_back(Nucleus("14N", 2, 14, A, Constants::g_14N));
  }
  return Experiment(nuclei);
}

Spins Experiment::spinSystem() const
{
  int spinHalfs = 1; // single unbound electron
  int spinOnes = 0;
  BOOST_FOREACH(const Nucleus& nucleus, nuclei) {
    switch(nucleus.twoJ) {
      case 1:
        spinHalfs++;
        break;
      case 2:
        spinOnes ++;
        break;
      default:
        cerr << "unhandled J = " << (nucleus.twoJ/2) << " for nucleus " << nucleus.name << endl;
        throw "unhandled J";
        break;
    }
  }
  Spins s(spinHalfs, spinOnes);
  if (s.states != dimension) {
    cerr << "spin states and experiment dimension differ" << endl;
    throw "spin system initialization error";
  }
  return s;
}

Vector3c Experiment::staticBField(const fp B) const
{
  if (getenv("USE_LABOR_Z")) {
    return (Vector3c() << 0, 0, B).finished();
  }
  Vector3c gz;
  const Vector3& E = m_gTensorEigenValues;
  const Matrix3c& V = m_gTensorEigenVectors;
  // |g_z> has the largest eigen value
  if (E(0) > E(1)) {
    if (E(0) > E(2)) {
      gz = V.col(0);
    } else {
      gz = V.col(2);
    }
  } else if (E(1) > E(2)) {
    gz = V.col(1);
  } else {
    gz = V.col(2);
  }
  return gz * B / gz.norm();
}

void Experiment::setGTensor(const Matrix3& gTensor)
{
  m_gTensor = gTensor;
  Eigen::EigenSolver<Matrix3> solver(gTensor);
  m_gTensorEigenValues = solver.eigenvalues().real();
  m_gTensorEigenVectors = solver.eigenvectors();
}

const Matrix3c& Experiment::gTensorEigenVectors() const
{
  return m_gTensorEigenVectors;
}

const Vector3& Experiment::gTensorEigenValues() const
{
  return m_gTensorEigenValues;
}
