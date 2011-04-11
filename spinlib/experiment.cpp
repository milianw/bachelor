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
#include <map>

#include <boost/foreach.hpp>

using namespace std;
using namespace Eigen;

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
  const Matrix3 A = Matrix3::Identity() * 1420;
  for(int i = 0; i < protons; ++i) {
    nuclei.push_back(Nucleus("1H", 1, 1, A, Constants::g_1H, Vector3::Zero()));
  }
  for(int i = 0; i < nitrogens; ++i) {
    nuclei.push_back(Nucleus("14N", 2, 14, A, Constants::g_14N, Vector3::Zero()));
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
  const Vector3c& gz = m_gTensorEigenVectors.col(2);
  return gz * B / gz.norm();
}

void Experiment::setGTensor(const Matrix3& gTensor)
{
  m_gTensor = gTensor;

  EigenSolver<Matrix3> solver(gTensor);

  multimap<double, int> sortedKeys;
  sortedKeys.insert(pair<double, int>(solver.eigenvalues()(0).real(), 0));
  sortedKeys.insert(pair<double, int>(solver.eigenvalues()(1).real(), 1));
  sortedKeys.insert(pair<double, int>(solver.eigenvalues()(2).real(), 2));

  multimap<double, int>::iterator it = sortedKeys.begin();
  m_gTensorEigenVectors.col(0) = solver.eigenvectors().col(it->second);
  m_gTensorEigenValues(0) = it->first;
  ++it;
  m_gTensorEigenVectors.col(1) = solver.eigenvectors().col(it->second);
  m_gTensorEigenValues(1) = it->first;
  ++it;
  m_gTensorEigenVectors.col(2) = solver.eigenvectors().col(it->second);
  m_gTensorEigenValues(2) = it->first;

}

const Matrix3c& Experiment::gTensorEigenVectors() const
{
  return m_gTensorEigenVectors;
}

const Vector3& Experiment::gTensorEigenValues() const
{
  return m_gTensorEigenValues;
}
