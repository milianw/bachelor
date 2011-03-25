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

#include "resonancefield.h"

#include "experiment.h"
#include "spinhamiltonian.h"
#include "constants.h"
#include "nucleus.h"

#include <stack>
#include <iostream>
#include <cmath>

#include <boost/foreach.hpp>

#define GNUPLOT_DEBUG(x)
#define DEBUG(x)

namespace {

// segment of B
struct Segment {
  Segment()
  : B_min(0), B_max(0)
  {
  }

  Segment(fp min, fp max)
  : B_min(min), B_max(max)
  {
  }

  bool operator==(const Segment& other) const
  {
    return B_min == other.B_min && B_max == other.B_max;
  }

  fp B_min;
  fp B_max;
};

static BisectAnswer notResonantAnswer(const fp from, const fp to)
{
  BisectAnswer answer;
  answer.status = BisectAnswer::NotResonant;
  answer.from = from;
  answer.to = to;
  return answer;
}

static BisectAnswer continueAnswer(const fp from, const fp to, const BisectNode& mid)
{
  BisectAnswer answer;
  answer.status = BisectAnswer::Continue;
  answer.from = from;
  answer.to = to;
  answer.mid = mid;
  return answer;
}

static BisectAnswer resonantAnswer(const fp from, const fp to, const BisectNode& mid)
{
  BisectAnswer answer;
  answer.status = BisectAnswer::Resonant;
  answer.from = from;
  answer.to = to;
  answer.mid = mid;
  return answer;
}

}

using namespace Constants;
using namespace Eigen;
using namespace std;

ResonanceField::ResonanceField(const Experiment& exp)
: m_exp(exp)
  // to atomic units:
, m_mwFreq(exp.mwFreqGHz * 1.0E9 * h)
, m_loopingResonanceCanOccur(checkForLoopingResonance())
, m_lambda(calculateLambda())
{
  DEBUG(cout << (m_loopingResonanceCanOccur ? "looping resonance can occur" : "no looping resonance can occur") << endl;)
}

bool ResonanceField::checkForLoopingResonance() const
{
  ///TODO: make sure they are always sorted in the correct order...
  VectorX eVals = SpinHamiltonian(0, m_exp).calculateEigenValues();
  DEBUG(cout << "maximum splitting at zero field: " << (eVals(m_exp.dimension - 1) - eVals(0))/1E6/h << "MHz" << endl;)
  return (eVals(m_exp.dimension - 1) - eVals(0)) >= m_mwFreq;
}

fp ResonanceField::calculateLambda() const
{
  //implementation of eq 14
  ///FIXME: generalize this once we support multiple electrons or nuclei with I != 1/2
  fp lambda = 0;

  ///TODO: is this correct?
  Vector3c n = m_exp.staticBFieldDirection / m_exp.staticBFieldDirection.norm();

  // we assume only a single electron
  // S = 0.5
  lambda += Bohrm * 0.5 * (n.transpose() * m_exp.gTensor).norm();

  BOOST_FOREACH(const Nucleus& nucleus, m_exp.nuclei) {
    lambda += NUC_MAGNETON * nucleus.g * (nucleus.twoJ / 2);
  }
  return lambda;
}

BisectNode ResonanceField::diagonalizeNode(const fp B) const
{
  SpinHamiltonian H(B, m_exp);
  SelfAdjointEigenSolver<MatrixXc> eigenSolver(H.hamiltonian());
  const VectorX& E = eigenSolver.eigenvalues();
  const MatrixXc& eigenVectors = eigenSolver.eigenvectors();

  MatrixXc G(m_exp.dimension, m_exp.dimension);
  G.setZero();
  H.addNuclearZeeman(G);
  H.addElectronZeeman(G);

  VectorX E_deriv(m_exp.dimension);
  for(int u = 0; u < m_exp.dimension; ++u) {
    // <u| G |u> => expectation value is always real
    E_deriv(u) = (eigenVectors.col(u).adjoint() * G * eigenVectors.col(u))(0, 0).real();
  }
  return BisectNode(B, E, E_deriv);
}

void ResonanceField::cleanupResonancyField(vector< fp >& field)
{
  if (field.size() < 2) {
    return;
  }

  sort(field.begin(), field.end());
  vector<fp>::iterator it = field.begin();
  while(it != field.end() - 1) {
    if (abs(*it / *(it+1) - 1) < 0.001) {
      field.erase(it + 1);
    } else {
      ++it;
    }
  }
}

vector< fp > ResonanceField::calculate(fp B_min, fp B_max)
{
  m_eVals.clear();
  vector< fp > field = findRoots(resonantSegments(B_min, B_max));

  if (field.empty()) {
    cerr << "ATTENTION: no resonant segments found in range [" << B_min << ", " << B_max << "] for mwFreq = " << m_exp.mwFreqGHz << endl;
    return field;
  }

  cleanupResonancyField(field);

  DEBUG(cout << "cleaned resonance field contains " << field.size() << " roots" << endl;)

  m_eVals.clear();
  return field;
}

BisectAnswer ResonanceField::checkSegment(const BisectNode& from, const BisectNode& to) const
{
  const fp B_diff = to.B - from.B;
  bool resonancePossible = false;
  if ((to.E(m_exp.dimension - 1) - to.E(0)) > m_mwFreq) {
    if (!m_loopingResonanceCanOccur) {
      // for all kombinations u,v do eq 13
      for(int u = 0; u < m_exp.dimension; ++u) {
        for(int v = u + 1; v < m_exp.dimension; ++v) {
          // R_{uv}(B_q) * R_{uv}(B_r) <= 0
          if (((from.E(v) - from.E(u) - m_mwFreq) * (to.E(v) - to.E(u) - m_mwFreq)) <= 0) {
            resonancePossible = true;
            break;
          }
        }
        if (resonancePossible) {
          break;
        }
      }
    } else {
      // else for all kombinations u,v do eq 15
      for(int u = 0; u < m_exp.dimension; ++u) {
        for(int v = u + 1; v < m_exp.dimension; ++v) {
          if (abs((from.E(v) - from.E(u) + to.E(v) - to.E(u)) * 0.5 - m_mwFreq) <= m_lambda * B_diff) {
            resonancePossible = true;
            break;
          }
        }
        if (resonancePossible) {
          break;
        }
      }
    }
  } // else ResonPossible = false
  if (resonancePossible) {
    // eq 11
    const BisectNode mid = diagonalizeNode((from.B + to.B) * 0.5);
    fp epsilon = 0;
    for(int u = 0; u < m_exp.dimension; ++u) {
      const fp E_u_tilde = 0.5 * (from.E(u) + to.E(u)) + B_diff / 8.0 * (from.E_deriv(u) - to.E_deriv(u));
      fp epsilon_u = abs(mid.E(u) - E_u_tilde);
      if (epsilon_u > epsilon) {
        epsilon = epsilon_u;
      }
    }
    epsilon *= 2;

    static const float threshold = getenv("RESFIELD_THRESHOLD") ? atof(getenv("RESFIELD_THRESHOLD")) : 1.0E-4;

    if (epsilon > threshold * m_mwFreq) {
      return continueAnswer(from.B, to.B, mid);
    } else {
      return resonantAnswer(from.B, to.B, mid);
    }
  } else {
    return notResonantAnswer(from.B, to.B);
  }
}

map< fp, fp > ResonanceField::resonantSegments(fp B_minStart, fp B_maxStart)
{
  //STEP 1: find knots
  map<fp, fp> resonantSegments;

  m_eVals[B_minStart] = diagonalizeNode(B_minStart);
  m_eVals[B_maxStart] = diagonalizeNode(B_maxStart);

  stack<Segment> segments;
  segments.push(Segment(B_minStart, B_maxStart));

  while (!segments.empty()) {
    const Segment s = segments.top();
    segments.pop();
    const BisectAnswer answer = checkSegment(m_eVals.at(s.B_min), m_eVals.at(s.B_max));
    switch(answer.status) {
      case BisectAnswer::NotResonant:
        // nothing to do
        break;
      case BisectAnswer::Continue:
        m_eVals[answer.mid.B] = answer.mid;
        segments.push(Segment(answer.from, answer.mid.B));
        segments.push(Segment(answer.mid.B, answer.to));
        break;
      case BisectAnswer::Resonant:
        m_eVals[answer.mid.B] = answer.mid;
        resonantSegments[answer.from] = answer.mid.B;
        resonantSegments[answer.mid.B] = answer.to;
        break;
    }
  }

  if (resonantSegments.empty()) {
    cerr << "ATTENTION: no resonant segments found in range [" << B_minStart << ", " << B_maxStart << "] for mwFreq = " << (m_exp.mwFreqGHz) << endl;
  }

  DEBUG(cout << "segmentation finished, " << resonantSegments.size() << " segments" << endl;)

  return resonantSegments;
}

fp evalPolynomial(const Vector4& p, const fp t, const fp mwFreq)
{
  return p(0) * t * t * t + p(1) * t * t + p(2) * t + p(3) - mwFreq;
}

fp evalDerivPolynomial(const Vector4& p, const fp t)
{
  return 3.0 * p(0) * t * t + 2.0 * p(1) * t + p(2);
}

fp newtonRaphson(const Vector4& p, const fp t, const fp mwFreq)
{
  return t - evalPolynomial(p, t, mwFreq) / evalDerivPolynomial(p, t);
}

vector<fp> ResonanceField::findRoots(const map<fp, fp>& resonantSegments)
{
  //STEP 2: find roots
  vector<fp> resonanceField;

  map<fp, fp>::const_iterator it = resonantSegments.begin();
  const map<fp, fp>::const_iterator end = resonantSegments.end();
  GNUPLOT_DEBUG(
    cout << "set xrange[" << it.key() * 0.9 << ":" << (end-1).value() * 1.1 << "]" << endl;
    cout << "u(x,min,max) = (x>=min)&&(x<max)? 1 : 1/0;" << endl;
  )
  while(it != end) {
    const vector<fp> roots = findRootsInSegment(m_eVals.at(it->first), m_eVals.at(it->second));
    copy(roots.begin(), roots.end(), back_inserter(resonanceField));
    ++it;
  }
  DEBUG(cout << "found " << resonanceField.size() << " roots in resonant segments" << endl;)
  return resonanceField;
}

// t_k from the wiki
fp trigonometricRoot(const int k, const fp p, const fp q)
{
    return 2.0 * sqrt(-p / 3.0) * cos(1.0/3.0 * acos(3.0 * q / (2.0 * p) * sqrt(-3.0 / p)) - k * 2.0 * M_PI / 3.0);
}

vector<fp> ResonanceField::findRootsInSegment(const BisectNode& from, const BisectNode& to) const
{
  vector<fp> roots;
  const fp B_diff = to.B - from.B;
  GNUPLOT_DEBUG(
    cout << "set arrow from " << B_min << ",-5e-26 to " << B_min << ",5e-26 ls 2" << endl;
    cout << "set arrow from " << B_max << ",-5e-26 to " << B_max << ",5e-26 ls 2" << endl;
  )
  // see eq 8
  static const Matrix4 M = (Matrix4() << 
                            2, -2, 1, 1,
                            -3, 3, -2, -1,
                            0, 0, 1, 0,
                            1, 0, 0, 0).finished();
  for(int u = 0; u < m_exp.dimension; ++u) {
    const Vector4 e_u = (Vector4() << from.E(u), to.E(u), B_diff * from.E_deriv(u), B_diff * to.E_deriv(u)).finished();
    for(int v = u + 1; v < m_exp.dimension; ++v) {
      const Vector4 e_v = (Vector4() << from.E(v), to.E(v), B_diff * from.E_deriv(v), B_diff * to.E_deriv(v)).finished();
      ///NOTE: paper has different notation: p(0) == p_3, p(1) == p_2, ...
      const Vector4 p = M * (e_v - e_u);
      fp root = 0;
      // apply newton raphson rootfinding algo if no looping resonance can occur or polynomial is monotonic
      if (!m_loopingResonanceCanOccur || (p(1) * p(1) - 3.0 * p(0) * p(2)) <= 0) {
        if (((from.E(v) - from.E(u) - m_mwFreq) * (to.E(v) - to.E(u) - m_mwFreq)) > 0) {
          continue;
        }
        // Newton-Raphson root finding
        fp t = newtonRaphson(p, 0.5, m_mwFreq); // x_0 = 0.5
        if (!isfinite(t)) {
          continue;
        }
        fp t2 = t;
        const int maxIt = 10000;
        int it = 0;
        do {
          t = t2;
          t2 = newtonRaphson(p, t, m_mwFreq);
        ///TODO: when to abort?
        } while(abs(t2/t - 1) > 1E-5 && (maxIt > ++it));
        root = from.B + t2 * B_diff;
      } else {
        ///FIXME
        cout << "EXPERIMENTAL - non-monotonic polynomial with three roots - CHECK RESULTS!" << endl;
        // http://en.wikipedia.org/wiki/Cubic_function#General_formula_of_roots
        // reduction to a monic trinomial
        const fp p_ = (3.0 * p(0) * p(2) - p(1) * p(1)) / (3.0 * p(0) * p(0));
        const fp q_ = (2.0 * p(1) * p(1) * p(1) - 9.0 * p(0) * p(1) * p(2) + 27.0 * p(0) * p(0) * p(3))
                        / (27.0 * p(0) * p(0) * p(0));
        // only real roots are interesting
        if (4.0 * p_ * p_ * p_ + 27.0 * q_ * q_ <= 0) {
            // trigonometric (and hyperbolic) method
            for(int i = 0; i < 3; ++i) {
                const fp t = trigonometricRoot(i, p_, q_);
                if (!isnan(t)) {
                    roots.push_back(t);
                }
            }
        }
        continue;
      }

      roots.push_back(root);
      GNUPLOT_DEBUG(
        cout << "replot u(x, " << B_min << "," << B_max << ") * ((((x-" << B_min << ")/" << B_diff << "))**3 * (" << p(0) << ") + (((x-" << B_min << ")/" << B_diff << "))**2 * (" << p(1) << ") + (((x-" << B_min << ")/" << B_diff << ")) * (" << p(2) << ") + (" << p(3) << " - " << mwFreq << ")) "
                  "title \"B_min = " << B_min << ", B_max = " << B_max << " || u = " << u << ", v = " << v << "\"" << endl;
        cout << "set arrow from " << root << ",-1e-24 to " << root << ",1e-24 ls 0" << endl;
        nodes << B_min << '\t' << (from.E(v) - from.E(u) - mwFreq) << endl;
        nodes << B_max << '\t' << (to.E(v) - to.E(u) - mwFreq) << endl;
      )
    }
  }
  return roots;
}
