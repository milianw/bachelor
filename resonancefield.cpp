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

#include <QtCore/QStack>
#include <QtCore/QDebug>

#define GNUPLOT_DEBUG(x)

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

}

using namespace Constants;

ResonanceField::ResonanceField(const Experiment& exp)
: m_exp(exp)
  // to atomic units:
, m_mwFreq(exp.mwFreqGHz * 1.0E9 * h)
, m_loopingResonanceCanOccur(checkForLoopingResonance())
, m_lambda(calculateLambda())
{
}

bool ResonanceField::checkForLoopingResonance() const
{
  ///TODO: make sure they are always sorted in the correct order...
  VectorX eVals = SpinHamiltonian(0, m_exp).calculateEigenValues();
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
  lambda += Bohrm * 0.5 * (n.transpose() * m_exp.gTensor).norm();

  // we assume n protons
  lambda += NUC_MAGNETON * m_exp.nProtons * g_1H * 0.4;
  return lambda;
}

BisectNode ResonanceField::diagonalizeNode(const fp B) const
{
  SpinHamiltonian H(B, m_exp);
  SelfAdjointEigenSolver<MatrixXc> eigenSolver(H.hamiltonian());
  VectorX E = eigenSolver.eigenvalues();
  const MatrixXc eigenVectors = eigenSolver.eigenvectors();
  const MatrixXc G = H.nuclearZeeman() + H.electronZeeman();
  VectorX E_deriv(m_exp.dimension);
  for(int u = 0; u < m_exp.dimension; ++u) {
    // <u| G |u> => expectation value is always real
    E_deriv(u) = (eigenVectors.col(u).adjoint() * G * eigenVectors.col(u))(0, 0).real();
  }
  return BisectNode(B, E, E_deriv);
}

QVector< fp > ResonanceField::calculate(fp B_min, fp B_max)
{
  m_eVals.clear();
  QVector< fp > field = findRoots(resonantSegments(B_min, B_max));

  if (field.isEmpty()) {
    qWarning() << "ATTENTION: no resonant segments found in range [" << B_min << ", " << B_max << "] for mwFreq = " << m_exp.mwFreqGHz;
    return field;
  }

  // cleanup field
  qSort(field);
  fp* it = field.begin();
  while(it != field.end() - 1) {
    if (abs(*it / *(it+1) - 1) < 1E-05) {
      it = field.erase(it);
    } else {
      ++it;
    }
  }

  m_eVals.clear();
  return field;
}

BisectAnswer ResonanceField::checkSegment(const BisectNode& min, const BisectNode& max) const
{
  const fp B_diff = max.B - min.B;
  bool resonancePossible = false;
  if ((max.E(m_exp.dimension - 1) - max.E(0)) > m_mwFreq) {
    if (!m_loopingResonanceCanOccur) {
      // for all kombinations u,v do eq 13
      for(int u = 0; u < m_exp.dimension; ++u) {
        for(int v = u + 1; v < m_exp.dimension; ++v) {
          // R_{uv}(B_q) * R_{uv}(B_r) <= 0
          if (((min.E(v) - min.E(u) - m_mwFreq) * (max.E(v) - max.E(u) - m_mwFreq)) <= 0) {
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
          if (abs((min.E(v) - min.E(u) + max.E(v) - max.E(u)) * 0.5 - m_mwFreq) <= m_lambda * B_diff) {
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
    const BisectNode mid = diagonalizeNode((min.B + max.B) * 0.5);
    fp epsilon = 0;
    for(int u = 0; u < m_exp.dimension; ++u) {
      const fp E_u_tilde = 0.5 * (min.E(u) + max.E(u)) + B_diff / 8.0 * (min.E_deriv(u) - max.E_deriv(u));
      fp epsilon_u = abs(mid.E(u) - E_u_tilde);
      if (epsilon_u > epsilon) {
        epsilon = epsilon_u;
      }
    }
    epsilon *= 2;

    if (epsilon > 1.0E-5 * m_mwFreq) {
      return BisectAnswer::continueAnswer(min.B, max.B, mid);
    } else {
      return BisectAnswer::resonantAnswer(min.B, max.B, mid);
    }
  } else {
    return BisectAnswer::notResonantAnswer(min.B, max.B);
  }
}

QMap< fp, fp > ResonanceField::resonantSegments(fp B_minStart, fp B_maxStart)
{
  //STEP 1: find knots
  QMap<fp, fp> resonantSegments;

  m_eVals[B_minStart] = diagonalizeNode(B_minStart);
  m_eVals[B_maxStart] = diagonalizeNode(B_maxStart);

  QStack<Segment> segments;
  segments << Segment(B_minStart, B_maxStart);

  ///FIXME: how to parallelize this?
  while (!segments.isEmpty()) {
    const Segment s = segments.pop();
    const BisectAnswer answer = checkSegment(m_eVals.value(s.B_min), m_eVals.value(s.B_max));
    switch(answer.status) {
      case BisectAnswer::NotResonant:
        // nothing to do
        break;
      case BisectAnswer::Continue:
        m_eVals[answer.mid.B] = answer.mid;
        segments << Segment(answer.from, answer.mid.B);
        segments << Segment(answer.mid.B, answer.to);
        break;
      case BisectAnswer::Resonant:
        m_eVals[answer.mid.B] = answer.mid;
        resonantSegments[answer.from] = answer.mid.B;
        resonantSegments[answer.mid.B] = answer.to;
        break;
    }
  }

  if (resonantSegments.isEmpty()) {
    qWarning() << "ATTENTION: no resonant segments found in range [" << B_minStart << ", " << B_maxStart << "] for mwFreq = " << (m_exp.mwFreqGHz);
  }

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

QVector<fp> ResonanceField::findRoots(const QMap<fp, fp>& resonantSegments)
{
  //STEP 2: find roots
  QVector<fp> resonanceField;

  // see eq 8
  const Matrix4 M = (Matrix4() << 
                      2, -2, 1, 1,
                      -3, 3, -2, -1,
                      0, 0, 1, 0,
                      1, 0, 0, 0).finished();

  QMap<fp, fp>::const_iterator it = resonantSegments.constBegin();
  const QMap<fp, fp>::const_iterator end = resonantSegments.constEnd();
  GNUPLOT_DEBUG(
    cout << "set xrange[" << it.key() * 0.9 << ":" << (end-1).value() * 1.1 << "]" << endl;
    cout << "u(x,min,max) = (x>=min)&&(x<max)? 1 : 1/0;" << endl;
    QString _nodes;
    QTextStream nodes(&_nodes);
  )
  while(it != end) {
    const fp B_min = it.key();
    const fp B_max = it.value();
    ++it;
    const BisectNode& min = m_eVals.value(B_min);
    const BisectNode& max = m_eVals.value(B_max);

    const fp B_diff = B_max - B_min;
    GNUPLOT_DEBUG(
      cout << "set arrow from " << B_min << ",-5e-26 to " << B_min << ",5e-26 ls 2" << endl;
      cout << "set arrow from " << B_max << ",-5e-26 to " << B_max << ",5e-26 ls 2" << endl;
    )
    for(int u = 0; u < m_exp.dimension; ++u) {
      const Vector4 e_u = (Vector4() << min.E(u), max.E(u), B_diff * min.E_deriv(u), B_diff * max.E_deriv(u)).finished();
      for(int v = u + 1; v < m_exp.dimension; ++v) {
        const Vector4 e_v = (Vector4() << min.E(v), max.E(v), B_diff * min.E_deriv(v), B_diff * max.E_deriv(v)).finished();
        ///NOTE: paper has different notation: p(0) == p_3, p(1) == p_2, ...
        const Vector4 p = M * (e_v - e_u);
        fp root = 0;
        if (!m_loopingResonanceCanOccur) {
          if (((min.E(v) - min.E(u) - m_mwFreq) * (max.E(v) - max.E(u) - m_mwFreq)) > 0) {
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
          root = B_min + t2 * B_diff;
        } else {
          qDebug() << "NOT IMPLEMENTED YET";
        }

        resonanceField << root;
        GNUPLOT_DEBUG(
          cout << "replot u(x, " << B_min << "," << B_max << ") * ((((x-" << B_min << ")/" << B_diff << "))**3 * (" << p(0) << ") + (((x-" << B_min << ")/" << B_diff << "))**2 * (" << p(1) << ") + (((x-" << B_min << ")/" << B_diff << ")) * (" << p(2) << ") + (" << p(3) << " - " << mwFreq << ")) "
                   "title \"B_min = " << B_min << ", B_max = " << B_max << " || u = " << u << ", v = " << v << "\"" << endl;
          cout << "set arrow from " << root << ",-1e-24 to " << root << ",1e-24 ls 0" << endl;
          nodes << B_min << '\t' << (min.E(v) - min.E(u) - mwFreq) << endl;
          nodes << B_max << '\t' << (max.E(v) - max.E(u) - mwFreq) << endl;
        )
      }
    }
  }
  GNUPLOT_DEBUG(
    cout << flush;
    qDebug() << _nodes;
  )

  return resonanceField;
}
