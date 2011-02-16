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


#ifndef MW_BACHELOR_RESONANCEFIELD_H
#define MW_BACHELOR_RESONANCEFIELD_H

#include <QtCore/QVector>
#include <QtCore/QStack>
#include <QtCore/QMap>
#include <QtCore/QDebug>

#include "spinhamiltonian.h"

/**
 * Implementation of "S. Stoll, A. Schweiger / Chemical Physics Letters 380 (2003) 464 - 470"
 *
 * NOTE: the paper uses the hamiltonian in terms of h, while we use atomic units. meaning:
 * paper * h => SI value
 * us / h => SI value
 * hence most algorithms in the paper above have to be multiplied by h *twice*
 *
 * NOTE: mwFreq is in GHz!
 *
 * NOTE: This is WIP
 */
class ResonanceField {
public:
  ResonanceField(const Experiment& exp);
  QVector<fp> findRoots(fp B_min, fp B_max, fp mwFreq);

private:
  bool needSignChangeCheck(fp mwFreq) const;
  fp calculateLambda() const;

  struct Knot {
    Knot()
    : B_min(0), B_max(0)
    {
    }

    Knot(fp min, fp max)
    : B_min(min), B_max(max)
    {
    }

    bool operator==(const Knot& other) const
    {
      return B_min == other.B_min && B_max == other.B_max;
    }

    fp B_min;
    fp B_max;
  };
  // stores E_u and E'_u for all u of a given B
  struct EigenValues {
    EigenValues()
    {
    }
    EigenValues(fp B, const Experiment& exp)
    : E_deriv(exp.dimension)
    {
      SpinHamiltonian H(B, exp);
      SelfAdjointEigenSolver<MatrixXc> eigenSolver(H.hamiltonian());
      E = eigenSolver.eigenvalues();
      const MatrixXc eigenVectors = eigenSolver.eigenvectors();
      const MatrixXc G = H.nuclearZeeman() + H.electronZeeman();
      for(int u = 0; u < exp.dimension; ++u) {
        // <u| G |u> => expectation value is always real
        E_deriv(u) = (eigenVectors.col(u).adjoint() * G * eigenVectors.col(u))(0, 0).real();
      }
    }
    VectorX E;
    VectorX E_deriv;
  };
  const Experiment& m_exp;
};

ResonanceField::ResonanceField(const Experiment& exp)
: m_exp(exp)
{

}

bool ResonanceField::needSignChangeCheck(fp mwFreq) const
{
  ///TODO: make sure they are always sorted in the correct order...
  VectorX eVals = SpinHamiltonian(0, m_exp).calculateEigenValues();
  qDebug() << "sign change:" << (eVals(m_exp.dimension - 1) - eVals(0)) << mwFreq;
  return (eVals(m_exp.dimension - 1) - eVals(0)) < mwFreq;
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
  qDebug() << "lambda:" << lambda;
  return lambda;
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

QVector<fp> ResonanceField::findRoots(fp in_B_min, fp in_B_max, fp _mwFreq)
{
  // to atomic units:
  fp mwFreq = _mwFreq * 1E9 * h;
  qDebug() << " " << mwFreq;

  //STEP 1: find knots
  QMap<fp, EigenValues> eVals;
  QMap<fp, fp> resonantSegments;

  ///TODO: rename method
  const bool loopingResonanceCanOccur = !needSignChangeCheck(mwFreq);
  {
    const fp lambda = calculateLambda();

    eVals[in_B_min] = EigenValues(in_B_min, m_exp);
    eVals[in_B_max] = EigenValues(in_B_max, m_exp);

    QStack<Knot> knots;
    knots << Knot(in_B_min, in_B_max);

    ///FIXME: how to parallelize this?
    while (!knots.isEmpty()) {
      Knot knot = knots.pop();
      const EigenValues& max = eVals.value(knot.B_max);
      const EigenValues& min = eVals.value(knot.B_min);
      const fp B_diff = knot.B_max - knot.B_min;
      bool resonancePossible = false;
//       qDebug() << max.E(m_exp.dimension - 1) << " - " << max.E(0) << " = " << (max.E(m_exp.dimension - 1) - max.E(0)) << " > " << mwFreq;
      if ((max.E(m_exp.dimension - 1) - max.E(0)) > mwFreq) {
        if (!loopingResonanceCanOccur) {
//           qDebug() << "it's bigger, yay and no looping";
          // for all kombinations u,v do eq 13
          for(int u = 0; u < m_exp.dimension; ++u) {
            for(int v = u + 1; v < m_exp.dimension; ++v) {
              // R_{uv}(B_q) * R_{uv}(B_r) <= 0
//               qDebug() << "(" << min.E(v) << " - " << min.E(u) << ") * (" << max.E(v) << " - " << max.E(u) << ") = " << (min.E(v) - min.E(u)) << " * " << (max.E(v) - max.E(u)) << " = " << ((min.E(v) - min.E(u)) * (max.E(v) - max.E(u)));
              if (((min.E(v) - min.E(u) - mwFreq) * (max.E(v) - max.E(u) - mwFreq)) <= 0) {
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
              if (abs((min.E(v) - min.E(u) + max.E(v) - max.E(u)) * 0.5 - mwFreq) <= lambda * B_diff) {
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
//       qDebug() << "knot.B_min = " << knot.B_min << ", knot.B_max = " << knot.B_max << " ~ resonanct segment?" << resonancePossible;
      if (resonancePossible) {
        const fp B_new = (knot.B_min + knot.B_max) * 0.5;

        // eq 11
        EigenValues _new(B_new, m_exp);
        fp epsilon = 0;
        for(int u = 0; u < m_exp.dimension; ++u) {
          const fp E_u_tilde = 0.5 * (min.E(u) + max.E(u)) + B_diff / 8.0 * (min.E_deriv(u) - max.E_deriv(u));
          fp epsilon_u = abs(_new.E(u) - E_u_tilde);
          if (epsilon_u > epsilon) {
            epsilon = epsilon_u;
          }
        }
        epsilon *= 2;

        eVals[B_new] = _new;
        if (epsilon > 1.0E-3 * mwFreq) {
//           qDebug() << "adding new knot:" << B_new << "epsilon:" << epsilon <<  "VS:" << 1.0E-3 * mwFreq;
          knots << Knot(knot.B_min, B_new) << Knot(B_new, knot.B_max);
        } else {
//           qDebug() << "resonant segment approximated:" << knot.B_min << knot.B_max;
//           resonantSegments[knot.B_min] = knot.B_max;
          resonantSegments[knot.B_min] = B_new;
          resonantSegments[B_new] = knot.B_max;
        }
      }
    }
  }
  qDebug() << "resonant segments:" << resonantSegments;

  if (resonantSegments.isEmpty()) {
    qWarning() << "ATTENTION: no roots found in range [" << in_B_min << ", " << in_B_max << "] for mwFreq = " << mwFreq;
    return QVector<fp>();
  }
  //STEP 2: find roots
  QVector<fp> resonanceField;

  // see eq 8
  const Matrix4 M = (Matrix4() << 
                      2, -2, 1, 1,
                      -3, 3, -2, -1,
                      0, 0, 1, 0,
                      1, 0, 0, 0).finished();

//   QMap<fp, EigenValues>::const_iterator it = eVals.constBegin();
//   const QMap<fp, EigenValues>::const_iterator end = eVals.constEnd();
  QMap<fp, fp>::const_iterator it = resonantSegments.constBegin();
  const QMap<fp, fp>::const_iterator end = resonantSegments.constEnd();
  cout << "set xrange[" << it.key() * 0.9 << ":" << (end-1).value() * 1.1 << "]" << endl;
  cout << "u(x,min,max) = (x>=min)&&(x<max)? 1 : 1/0;" << endl;
  while(it != end) {
    const fp B_min = it.key();
    const fp B_max = it.value();
    ++it;
    const EigenValues& min = eVals.value(B_min);
    const EigenValues& max = eVals.value(B_max);
    /*
    const EigenValues& min = it.value();
    const fp B_min = it.key();
    ++it;
    const EigenValues& max = it.value();
    const fp B_max = it.key();
    */

    const fp B_diff = B_max - B_min;
//     qDebug() << "find roots between:" << B_min << B_max << B_diff;
    cout << "set arrow from " << B_min << ",-5e-26 to " << B_min << ",5e-26 ls 2" << endl;
    cout << "set arrow from " << B_max << ",-5e-26 to " << B_max << ",5e-26 ls 2" << endl;

    for(int u = 0; u < m_exp.dimension; ++u) {
      const Vector4 e_u = (Vector4() << min.E(u), max.E(u), B_diff * min.E_deriv(u), B_diff * max.E_deriv(u)).finished();
      for(int v = u + 1; v < m_exp.dimension; ++v) {
        const Vector4 e_v = (Vector4() << min.E(v), max.E(v), B_diff * min.E_deriv(v), B_diff * max.E_deriv(v)).finished();
        ///NOTE: paper has different notation: p(0) == p_3, p(1) == p_2, ...
        const Vector4 p = M * (e_v - e_u);
        fp root = 0;
        if (!loopingResonanceCanOccur) {
          if (((min.E(v) - min.E(u) - mwFreq) * (max.E(v) - max.E(u) - mwFreq)) > 0) {
            continue;
          }
          // Newton-Raphson root finding
          /* TODO: how is this supposed to work?
          fp E_diff_min = (min.E(v) - min.E(u));
          if (E_diff_min < 1E-24) {
            continue;
          }
          fp E_diff_max = (max.E(v) - max.E(u));
//           qDebug() << "Delta_uv(B_q) =" << E_diff_min << ", Delta_uv(B_r)" << E_diff_max;
          fp t = - E_diff_min / (E_diff_max - E_diff_min);
          */
//           fp t = newtonRaphson(p, 0, mwFreq); // x_0 = 0
          fp t = newtonRaphson(p, 0.5, mwFreq); // x_0 = 0.5
          if (!isfinite(t)) {
            continue;
          }
//           qDebug() << "initial root:" << t << B_min + t * B_diff;
          fp t2 = t;
          const int maxIt = 10000;
          int it = 0;
          do {
            t = t2;
            t2 = newtonRaphson(p, t, mwFreq);
//             qDebug() << t << t2;
          ///TODO: when to abort?
          } while(abs(t2/t - 1) > 1E-5 && (maxIt > ++it));
          root = B_min + t2 * B_diff;
        } else {
          qDebug() << "NOT IMPLEMENTED YET";
        }

//           qDebug() << "found root:" << B_min << " + " << t2 << " * " << B_diff << " = " << B_min << " + " << (t2 * B_diff) << " = " << root;
        resonanceField << root;
//         cout << "e = " << (e_v - e_u).transpose() << endl << "p = " << p.transpose() << endl;
        cout << "replot u(x, " << B_min << "," << B_max << ") * ((((x-" << B_min << ")/" << B_diff << "))**3 * (" << p(0) << ") + (((x-" << B_min << ")/" << B_diff << "))**2 * (" << p(1) << ") + (((x-" << B_min << ")/" << B_diff << ")) * (" << p(2) << ") + (" << p(3) << " - " << mwFreq << ")) "
                     "title \"B_min = " << B_min << ", B_max = " << B_max << " || u = " << u << ", v = " << v << "\"" << endl;
        cout << "set arrow from " << root << ",-1e-24 to " << root << ",1e-24 ls 0" << endl;
      }
    }
  }
  cout << flush;
  qSort(resonanceField);
  qDebug() << "resonance field:" << resonanceField;

  return resonanceField;
}


#endif // MW_BACHELOR_RESONANCEFIELD_H
