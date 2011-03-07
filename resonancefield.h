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
#include <QtCore/QMap>

#include "types.h"

class Experiment;

/**
 * Implementation of "S. Stoll, A. Schweiger / Chemical Physics Letters 380 (2003) 464 - 470"
 *
 * NOTE: the paper uses the hamiltonian in terms of h, while we use atomic units. meaning:
 * paper * h => SI value
 * us / h => SI value
 * hence most algorithms in the paper above have to be multiplied by h *twice*
 */
class ResonanceField {
public:
  explicit ResonanceField(const Experiment& exp);
  /**
   * Calculate resonance Field for given B-range and micro wave frequency.
   *
   * @p B_min minimum static B-Field in Tesla
   * @p B_max maximum static B-Field in Tesla
   */
  QVector<fp> calculate(fp B_min, fp B_max);

  /**
   * Diagonalize spin hamiltonian for @p B and save eigen values to @p E
   * and derivatives to @p E_deriv
   */
  BisectNode diagonalizeNode(const fp B) const;

  BisectAnswer checkSegment(const BisectNode& from, const BisectNode& to) const;

private:
  const Experiment& m_exp;
  // micro wave frequency in atomic units
  const fp m_mwFreq;
  const bool m_loopingResonanceCanOccur;
  const fp m_lambda;

  fp calculateLambda() const;

  /// @return true if looping resonance can occur, false otherwise
  bool checkForLoopingResonance() const;

  QMap<fp, fp> resonantSegments(fp B_minStart, fp B_maxStart);
  QVector<fp> findRoots(const QMap<fp, fp>& resonantSegments);

  QMap<fp, BisectNode> m_eVals;
};


#endif // MW_BACHELOR_RESONANCEFIELD_H
