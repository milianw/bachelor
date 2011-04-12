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

#include <vector>
#include <map>

#include "eigentypes.h"

class Experiment;

class BisectNode;
class BisectAnswer;

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
  std::vector<fp> calculate(fp B_min, fp B_max);

  /**
   * Diagonalize spin hamiltonian for @p B and save eigen values to @p E
   * and derivatives to @p E_deriv
   */
  BisectNode diagonalizeNode(const fp B) const;

  /**
   * Check B-segment @p from @p to and decide whether it is resonant or not.
   */
  BisectAnswer checkSegment(const BisectNode& from, const BisectNode& to) const;

  /**
   * Find roots in given segment @p from @p to.
   */
  std::vector<fp> findRootsInSegment(const BisectNode& from, const BisectNode& to) const;

  /**
   * Cleanup resonancy field, removes nearby B ranges below resolution threshold.
   */
  void cleanupResonancyField(std::vector<fp>& field) const;

private:
  const Experiment& m_exp;
  // micro wave frequency in atomic units
  const fp m_mwFreq;
  const fp m_lambda;

  const fp m_cleanupThreshold;
  const fp m_resfieldThreshold;

  fp calculateLambda() const;

  /// @return true if looping resonance can occur, false otherwise
  bool checkForLoopingResonance() const;
  mutable bool m_loopingResonanceChecked;
  mutable bool m_loopingResonanceCanOccur;

  std::map<fp, fp> resonantSegments(fp B_minStart, fp B_maxStart);
  std::vector<fp> findRoots(const std::map<fp, fp>& resonantSegments);

  std::map<fp, BisectNode> m_eVals;
};

class BisectNode {
public:
  BisectNode()
  { }

  BisectNode(fp _B, const VectorX& _E, const VectorX& _E_deriv)
  : B(_B), E(_E), E_deriv(_E_deriv)
  { }

  fp B;
  VectorX E;
  VectorX E_deriv;
};

class BisectAnswer {
public:
  BisectAnswer()
  { }

  enum Status {
    Continue, Resonant, NotResonant
  };

  Status status;
  fp from;
  fp to;

  BisectNode mid;
};
#endif // MW_BACHELOR_RESONANCEFIELD_H
