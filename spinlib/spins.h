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

#ifndef MW_BACHELOR_SPINS_H
#define MW_BACHELOR_SPINS_H

#include <iostream>
#include "eigentypes.h"

// pauli matrices for J = 1/2
namespace PauliMatrix_J_half {
  const Matrix2c X = (Matrix2c() << 0, 0.5,
                                    0.5, 0).finished();
  const Matrix2c Y = (Matrix2c() << 0, c_fp(0, -0.5),
                                    c_fp(0, 0.5), 0).finished();
  const Matrix2c Z = (Matrix2c() << 0.5, 0,
                                    0, -0.5).finished();
}
// pauli matrices for J = 1
namespace PauliMatrix_J_one {
  const Matrix3c X = (Matrix3c() << 0, 0.5 * sqrt(2), 0,
                                    0.5 * sqrt(2), 0, 0.5 * sqrt(2),
                                    0, 0.5 * sqrt(2), 0).finished();
  const Matrix3c Y = (Matrix3c() << 0, c_fp(0, -0.5 * sqrt(2)), 0,
                                    c_fp(0, 0.5 * sqrt(2)), 0, c_fp(0, -0.5 * sqrt(2)),
                                    0, c_fp(0, 0.5 * sqrt(2)), 0).finished();
  const Matrix3c Z = (Matrix3c() << 1, 0, 0,
                                    0, 0, 0,
                                    0, 0, -1).finished();
}

/**
 * Class that holds info about spins and possible states.
 *
 * \code
 * // assume 2 x J = 1/2, 1 x J = 1
 * // e.g. one electron, one nucleus with J = 1/2 and one with  J = 1
 * Spins system(2, 1);
 * // iterate over all states
 * for (int bra = 0; bra < system.states; ++bra) {
 *   for (int ket = 0; ket < system.states; ++ket) {
 *     // electron spin vector:
 *     Vector3c S = system.spinVector(bra, ket, 0);
 *     // J = 1/2 nucleus vector
 *     Vector3c I_half = system.spinVector(bra, ket, 1);
 *     // J = 1 nucleus vector
 *     Vector3c I_one = system.spinVector(bra, ket, 2);
 *   }
 * }
 * \endcode
 */
struct Spins
{
  /**
   * @p _spinHalfs total number of J = 1/2 in system
   * @p _spinOnes total number of J = 1 in system
   */
  explicit Spins(int _spinHalfs, int _spinOnes)
  : spinHalfs(_spinHalfs)
  , spinOnes(_spinOnes)
  , elements(spinHalfs + spinOnes)
  , states(pow(2, spinHalfs) * pow(3, spinOnes))
  , dimPow_one(pow(3, spinOnes))
  , elemPows_one(cacheElemPows(3, spinOnes))
  {

  }

  /// @return spin of @p spinId in @p state
  int spinInState(const int spinId, const int state) const
  {
    // decide which spin this spinId has
    // actually we calc 2J
    const int twoJ = spinId < spinHalfs ? 1 : 2;
    if (twoJ == 1) {
      /// optimization for twoJ = 1, i.e. spin 1/2
      return bool(state & (1 << spinId));
    }
    // number of spins with same dimension
    const int relatedSpins = spinOnes;
    // position of spinId in the related spins
    const int posInRelatedSpins = spinId - spinHalfs;
    if (posInRelatedSpins >= relatedSpins) {
      std::cerr << "invalid spin pos:" << posInRelatedSpins << "for spin id:" << spinId << ", elements:" << elements << ", spin halfs:" << spinHalfs << ", spin ones:" << spinOnes;
      throw 1;
    }

    // dimension is 2J + 1
    const int dim = twoJ + 1;
    // number of states created by spins of this type
    // generic version: const int dimPow = pow(dim, relatedSpins);
    // cached version for J = 1
    const int dimPow = dimPow_one;
    // modulate state by dimPow to get a number representing
    // the sate of spins of the current type
    const int dimState = state % dimPow;
    // now interpret dimstate as number in base dim:
    // dimState = dim^0 * a_0 + dim^1 * a_1 + ... dim^N * a_n
    // and return a_{posInRelatedSpin},
    // for dim = 2 this is equal to dimState & posInRelatedSpins
    // for arbitrary dims we first calculate dim^{posInRelatedSpins}
    // and devide dimState by it and modulate again by dim
    // generic version: const int elemPow = pow(dim, posInRelatedSpins)
    // cached version for J = 1
    const int elemPow = elemPows_one[posInRelatedSpins];
    return dimState / elemPow % dim;
  }

  /// @return spin vector for @p spinId in @p bra and @p ket configuration
  Vector3c spinVector(const int bra, const int ket, const int spinId) const
  {
    const int braSpin = spinInState(spinId, bra);
    const int ketSpin = spinInState(spinId, ket);
    if (spinId < spinHalfs) {
      // J = 1/2
      return (Vector3c() << PauliMatrix_J_half::X(braSpin, ketSpin), PauliMatrix_J_half::Y(braSpin, ketSpin), PauliMatrix_J_half::Z(braSpin, ketSpin)).finished();
    } else {
      // J = 1
      return (Vector3c() << PauliMatrix_J_one::X(braSpin, ketSpin), PauliMatrix_J_one::Y(braSpin, ketSpin), PauliMatrix_J_one::Z(braSpin, ketSpin)).finished();
    }
  }

  const int spinHalfs;
  const int spinOnes;
  const int elements;
  const int states;
  /// caches for pow calls
  const int dimPow_one;
  // size: relatedSpins, index: posInRelatedSpins
  const int* const elemPows_one;

private:
  const int* cacheElemPows(const int dim, const int size) const
  {
    int* ret = new int[size];
    int pow = 1;
    for(int i = 0; i < size; ++i) {
      ret[i] = pow;
      pow *= dim;
    }
    return ret;
  }
};

#endif // MW_BACHELOR_SPINS_H