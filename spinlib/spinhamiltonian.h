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

#include "eigentypes.h"
#include "spins.h"

class Experiment;

/**
 * TODO: check whether inlining some parts is noticeables
 *
 * NOTE: size of MatrixXc is (sizeof(complex< fp >) bytes * (2^(nProtons + 1)*3^(nNitrogens))^2)
 *                            ^ == 8 for float, 16 for double
 *
 * NOTE: if you have a dot product between two (spin-)operators, it is *not* the same as
 *       the dot product between two (complex) vectors!
 *       complex vectors: a . b = a.adjoint() * b
 *       spin operators: a . b = a.transpose() * b
 */
class SpinHamiltonian {
  public:
    /// @p B magnetic field in Tesla
    explicit SpinHamiltonian(const fp B, const Experiment& experiment);
    ~SpinHamiltonian();

    /// list all possible transitions with their frequency for the given B field
    void calculateTransitions() const;
    /// calculate total intensity of all transitions that are valid for the
    /// incoming microwave frequency (as set in @c Experiment)
    fp calculateIntensity() const;

    VectorX calculateEigenValues() const;

    /// return the complete spin hamiltonian
    MatrixXc hamiltonian() const;
    /// add the nuclear Zeeman Hamiltonian component to @p H
    void addNuclearZeeman(MatrixXc& H) const;
    /// add the hyper fine Hamiltonian component to @p H
    void addHyperFine(MatrixXc& H) const;
    /// electron Zeeman Hamiltonian component to @p H
    void addElectronZeeman(MatrixXc& H) const;

  private:
    /// interprets @p i as binary number and returns the k-th bit of it
    inline bool spinState(int state, int k) const;

    /// return spin vector from pauli matrices
    /// NOTE: don't .dot() them, use .cwiseProduct().sum()!
    inline Vector3c spinVector(int bra, int ket, int k) const;

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
    inline bool stateContributes(int bra, int ket, int k, bool ignoreElectron = IgnoreElectron) const;

    /// @return matrix of the |g_x> components of the magnetic moments in the system
    MatrixXc magneticMoments() const;

    /// intensitiy matrix with coefficients (i, j) = |< psi_j | M | psi_i>|^2
    /// psi_i being the i-th eigen vector
    /// M being the magnetic moment matrix
    MatrixX intensityMatrix(const MatrixXc& eigenVectors) const;

    const fp m_B;
    const Experiment& m_exp;
    const Spins m_spins;
    const Vector3c m_staticBField;
};

#endif // MW_BACHELOR_SPINHAMILTONIAN_H
