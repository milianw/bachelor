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

#include "types.h"

class QTextStream;
class Experiment;

/**
 * TODO: check whether inlining some parts is noticeable
 *
 * notes on porting:
 *
 * ~~~~~~~~~~~
 * gsl_blas_zdotu (x, y, dotu)
 * => eigen: dotu = x.conjugate().dot(y)
 *
 * from eigen docs about .dot():
 * Note: If the scalar type is complex numbers, then this function returns the hermitian (sesquilinear) dot product,
 *       conjugate-linear in the first variable and linear in the second variable.
 *       but gsl_blas_zdotu seems to differ from this, hence use a.conjugate().dot(b) instead of a.dot(b)
 * NOTE: <orzel> in this very specific case, it might be that (a.transpose()*b)(0,0) is faster (taking the only element of the 1x1 matrix x^T.y
 * ~~~~~~~~~~~
 * gsl_blas_zgemv(CblasTrans, gsl_complex_rect(1,0), atensor, I, gsl_complex_rect(0,0), atensor_I);
 * => eigen: atensor_I = atensor.transpose() * I
 * ~~~~~~~~~~~
 * gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1,0),
 *                eigenvectors, moments, gsl_complex_rect(1,0), intermediate);
 * docs: C = \alpha op(A) op(B) + \beta C
 * => eigen: intermediate = eigenvectors.adjoint() * moments + intermediate
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
    void calculateIntensity(QTextStream* out) const;

    VectorX calculateEigenValues() const;

    /// the complete hamiltonian
    MatrixXc hamiltonian() const;
    /// nuclear Zeeman Hamiltonian component
    MatrixXc nuclearZeeman() const;
    /// hyper fine Hamiltonian component
    MatrixXc hyperFine() const;
    /// electron Zeeman Hamiltonian component
    MatrixXc electronZeeman() const;

  private:
    /// interprets @p i as binary number and returns the k-th bit of it
    inline bool spinState(int i, int k) const;

    /// return spin vector from pauli matrices
    inline Vector3c spinVector(int i, int j, int k) const;

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
    inline bool stateContributes(int i, int j, int k, bool ignoreElectron = IgnoreElectron) const;

    /// moments
    /// TODO: better document
    MatrixXc magneticMoments() const;

    inline c_fp magneticMoment(const int i, const int j) const;

    /// intensitiy matrix with coefficients (i, j) = |< psi_j | M | psi_i>|^2
    /// psi_i being the i-th eigen vector
    /// M being the magnetic moment matrix
    MatrixX intensityMatrix(const MatrixXc& eigenVectors) const;

    const fp m_B;
    const Experiment& m_exp;
    const Vector3c m_staticBField;
};

#endif // MW_BACHELOR_SPINHAMILTONIAN_H
