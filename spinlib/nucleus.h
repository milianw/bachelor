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

#ifndef MW_BACHELOR_NUCLEUS_H
#define MW_BACHELOR_NUCLEUS_H

#include "eigentypes.h"

/**
 * defines a nucleus in the experiment
 */
struct Nucleus {
  Nucleus() {}
  Nucleus(const std::string& name_, const int twoJ_, const int isotope_,
          const Matrix3& A_, const fp g_, const Vector3 Q_, const Matrix3 EFG_)
  : name(name_), twoJ(twoJ_), isotope(isotope_), A(A_), g(g_)
  , Q(Q_), EFG(EFG_)
  { }

  // arbitrary name identifying the nucleus
  std::string name;
  // 2*J (positive spin maximum)
  int twoJ;
  // isotope number
  int isotope;
  // hyperfine coupling matrix in MHz
  Matrix3 A;
  // g value
  fp g;
  // quadrupole tensor eigenvalues in MHz
  // i.e. the diagonal representation of the SH term I*Q*I = e**2qQ/(4I(2I-1))*[-(1-eta),-(1+eta),2]
  Vector3 Q;
  // electric field gradient matrix
  Matrix3 EFG;

  /// return full quadrupole coupling matrix
  /// if EFG is not zero this will re-construct the original
  /// Q-matrix from the principal values in Q
  /// otherwise this will just return a diagonal matrix
  /// with the principal values
  Matrix3c quadrupoleCouplingMatrix() const
  {
    Matrix3c ret = Matrix3c::Zero();
    if (EFG.isZero()) {
      // use principal values only
      ret.diagonal() = Q.cast<c_fp>();
      return ret;
    }

    // full Q tensor from EFG matrix
    // A = e^2qQ / (2I*(2I-1)) = Q(2)
    const fp A = Q(2);
    Vector3c V = EFG.eigenvalues();
    c_fp V_zz = V(0);
    for(int i = 1; i < 3; ++i) {
      if (abs(V(i)) > abs(V_zz)) {
        V_zz = V(i);
      }
    }
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        // Q_ij = A * V_ij/V_zz
        ret(i, j) = A * EFG(i, j) / V_zz;
      }
    }

    return ret;
  }
};


#endif // MW_BACHELOR_NUCLEUS_H
