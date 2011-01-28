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

#ifndef MW_BACHELOR_EXPERIMENT_H
#define MW_BACHELOR_EXPERIMENT_H

/**
 * Experiment data and physical constants
 */
namespace Experiment {
#ifndef PROTONS
const int nprotons = 7;
#else
const int nprotons = PROTONS;
#endif

const int dimension = pow(2, nprotons + 1);

//gtensor=============================
static inline Matrix3cd gTensor()
{
  return Matrix3cd::Identity() * 2.0022838;
}

//static field========================
static inline Vector3cd staticBField(const double B)
{
  const double B0_norm = B; //0.2839; //field in Tesla
  double B0[3] = {0, 0, 1}; //field direction
  const double B_temp = sqrt( B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2] );

  for (int i = 0; i<3; i++) {
    B0[i] = B0[i]*B0_norm/B_temp;
  }

  //cout << endl << "B0:" << " " << B0[0] << '\t' << B0[1] << '\t' << B0[2] << endl;
  //cout.precision(5);
  //cout << "\nStatic Field (T): " << sqrt(B0[0]*B0[0]+B0[1]*B0[1]+B0[2]*B0[2]) << endl;

  return Vector3cd(B0[0], B0[1], B0[2]);
}

//arbitrary A-tensor==================
static inline Matrix3cd aTensor()
{
  return Matrix3cd::Identity() * 1420;  //proton hyperfine coupling in T
}

};

#endif // MW_BACHELOR_EXPERIMENT_H
