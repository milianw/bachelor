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

#ifndef MW_BACHELOR_TYPES_H
#define MW_BACHELOR_TYPES_H

#include <complex>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

typedef double fp;
// typedef float fp;

typedef complex<fp> c_fp;

typedef Matrix<fp, 2, 2> Matrix2;
typedef Matrix<c_fp, 2, 2> Matrix2c;

typedef Matrix<fp, 3, 3> Matrix3;
typedef Matrix<c_fp, 3, 3> Matrix3c;

typedef Matrix<fp, 4, 4> Matrix4;
typedef Matrix<c_fp, 4, 4> Matrix4c;

typedef Matrix<fp, Dynamic, Dynamic> MatrixX;
typedef Matrix<c_fp, Dynamic, Dynamic> MatrixXc;

typedef Matrix<fp, 3, 1> Vector3;
typedef Matrix<c_fp, 3, 1> Vector3c;

typedef Matrix<fp, 4, 1> Vector4;
typedef Matrix<c_fp, 4, 1> Vector4c;

typedef Matrix<fp, Dynamic, 1> VectorX;
typedef Matrix<c_fp, Dynamic, 1> VectorXc;

#endif // MW_BACHELOR_TYPES_H
