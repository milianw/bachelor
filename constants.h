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

#ifndef MW_BACHELOR_CONSTANTS_H
#define MW_BACHELOR_CONSTANTS_H

#include <gsl/gsl_const_mksa.h>

/**
 * Physical constants, all in MKS units
 *
 * TODO: rename and make all uppercase
 */
namespace Constants {
const double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
const double h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H;
const double NUC_MAGNETON = GSL_CONST_MKSA_NUCLEAR_MAGNETON; // (J/T)
const double g_1H = 5.585694701 ;            // proton g factor
const double g_14N = 0.403761;               // N14 g factor
const double GAMMA_1H = 267.513E6 ;          // (rad/s.T)
const double GAMMA_14N =  19.331E6 ;         // (rad/s.T)
const double Bohrm = 9.27400949E-24;         // Bohr magneton in J/T
// Reminder: GAMMA_1H * hbar =  NUC_MAGNETON * g_1H
//const double B = 1;                          // Static magnetic field (T)
//const double B2 = ;                        // RF field
//const double LARMOR_1H = B * GAMMA_1H /2/M_PI / 1.0E6;
}

#endif // MW_BACHELOR_CONSTANTS_H
