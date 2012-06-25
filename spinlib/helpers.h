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

#ifndef MW_BACHELOR_HELPERS_H
#define MW_BACHELOR_HELPERS_H

#include <ostream>

#include "fptype.h"
#include "eigentypes.h"

class Experiment;

Experiment getExperiment(const string& orcaInput, int cutoffcount, int protons, int nitrogens);

std::string formatSize(long long unsigned int size);

void printExperiment(std::ostream& out, const Experiment& exp);

std::string identifierForExperiment(const Experiment& exp);

std::string guessPeakMemConsumption(const Experiment& exp);

void guessBRange(const Experiment& exp, fp& from, fp& to);

/**
 * generates rotation matrix to transform from basis @p from to basis @p to.
 *
 * to rotate a vector into the basis: v' = R * v;
 * to rotate a tensor into the basis: a' = R * a * R.transpose()
 */
Matrix3c rotationMatrix(const Matrix3c& to, const Matrix3c& from);

#endif // MW_BACHELOR_HELPERS_H
