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

#include <cmath>

#include "constants.h"
#include "eigentypes.h"

/**
 * Experiment data / simulation parameters
 */
struct Experiment {
  Experiment(int _nProtons)
  : nProtons(_nProtons)
  , dimension(pow(2, nProtons + 1))
  , mwFreqGHz(0)
  {
    gTensor = Matrix3c::Identity() * Constants::g_E;
    aTensor = Matrix3c::Identity() * 1420;
    staticBFieldDirection << 0, 0, 1;
  }

  const int nProtons;
  const int dimension;

  Matrix3c gTensor;
  Matrix3c aTensor;
  Vector3c staticBFieldDirection;
  fp mwFreqGHz;

  /**
   * static field of given strength in direction of @c staticBFieldDirection
   *
   * @p B field strength in tesla
   */
  inline Vector3c staticBField(const fp B) const
  {
    return staticBFieldDirection * B / staticBFieldDirection.norm();
  }
};

#endif // MW_BACHELOR_EXPERIMENT_H
