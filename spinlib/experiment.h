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

#include "eigentypes.h"
#include "nucleus.h"

#include <vector>

class Spins;

/**
 * Experiment data / simulation parameters
 *
 * so far assumes a single unpaired electron
 */
struct Experiment {
  /**
   * construct experiment out of list of nuclei
   *
   * electron gTensor defaults to Constants::g_E * Matrix3c::Identity();
   * static B Field direction defaults to [0 0 1], i.E. only in z-Direction
   * mwFreqGHz defaults to zero
   */
  explicit Experiment(const std::vector<Nucleus>& nuclei);

  /**
   * generate experiment with dummy data
   *
   * all nuclei will have an A-Tensor of 1420 * Matrix3c::Identity()
   * for the other defaults see the Experiment() ctor
   *
   * @p protons number of 1H nuclei (J = 0.5, g = Constants::g_1H)
   * @p nitrogens number of 14N nuclei (J = 0.5, g = Constants::g_14N)
   */
  static Experiment generateDummy(int protons, int nitrogens);

  /**
   * generate spin system for this experiment
   */
  Spins spinSystem() const;

  // nuclei in the experiment
  // guaranteed sort order by Nucleus.twoJ (ascending)
  const std::vector<Nucleus> nuclei;
  // dimensionality of the experiment, basically this is:
  // \Pi_k (2J_k + 1)
  // in our case: 2^protons * 3^nitrogens
  const int dimension;

  // electron g Tensor
  Matrix3c gTensor;
  // direction of static B field
  Vector3c staticBFieldDirection;
  // incident micro wave frequency in GHz
  fp mwFreqGHz;

  /**
   * static field of given strength in direction of |g_z>
   *
   * @p B field strength in tesla
   */
  Vector3c staticBField(const fp B) const;
};

#endif // MW_BACHELOR_EXPERIMENT_H
