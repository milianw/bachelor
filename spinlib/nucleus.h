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
 *
 * NOTE: use MKS units!
 */
struct Nucleus {
  Nucleus() {}
  Nucleus(const std::string& name_, const int twoJ_,
          const Matrix3c& A_, const fp g_)
  : name(name_), twoJ(twoJ_), A(A_), g(g_)
  { }

  std::string name;
  int twoJ;
  ///TODO: complex data
  Matrix3c A;
  fp g;
};


#endif // MW_BACHELOR_NUCLEUS_H
