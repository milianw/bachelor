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

#ifndef MW_MPI_IFACE_H
#define MW_MPI_IFACE_H

const int MASTER_RANK = 0;

enum Commands {
  CMD_CLOSE,
  CMD_BISECT,
  CMD_DIAGONALIZE,
  CMD_FINDROOTS
//   CMD_INTENSITY
};

enum Tags {
  TAG_CMD,
  TAG_DIAGONALIZE_INPUT,
  TAG_DIAGONALIZE_RESULT,
  TAG_BISECT_INPUT,
  TAG_BISECT_RESULT,
  TAG_FINDROOTS_INPUT,
  TAG_FINDROOTS_RESULT
};

#endif // MW_MPI_IFACE_H
