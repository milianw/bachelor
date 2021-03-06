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

#include "mpi_iface.h"

const char* stringifyCommand(Commands cmd)
{
  switch(cmd) {
    case CMD_CLOSE:
      return "Close";
    case CMD_BISECT:
      return "Bisect";
    case CMD_DIAGONALIZE:
      return "Diagonalize";
    case CMD_FINDROOTS:
      return "Findroots";
    case CMD_INTENSITY:
      return "Intensity";
  }
  return "UNKNOWN COMMAND";
}
