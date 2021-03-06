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
#ifndef MW_MPI_SLAVE_H
#define MW_MPI_SLAVE_H

#include "mpi_types.h"

#include "spinlib/experiment.h"

class MPISlave {
public:
  MPISlave(const mpi::communicator& comm, const Experiment& exp);
  ~MPISlave();

  void work();

private:
  const mpi::communicator& m_comm;
  Experiment m_exp;
  const ResonanceField m_resonanceField;
};

#endif // MW_MPI_SLAVE_H
