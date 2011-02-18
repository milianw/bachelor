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

#ifndef MW_MPI_MASTER_H
#define MW_MPI_MASTER_H

#include "types.h"

/**
 * MPI master node that delegates work to slave nodes
 */
class MPIMaster {
public:
  MPIMaster(const mpi::communicator& comm);
  ~MPIMaster();

  void startBisect(const fp from, const fp to, const fp mwFreqGHz);

private:
  const mpi::communicator& m_comm;
  vector<int> m_slaves;
  vector<int> m_availableSlaves;

  typedef pair< mpi::status, vector< mpi::request >::iterator> ResponsePair;
  void handleBisectResponse(ResponsePair response);

  typedef pair<fp, fp> BRange;
  vector<BRange> m_pendingSegments;
  vector<mpi::request> m_pendingRequests;
  vector<BisectAnswer> m_responses;
  vector<BRange> m_resonantSegments;
  vector<BisectNode> m_bisectNodes;
};

#endif // MW_MPI_SERVER_H
