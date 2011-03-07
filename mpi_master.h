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

#include <boost/function.hpp>

#include <fstream>

class Experiment;

/**
 * MPI master node that delegates work to slave nodes
 */
class MPIMaster {
public:
  MPIMaster(const mpi::communicator& comm, const Experiment& exp,
            const string& outputDir);
  ~MPIMaster();

  void startBisect(const fp from, const fp to);

private:
  const mpi::communicator& m_comm;
  const Experiment& m_exp;
  vector<int> m_slaves;
  vector<int> m_availableSlaves;
  const string& m_outputDir;

  typedef pair<fp, fp> BRange;
  vector<BRange> m_pendingSegments;
  vector<BisectAnswer> m_bisectResponses;
  vector<BRange> m_resonantSegments;
  map<fp, BisectNode> m_bisectNodes;

  vector< vector<fp> > m_findRootResponses;
  vector<fp> m_resonancyField;

  string m_intensityOutputFile;
  ofstream m_intensityOutput;
  vector<IntensityAnswer> m_intensityResponses;

  typedef vector<mpi::request> RequestList;
  RequestList m_pendingRequests;
  typedef pair<mpi::status, RequestList::iterator> ResponsePair;
  typedef boost::function<void(int)> ResponseHandler;
  void checkResponses(ResponseHandler handler);
  void handleResponseGeneric(ResponseHandler handler, const ResponsePair& response);
  void handleBisectResponse(int slave);
  void handleFindRootResponse(int slave);
  void handleIntensityResponse(int slave);
};

#endif // MW_MPI_SERVER_H
