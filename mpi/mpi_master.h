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

#include "mpi_types.h"

#include <boost/function.hpp>

#include <fstream>

class Experiment;

/**
 * MPI master node that delegates work to slave nodes
 */
class MPIMaster {
public:
  MPIMaster(const mpi::communicator& comm, const Experiment& exp,
            const std::string& outputDir);
  ~MPIMaster();

  void startBisect(const fp from, const fp to);

private:
  const mpi::communicator& m_comm;
  const Experiment& m_exp;
  std::vector<int> m_slaves;
  std::vector<int> m_availableSlaves;
  const std::string& m_outputDir;

  typedef std::pair<fp, fp> BRange;
  std::vector<BRange> m_pendingSegments;
  std::vector<BisectAnswer> m_bisectResponses;
  std::vector<BRange> m_resonantSegments;
  std::map<fp, BisectNode> m_bisectNodes;

  std::vector< std::vector<fp> > m_findRootResponses;
  std::vector<fp> m_resonancyField;

  std::string m_intensityOutputFile;
  std::ofstream m_intensityOutput;
  std::vector<IntensityAnswer> m_intensityResponses;

  typedef std::vector<mpi::request> RequestList;
  RequestList m_pendingRequests;
  typedef std::pair<mpi::status, RequestList::iterator> ResponsePair;
  typedef boost::function<void(int)> ResponseHandler;
  void checkResponses(ResponseHandler handler);
  void handleResponseGeneric(ResponseHandler handler, const ResponsePair& response);
  void handleBisectResponse(int slave);
  void handleFindRootResponse(int slave);
  void handleIntensityResponse(int slave);
};

#endif // MW_MPI_SERVER_H
