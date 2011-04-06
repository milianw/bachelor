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

#include "mpi_iface.h"
#include "mpi_types.h"

#include <fstream>
#include <queue>

class Experiment;
class MPIJob;

/**
 * MPI master node that delegates work to slave nodes
 */
class MPIMaster {
public:
  MPIMaster(const mpi::communicator& comm, const Experiment& exp,
            const std::string& outputDir);
  ~MPIMaster();

  void calculateIntensity(const fp from, const fp to, const int steps);

  template<typename InputT, typename OutputT>
  int runCommand(MPIJob* job, Commands cmd,
                 Tags inputTag, const InputT& input,
                 Tags outputTag, OutputT& output);

  void enqueueJob(MPIJob *job);

  std::ofstream& intensityOutputFile();

  int availableSlaves() const;

private:
  const mpi::communicator& m_comm;
  const Experiment& m_exp;
  std::vector<int> m_slaves;
  std::vector<int> m_availableSlaves;
  const std::string& m_outputDir;

  std::queue<MPIJob *> m_jobQueue;

  std::string m_intensityOutputFile;
  std::ofstream m_intensityOutput;

  // slave <-> job
  std::map<int, MPIJob *> m_runningJobs;

  // running requests
  typedef std::vector<mpi::request> RequestList;
  RequestList m_pendingRequests;
  typedef std::pair<mpi::status, RequestList::iterator> ResponsePair;
  void handleResponse(const ResponsePair& response);
};

template<typename InputT, typename OutputT>
int MPIMaster::runCommand(MPIJob* job, Commands cmd,
                          Tags inputTag, const InputT& input,
                          Tags outputTag, OutputT& output)
{
  const int slave = m_availableSlaves.back();
//   std::cout << "running job " << job << " on slave " << slave << std::endl;
  m_availableSlaves.pop_back();
  m_runningJobs[slave] = job;

  m_comm.isend(slave, TAG_CMD, cmd);
  m_comm.isend(slave, inputTag, input);
  m_pendingRequests.push_back(m_comm.irecv(slave, outputTag, output));

  return slave;
}

#endif // MW_MPI_MASTER_H
