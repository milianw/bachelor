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
#include <vector>

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

  void calculateIntensity(const fp from, const fp to, const int steps, const std::vector<Orientation>& orientations);

  template<typename InputT, typename OutputT>
  int runCommand(MPIJob* job, Commands cmd,
                 Tags inputTag, const InputT& input,
                 Tags outputTag, OutputT& output);

  void enqueueJob(MPIJob *job);

  std::ofstream& intensityOutputFile();

  std::string timeStamp() const;

private:
  std::string jobFile(unsigned int jobId) const;
  bool readdJob(const std::string& jobFile);

  const mpi::communicator& m_comm;
  const Experiment& m_exp;
  std::vector<int> m_slaves;
  std::vector<int> m_availableSlaves;
  const std::string& m_outputDir;

  std::vector<MPIJob *> m_jobQueue;
  unsigned int m_lastJobId;

  std::string m_intensityOutputFile;
  std::ofstream m_intensityOutput;

  // slave <-> job
  std::map<int, MPIJob *> m_runningJobs;

  // running requests
  typedef std::vector<mpi::request> RequestList;
  RequestList m_pendingRequests;
  typedef std::pair<mpi::request, mpi::request> SendRequestPair;
  std::map<int, SendRequestPair> m_sendRequests;
  typedef std::pair<mpi::status, RequestList::iterator> ResponsePair;
  void handleResponse(const ResponsePair& response);
};

template<typename InputT, typename OutputT>
int MPIMaster::runCommand(MPIJob* job, Commands cmd,
                          Tags inputTag, const InputT& input,
                          Tags outputTag, OutputT& output)
{
  const int slave = m_availableSlaves.back();
  std::cout << timeStamp() << "RUN: " << stringifyCommand(cmd) << ", slave: " << slave << std::endl;
  m_availableSlaves.pop_back();
  m_runningJobs[slave] = job;

  ///NOTE: if we don't store the isend requests, we can/will encounter crashes for custom (non-MPI) data types
  ///      since the parameter passed to the MPI call is actually a reference to a member of that request.
  ///      Hence we need to store it and clean it up afterwards!
  /// reproducible on my machine with: mpirun -np 3 ./mpi/hs-mpi -p 5 -n 1 -f 0 -t 1 -m 9.5
  mpi::request cmdReq = m_comm.isend(slave, TAG_CMD, cmd);
  mpi::request inputReq = m_comm.isend(slave, inputTag, input);
  m_sendRequests[slave] = SendRequestPair(cmdReq, inputReq);
  m_pendingRequests.push_back(m_comm.irecv(slave, outputTag, output));

  return slave;
}

#endif // MW_MPI_MASTER_H
