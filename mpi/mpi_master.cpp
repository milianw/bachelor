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

#include "mpi_master.h"

#include "mpi_job.h"

#include "spinlib/experiment.h"

#include <boost/foreach.hpp>

#include <sstream>

using namespace std;

string intensityOutputFile(const Experiment& exp, const string& outputDir, const fp from, const fp to)
{
  stringstream stream;
  stream << outputDir << '/' << exp.nProtons << ':' << from << '-' << to << ":auto:" << exp.mwFreqGHz << ":mpi";
  return stream.str();
}

MPIMaster::MPIMaster(const mpi::communicator& comm, const Experiment& exp,
                     const string& outputDir)
: m_comm(comm)
, m_exp(exp)
, m_outputDir(outputDir)
{
  if (comm.rank() != MASTER_RANK) {
    cerr << "MPIMaster created outside of master rank!" << endl;
    comm.abort(1);
    return;
  }
  if (comm.size() - 1 < 2) {
    cerr << "need at least two slaves to operate properly" << endl;
    comm.abort(2);
    return;
  }
  for (int i = 0; i < comm.size(); ++i) {
    if (i != MASTER_RANK) {
      m_slaves.push_back(i);
    }
  }
  m_availableSlaves = m_slaves;
}

MPIMaster::~MPIMaster()
{
  BOOST_FOREACH(int slave, m_slaves) {
    m_comm.isend(slave, TAG_CMD, CMD_CLOSE);
  }
  cout << "intensity data written to file:" << endl << m_intensityOutputFile  << endl;
}

void MPIMaster::startBisect(const fp from, const fp to)
{
  m_intensityOutputFile = ::intensityOutputFile(m_exp, m_outputDir, from, to);
  m_intensityOutput.open(m_intensityOutputFile.data());

  enqueueJob(new BisectStartJob(this, from, to));

  while(!m_jobQueue.empty() || !m_pendingRequests.empty()) {
    // check for finished requests
    while(!m_pendingRequests.empty()) {
      boost::optional<ResponsePair> status = mpi::test_any(m_pendingRequests.begin(), m_pendingRequests.end());
      if (status) {
        handleResponse(*status);
      } else {
        break;
      }
    }

    if (m_availableSlaves.empty()) {
//       cout << "all slaves working, waiting for any to finish before assigning new work..." << endl;
      handleResponse(mpi::wait_any(m_pendingRequests.begin(), m_pendingRequests.end()));
    }

    if (!m_jobQueue.empty()) {
      MPIJob* job = m_jobQueue.front();
      m_jobQueue.pop();
      job->start();
    }
  }
}

void MPIMaster::enqueueJob(MPIJob* job)
{
  m_jobQueue.push(job);
}

void MPIMaster::handleResponse(const ResponsePair& response)
{
  m_pendingRequests.erase(response.second);

  const int slave = response.first.source();
  m_availableSlaves.push_back(slave);

  MPIJob* job = m_runningJobs.at(slave);
  job->handleResult();
  m_runningJobs.erase(slave);

  bool stillRunning = false;
  std::map< int, MPIJob* >::const_iterator it = m_runningJobs.begin();
  std::map< int, MPIJob* >::const_iterator end = m_runningJobs.end();
  while(it != end) {
    if (it->second == job) {
      stillRunning = true;
      break;
    }
    ++it;
  }
  if (!stillRunning) {
    delete job;
  }
//   cout << "slave " << slave << " finished work on job " << job << (stillRunning ? ", job still running" : "") << endl;
}

ofstream& MPIMaster::intensityOutputFile()
{
  return m_intensityOutput;
}
