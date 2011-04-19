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
#include "spinlib/helpers.h"

#include <boost/foreach.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <typeinfo>
#include <sstream>

using namespace std;
namespace ar = boost::archive;
namespace fs = boost::filesystem;
namespace pt = boost::posix_time;

MPIMaster::MPIMaster(const mpi::communicator& comm, const Experiment& exp,
                     const string& outputDir)
: m_comm(comm)
, m_exp(exp)
, m_outputDir(outputDir)
, m_lastJobId(0)
{
  if (comm.rank() != MASTER_RANK) {
    cerr << timeStamp() << "MPIMaster created outside of master rank!" << endl;
    comm.abort(1);
    return;
  }
  if (comm.size() - 1 < 2) {
    cerr << timeStamp() << "need at least two slaves to operate properly" << endl;
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
  cout << timeStamp() << "intensity data written to file:" << endl << m_intensityOutputFile  << endl;
}

void MPIMaster::calculateIntensity(const fp from, const fp to, const int steps, const std::vector< Orientation >& orientations)
{
  bool continuingJobs = false;
  {
  fs::path outputPath(m_outputDir);
  fs::directory_iterator it(outputPath);
  fs::directory_iterator end;
  while(it != end) {
    if (fs::is_regular_file(it->status())) {
      if (boost::ends_with(it->path().filename(), ".job")) {
        m_lastJobId = max(boost::lexical_cast<unsigned int>(it->path().stem()), m_lastJobId);
        bool added = readdJob(it->path().file_string());
        continuingJobs = continuingJobs || added;
      }
    }
    ++it;
  }
  }

  if (!continuingJobs) {
    if (steps > 0) {
      const fp stepSize = (to - from) / steps;
      BOOST_FOREACH(const Orientation& orientation, orientations) {
        fp B = from;
        for(int i = 0; i < steps; ++i) {
          enqueueJob(new IntensityJob(this, IntensityInput(B, orientation)));
          B += stepSize;
        }
      }
    } else {
      BOOST_FOREACH(const Orientation& orientation, orientations) {
        enqueueJob(new BisectStartJob(this, from, to, orientation));
      }
    }
    cout << timeStamp() << "starting jobs from scratch" << endl;
  } else {
    cout << timeStamp() << "continuing with " << m_jobQueue.size() << " jobs from last run" << endl;
  }

  m_intensityOutputFile = m_outputDir + "/intensity.data";
  // truncate data file if we start from scratch
  m_intensityOutput.open(m_intensityOutputFile.data(), ios_base::out | (continuingJobs ? ios_base::app : ios_base::trunc));
  if (!m_intensityOutput.is_open()) {
    cerr << timeStamp() << "could not open output file" << m_intensityOutputFile << endl;
    m_comm.abort(3);
  }

  bool forceWaitJob = false;
  while(!m_jobQueue.empty() || !m_pendingRequests.empty()) {
    cout << timeStamp() << "QUEUE: free slaves: " << m_availableSlaves.size() << ", pending jobs: " << m_jobQueue.size() << endl;

    // check for finished requests
    while(!m_pendingRequests.empty()) {
      boost::optional<ResponsePair> status = mpi::test_any(m_pendingRequests.begin(), m_pendingRequests.end());
      if (status) {
        handleResponse(*status);
        forceWaitJob = false;
      } else {
        break;
      }
    }

    if ((m_availableSlaves.empty() || m_jobQueue.empty() || forceWaitJob) && !m_pendingRequests.empty()) {
//       cout << timeStamp() << "all slaves working, waiting for any to finish before assigning new work..." << endl;
      handleResponse(mpi::wait_any(m_pendingRequests.begin(), m_pendingRequests.end()));
      forceWaitJob = false;
    }

    while (!m_jobQueue.empty() && !m_availableSlaves.empty()) {
      std::vector<MPIJob*>::reverse_iterator it = m_jobQueue.rbegin();
      bool foundJob = false;
      while(it != m_jobQueue.rend()) {
        if (dynamic_cast<BisectStartJob*>(*it) && m_availableSlaves.size() < 2) {
          ///TODO: put this into the MPIJob API?
          ++it;
          continue;
        }
        foundJob = true;
        break;
      }
      if (!foundJob) {
        cout << timeStamp() << "not enough slaves to continue jobs, waiting" << endl;
        forceWaitJob = true;
        break;
      }
      MPIJob* job = *it;
      cout << timeStamp() << "START: " << job->name() << endl;
      m_jobQueue.erase(--it.base());
      job->start();
    }
  }
}

string MPIMaster::jobFile(unsigned int jobId) const
{
  stringstream stream;
  stream << m_outputDir << '/' << jobId << ".job";
  return stream.str();
}

void MPIMaster::enqueueJob(MPIJob* job)
{
  unsigned int id = ++m_lastJobId;
  job->setJobId(id);

  // make jobs restartable, save their data
  const string file = jobFile(id);
  ofstream stream(file.c_str());
  if (!stream.is_open()) {
    cerr << timeStamp() << "could not open job file for writing: " << file << endl;
    m_comm.abort(4);
    return;
  }
  ar::binary_oarchive archive(stream);
  int type = static_cast<int>(job->type());
  archive << type;
  archive << id;
  job->saveTo(archive);

  m_jobQueue.push_back(job);
}

bool MPIMaster::readdJob(const std::string& jobFile)
{
  ifstream stream(jobFile.c_str());
  if (!stream.is_open()) {
    cerr << timeStamp() << "could not open job file for reading, skipping: " << jobFile << endl;
    return false;
  }

  try {
    ar::binary_iarchive archive(stream);
    int type;
    unsigned int id;
    archive >> type;
    archive >> id;

    MPIJob* job = 0;
    switch(static_cast<MPIJob::JobType>(type)) {
      case MPIJob::BisectStart:
        job = BisectStartJob::constructFrom(archive, this);
        break;
      case MPIJob::Bisect:
        job = BisectJob::constructFrom(archive, this);
        break;
      case MPIJob::FindRoots:
        job = FindRootsJob::constructFrom(archive, this);
        break;
      case MPIJob::Intensity:
        job = IntensityJob::constructFrom(archive, this);
        break;
    }

    // reuse last id
    job->setJobId(id);
    m_lastJobId = max(id, m_lastJobId);

    m_jobQueue.push_back(job);
    return true;
  } catch(boost::archive::archive_exception e) {
    cerr << timeStamp() << "could not load job file, skipping:" << jobFile << endl << e.what() << endl;
    return false;
  }
}

void MPIMaster::handleResponse(const ResponsePair& response)
{
  m_pendingRequests.erase(response.second);

  const int slave = response.first.source();
  m_availableSlaves.push_back(slave);

  m_sendRequests.erase(slave);

  MPIJob* job = m_runningJobs.at(slave);
  cout << timeStamp() << "RESPONSE: " << job->name() << ", slave:" << slave << endl;
  job->handleResult(slave);
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
    if (!fs::remove(jobFile(job->jobId()))) {
      cerr << timeStamp() << "could not remove job file: " << jobFile(job->jobId()) << endl;
    }
    delete job;
  }
//   cout << timeStamp() << "slave " << slave << " finished work on job " << job << (stillRunning ? ", job still running" : "") << endl;
}

ofstream& MPIMaster::intensityOutputFile()
{
  return m_intensityOutput;
}

std::string MPIMaster::timeStamp() const
{
  stringstream stream;
  stream << '[' << pt::second_clock::local_time() << "] ";
  return stream.str();
}
