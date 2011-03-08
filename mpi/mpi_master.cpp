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

#include "mpi_iface.h"

#include "spinlib/experiment.h"

#include <boost/foreach.hpp>
#include <boost/bind.hpp>

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
  /// find resonant segments
  // diagonalize the hamiltonian at the edges
  {
    m_bisectNodes[from] = BisectNode();
    m_bisectNodes[to] = BisectNode();
    vector<mpi::request> diagRequests(2);

    for (int i = 0; i < 2; ++i) {
      const fp B = (i == 0) ? from : to;
      m_comm.isend(m_slaves.at(i), TAG_CMD, CMD_DIAGONALIZE);

      m_comm.isend(m_slaves.at(i), TAG_DIAGONALIZE_INPUT, B);

      diagRequests[i] = m_comm.irecv(m_slaves.at(i), TAG_DIAGONALIZE_RESULT, m_bisectNodes[B]);
    }

    mpi::wait_all(diagRequests.begin(), diagRequests.end());
  }

  m_pendingSegments.push_back(BRange(from, to));
  // dispatch commands to nodes
  m_bisectResponses.resize(m_comm.size());
  while(!m_pendingSegments.empty() || !m_pendingRequests.empty()) {
    checkResponses(boost::bind(&MPIMaster::handleBisectResponse, this, _1));

    if (!m_pendingSegments.empty()) {
      const BRange segment = m_pendingSegments.back();
      m_pendingSegments.pop_back();
      int slave = m_availableSlaves.back();
      m_availableSlaves.pop_back();
      cout << "assigning bisect work to slave " << slave << " for range " << segment.first << " to " << segment.second << endl;
      m_comm.isend(slave, TAG_CMD, CMD_BISECT);
      m_comm.isend(slave, TAG_BISECT_INPUT, BisectInput(m_bisectNodes.at(segment.first), m_bisectNodes.at(segment.second)));
      m_pendingRequests.push_back(m_comm.irecv(slave, TAG_BISECT_RESULT, m_bisectResponses.at(slave)));
    }
  }
//   m_bisectResponses.clear();
  ///TODO: clear unneeded parts of m_bisectNodes

  if (m_resonantSegments.empty()) {
    cerr << "could not find any resonant segments in B range " << from << " to " << to << "T" << endl;
    return;
  }

  cout << "found resonant segments:" << m_resonantSegments.size() << endl;

  /// find roots in resonant segments
  m_findRootResponses.resize(m_comm.size());
  while(!m_resonantSegments.empty() || !m_pendingRequests.empty()) {
    // check for finished requests
    checkResponses(boost::bind(&MPIMaster::handleFindRootResponse, this, _1));

    if (!m_resonantSegments.empty()) {
      const BRange segment = m_resonantSegments.back();
      m_resonantSegments.pop_back();
      int slave = m_availableSlaves.back();
      m_availableSlaves.pop_back();
      cout << "assigning root-finding work to slave " << slave << " for range " << segment.first << " to " << segment.second << endl;
      m_comm.isend(slave, TAG_CMD, CMD_FINDROOTS);
      m_comm.isend(slave, TAG_FINDROOTS_INPUT, BisectInput(m_bisectNodes.at(segment.first), m_bisectNodes.at(segment.second)));
      m_pendingRequests.push_back(m_comm.irecv(slave, TAG_FINDROOTS_RESULT, m_findRootResponses.at(slave)));
    }
  }
  m_bisectNodes.clear();
  m_findRootResponses.clear();

  ResonanceField::cleanupResonancyField(m_resonancyField);

  cout << "found resonancy field:" << m_resonancyField.size() << endl;

  /// calculate intensity for found roots
  m_intensityOutputFile = intensityOutputFile(m_exp, m_outputDir, from, to);
  m_intensityOutput.open(m_intensityOutputFile.data());
  m_intensityResponses.resize(m_comm.size());
  while(!m_resonancyField.empty() || !m_pendingRequests.empty()) {
    // check for finished requests
    checkResponses(boost::bind(&MPIMaster::handleIntensityResponse, this, _1));

    if (!m_resonancyField.empty()) {
      const fp B = m_resonancyField.back();
      m_resonancyField.pop_back();
      int slave = m_availableSlaves.back();
      m_availableSlaves.pop_back();
      cout << "assigning intensity-calculation work to slave " << slave << " at B = " << B << "T" << endl;
      m_comm.isend(slave, TAG_CMD, CMD_INTENSITY);
      m_comm.isend(slave, TAG_INTENSITY_INPUT, B);
      m_pendingRequests.push_back(m_comm.irecv(slave, TAG_INTENSITY_RESULT, m_intensityResponses.at(slave)));
    }
  }
  m_intensityResponses.clear();
}

void MPIMaster::checkResponses(ResponseHandler handler)
{
  // check for finished requests
  while(!m_pendingRequests.empty()) {
    boost::optional<ResponsePair> status = mpi::test_any(m_pendingRequests.begin(), m_pendingRequests.end());
    if (status) {
      handleResponseGeneric(handler, *status);
    } else {
      break;
    }
  }

  if (m_availableSlaves.empty()) {
    cout << "all slaves working, waiting for any to finish before assigning new work..." << endl;
    handleResponseGeneric(handler, mpi::wait_any(m_pendingRequests.begin(), m_pendingRequests.end()));
  }
}

void MPIMaster::handleResponseGeneric(ResponseHandler handler, const ResponsePair& response)
{
  int slave = response.first.source();
  m_availableSlaves.push_back(slave);
  handler(slave);
  m_pendingRequests.erase(response.second);
}

void MPIMaster::handleBisectResponse(int slave)
{
  const BisectAnswer& answer = m_bisectResponses.at(slave);
  switch (answer.status) {
    case BisectAnswer::Continue:
      m_pendingSegments.push_back(BRange(answer.from, answer.mid.B));
      m_pendingSegments.push_back(BRange(answer.mid.B, answer.to));
      m_bisectNodes[answer.mid.B] = answer.mid;
      break;
    case BisectAnswer::Resonant:
      m_resonantSegments.push_back(BRange(answer.from, answer.mid.B));
      m_resonantSegments.push_back(BRange(answer.mid.B, answer.to));
      m_bisectNodes[answer.mid.B] = answer.mid;
      break;
    case BisectAnswer::NotResonant:
      // nothing to do
      break;
  }
}

void MPIMaster::handleFindRootResponse(int slave)
{
  const vector<fp>& roots = m_findRootResponses.at(slave);
  copy(roots.begin(), roots.end(), back_inserter(m_resonancyField));
}

void MPIMaster::handleIntensityResponse(int slave)
{
  const IntensityAnswer& answer = m_intensityResponses.at(slave);
  m_intensityOutput << answer.B << '\t' << answer.intensity << endl;
}
