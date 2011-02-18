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

#include <boost/foreach.hpp>

MPIMaster::MPIMaster(const mpi::communicator& comm)
: m_comm(comm)
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
  m_responses.resize(comm.size());
}

MPIMaster::~MPIMaster()
{
  BOOST_FOREACH(int slave, m_slaves) {
    m_comm.isend(slave, TAG_CMD, CMD_CLOSE);
  }
}

void MPIMaster::startBisect(const fp from, const fp to, const fp mwFreqGHz)
{
  // notify slaves about work conditions

  // diagonalize the hamiltonian at the edges
  {
    m_bisectNodes.resize(2);
    vector<mpi::request> diagRequests(2);

    for (int i = 0; i < 2; ++i) {
      m_comm.isend(m_slaves.at(i), TAG_CMD, CMD_DIAGONALIZE);

      m_comm.isend(m_slaves.at(i), TAG_DIAGONALIZE_INPUT, i == 0 ? from : to);

      diagRequests[i] = m_comm.irecv(m_slaves.at(i), TAG_DIAGONALIZE_RESULT, m_bisectNodes.at(i));
    }

    mpi::wait_all(diagRequests.begin(), diagRequests.end());
  }

  m_pendingSegments.push_back(BRange(from, to));
  // dispatch commands to nodes
  while(!m_pendingSegments.empty() || !m_pendingRequests.empty()) {
    if (!m_pendingRequests.empty()) {
      // check for finished requests
      while(true) {
        boost::optional<ResponsePair> status = mpi::test_any(m_pendingRequests.begin(), m_pendingRequests.end());
        if (status) {
          handleBisectResponse(*status);
        } else {
          break;
        }
      }
    }
    if (m_availableSlaves.empty()) {
      cout << "all slaves working, waiting for any to finish before assigning new work..." << endl;
      ResponsePair status = mpi::wait_any(m_pendingRequests.begin(), m_pendingRequests.end());
      handleBisectResponse(status);
    }
    if (!m_pendingSegments.empty()) {
      const BRange segment = m_pendingSegments.back();
      m_pendingSegments.pop_back();
      int slave = m_availableSlaves.back();
      m_availableSlaves.pop_back();
      cout << "assigning bisect work to slave " << slave << " for range " << segment.first << " to " << segment.second << endl;
      m_comm.isend(slave, TAG_CMD, CMD_BISECT);
      m_comm.isend(slave, TAG_BISECT_INPUT, BisectInput(segment.first, segment.second, mwFreqGHz));
      m_pendingRequests.push_back(m_comm.irecv(slave, TAG_BISECT_RESULT, m_responses.at(slave)));
    }
  }

  m_bisectNodes.clear();
}

void MPIMaster::handleBisectResponse(ResponsePair response)
{
  int slave = response.first.source();
  m_availableSlaves.push_back(slave);
  BisectAnswer answer = m_responses.at(slave);
  cout << "got bisect answer:" << slave << '\t' << answer.from << '\t' << answer.mid << '\t' << answer.to << endl;
  switch (answer.status) {
    case BisectAnswer::Continue:
      m_pendingSegments.push_back(BRange(answer.from, answer.mid));
      m_pendingSegments.push_back(BRange(answer.mid, answer.to));
      break;
    case BisectAnswer::Resonant:
      m_resonantSegments.push_back(BRange(answer.from, answer.mid));
      m_resonantSegments.push_back(BRange(answer.mid, answer.to));
      break;
    case BisectAnswer::NotResonant:
      break;
  }
  m_pendingRequests.erase(response.second);
}



