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

#include "mpi_slave.h"

#include "mpi_iface.h"

MPISlave::MPISlave(const boost::mpi::communicator& comm)
: m_comm(comm)
{

}

MPISlave::~MPISlave()
{

}

void MPISlave::work()
{
  int cmd;
  while(true) {
    m_comm.recv(MASTER_RANK, TAG_CMD, cmd);
    if (cmd == CMD_CLOSE) {
      break;
    } else if (cmd == CMD_BISECT) {
      BisectInput input;
      m_comm.recv(MASTER_RANK, TAG_BISECT_INPUT, input);
      if (abs(input.from - input.to) > 0.001) {
        m_comm.isend(MASTER_RANK, TAG_BISECT_RESULT, BisectAnswer(BisectAnswer::Continue, input.from, (input.from + input.to)/2, input.to));
      } else {
        m_comm.isend(MASTER_RANK, TAG_BISECT_RESULT, BisectAnswer(BisectAnswer::NotResonant, input.from, (input.from + input.to)/2, input.to));
      }
    }
  }
}
