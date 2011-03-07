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

#include "spinhamiltonian.h"
#include "experiment.h"

MPISlave::MPISlave(const mpi::communicator& comm, const Experiment& exp)
: m_comm(comm)
, m_exp(exp)
{
  if (comm.rank() == MASTER_RANK) {
    cerr << "MPISlave created in master rank!" << endl;
    comm.abort(1);
    return;
  }

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
    } else if (cmd == CMD_DIAGONALIZE) {
      fp B;
      m_comm.recv(MASTER_RANK, TAG_DIAGONALIZE_INPUT, B);
      SpinHamiltonian H(B, m_exp);
      SelfAdjointEigenSolver<MatrixXc> eigenSolver(H.hamiltonian());
      VectorX E = eigenSolver.eigenvalues();
      const MatrixXc eigenVectors = eigenSolver.eigenvectors();
      const MatrixXc G = H.nuclearZeeman() + H.electronZeeman();
      VectorX E_deriv(m_exp.dimension);
      for(int u = 0; u < m_exp.dimension; ++u) {
        // <u| G |u> => expectation value is always real
        E_deriv(u) = (eigenVectors.col(u).adjoint() * G * eigenVectors.col(u))(0, 0).real();
      }
      m_comm.isend(MASTER_RANK, TAG_DIAGONALIZE_RESULT, BisectNode(B, E, E_deriv));
    }
  }
}
