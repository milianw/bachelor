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

#include "mpi_job.h"

#include "mpi_master.h"

#include <iostream>

using namespace std;

#define DEBUG_FUNC
// #define DEBUG_FUNC cout << __FUNCTION__ << '\t' << this << endl;

//BEGIN MPIJob

MPIJob::MPIJob(MPIMaster* master)
: m_master(master)
{

}

MPIJob::~MPIJob()
{

}

//BEGIN BisectStartJob
BisectStartJob::BisectStartJob(MPIMaster* master, const fp from, const fp to)
: MPIJob(master)
, m_from(from)
, m_to(to)
{
  DEBUG_FUNC
}

BisectStartJob::~BisectStartJob()
{
  DEBUG_FUNC
}

void BisectStartJob::start()
{
  const int slavesToUse = m_master->availableSlaves();
  if (slavesToUse < 2) {
    cerr << "need at least two slaves to work properly" << endl;
    exit(1);
  }

  // use as many slaves as possible from the start
  m_results.resize(slavesToUse);
  const fp stepSize = (m_to - m_from) / slavesToUse;
  fp B = m_from;
  for(int i = 0; i < slavesToUse; ++i) {
    int slave = m_master->runCommand(this, CMD_DIAGONALIZE,
                                     TAG_DIAGONALIZE_INPUT, B,
                                     TAG_DIAGONALIZE_RESULT, m_results.at(i));
    m_slaveToResult[slave] = i;
    B += stepSize;
  }
}

void BisectStartJob::handleResult(const int slave)
{
  DEBUG_FUNC

  const int idx = m_slaveToResult.at(slave);

  // check previous result
  if (idx > 0 && m_results.at(idx - 1).E.size()) {
    // previous has also finished, we can bisect this segment further
    m_master->enqueueJob(new BisectJob(m_master, m_results.at(idx - 1), m_results.at(idx)));
  }
  // check next result
  if (idx + 1 < m_results.size() && m_results.at(idx + 1).E.size()) {
    // next has also finished, we can bisect this segment further
    m_master->enqueueJob(new BisectJob(m_master, m_results.at(idx), m_results.at(idx + 1)));
  }
}

//BEGIN BisectJob

BisectJob::BisectJob(MPIMaster* master, const BisectNode& from, const BisectNode& to)
: MPIJob(master)
, m_from(from)
, m_to(to)
{
  DEBUG_FUNC
}

BisectJob::~BisectJob()
{
  DEBUG_FUNC
}

void BisectJob::start()
{
  m_master->runCommand(this, CMD_BISECT,
                       TAG_BISECT_INPUT, BisectInput(m_from, m_to),
                       TAG_BISECT_RESULT, m_answer);
}

void BisectJob::handleResult(const int /*slave*/)
{
  switch (m_answer.status) {
    case BisectAnswer::Continue:
      m_master->enqueueJob(new BisectJob(m_master, m_from, m_answer.mid));
      m_master->enqueueJob(new BisectJob(m_master, m_answer.mid, m_to));
      break;
    case BisectAnswer::Resonant:
      m_master->enqueueJob(new FindRootsJob(m_master, m_from, m_answer.mid));
      m_master->enqueueJob(new FindRootsJob(m_master, m_answer.mid, m_to));
      break;
    case BisectAnswer::NotResonant:
      // nothing to do
      break;
  }
}

//BEGIN FindRootsJob

FindRootsJob::FindRootsJob(MPIMaster* master, const BisectNode& from, const BisectNode& to)
: MPIJob(master)
, m_from(from)
, m_to(to)
{
  DEBUG_FUNC
}

FindRootsJob::~FindRootsJob()
{
  DEBUG_FUNC
}

void FindRootsJob::start()
{
  m_master->runCommand(this, CMD_FINDROOTS,
                       TAG_FINDROOTS_INPUT, BisectInput(m_from, m_to),
                       TAG_FINDROOTS_RESULT, m_answer);
}

void FindRootsJob::handleResult(const int /*slave*/)
{
  ResonanceField::cleanupResonancyField(m_answer);
  for(int i = 0, c = m_answer.size(); i < c; ++i) {
    m_master->enqueueJob(new IntensityJob(m_master, m_answer.at(i)));
  }
}

//BEGIN IntensityJob

IntensityJob::IntensityJob(MPIMaster* master, const fp B)
: MPIJob(master)
, m_B(B)
{
  DEBUG_FUNC
}

IntensityJob::~IntensityJob()
{
  DEBUG_FUNC
}

void IntensityJob::start()
{
  m_master->runCommand(this, CMD_INTENSITY,
                       TAG_INTENSITY_INPUT, m_B,
                       TAG_INTENSITY_RESULT, m_answer);
}

void IntensityJob::handleResult(const int /*slave*/)
{
  m_master->intensityOutputFile() << m_B << '\t' << m_answer << endl;
}
