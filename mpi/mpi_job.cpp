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
using namespace boost::archive;

#define DEBUG_FUNC(args)
// #define DEBUG_FUNC(args) cout << this << ' ' << __PRETTY_FUNCTION__ << '\t' args << endl;

//BEGIN MPIJob

MPIJob::MPIJob(MPIMaster* master)
: m_master(master), m_id(0)
{

}

MPIJob::~MPIJob()
{

}

void MPIJob::setJobId(unsigned int id)
{
  m_id = id;
}

unsigned int MPIJob::jobId() const
{
  return m_id;
}

//BEGIN BisectStartJob
BisectStartJob::BisectStartJob(MPIMaster* master, const fp from, const fp to)
: MPIJob(master)
, m_from(from)
, m_to(to)
{
  DEBUG_FUNC(<< from << '\t' << to)
}

BisectStartJob::~BisectStartJob()
{
  DEBUG_FUNC(<< m_from << '\t' << m_to)
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
  const fp stepSize = (m_to - m_from) / (slavesToUse-1);
  fp B = m_from;
  for(int i = 0; i < slavesToUse; ++i) {
    int slave = m_master->runCommand(this, CMD_DIAGONALIZE,
                                     TAG_DIAGONALIZE_INPUT, B,
                                     TAG_DIAGONALIZE_RESULT, m_results.at(i));
    DEBUG_FUNC(<< slave << '\t' << B)
    m_slaveToResult[slave] = i;
    B += stepSize;
  }
}

void BisectStartJob::handleResult(const int slave)
{
  DEBUG_FUNC(<< slave)

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

void BisectStartJob::saveTo(binary_oarchive& archive) const
{
  archive << m_from << m_to;
}

MPIJob* BisectStartJob::constructFrom(boost::archive::binary_iarchive& archive, MPIMaster* master)
{
  fp from, to;
  archive >> from >> to;
  return new BisectStartJob(master, from, to);
}

MPIJob::JobType BisectStartJob::type() const
{
  return BisectStart;
}

//BEGIN BisectJob

BisectJob::BisectJob(MPIMaster* master, const BisectNode& from, const BisectNode& to)
: MPIJob(master)
, m_from(from)
, m_to(to)
{
  DEBUG_FUNC(<< from.B << '\t' << to.B)
}

BisectJob::~BisectJob()
{
  DEBUG_FUNC(<< m_from.B << '\t' << m_to.B)
}

void BisectJob::start()
{
  int slave = m_master->runCommand(this, CMD_BISECT,
                                   TAG_BISECT_INPUT, BisectInput(m_from, m_to),
                                   TAG_BISECT_RESULT, m_answer);

  DEBUG_FUNC(<< slave << '\t' << m_from.B << '\t' << m_to.B)
}

void BisectJob::handleResult(const int slave)
{
  DEBUG_FUNC(<< slave << '\t' << m_from.B << '\t' << m_to.B << '\t' << m_answer.status)
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

void BisectJob::saveTo(binary_oarchive& archive) const
{
  archive << m_from << m_to;
}

MPIJob* BisectJob::constructFrom(boost::archive::binary_iarchive& archive, MPIMaster* master)
{
  BisectNode from, to;
  archive >> from >> to;
  return new BisectJob(master, from, to);
}

MPIJob::JobType BisectJob::type() const
{
  return Bisect;
}

//BEGIN FindRootsJob

FindRootsJob::FindRootsJob(MPIMaster* master, const BisectNode& from, const BisectNode& to)
: MPIJob(master)
, m_from(from)
, m_to(to)
{
  DEBUG_FUNC(<< from.B << '\t' << to.B)
}

FindRootsJob::~FindRootsJob()
{
  DEBUG_FUNC(<< m_from.B << '\t' << m_to.B)
}

void FindRootsJob::start()
{
  int slave = m_master->runCommand(this, CMD_FINDROOTS,
                                   TAG_FINDROOTS_INPUT, BisectInput(m_from, m_to),
                                   TAG_FINDROOTS_RESULT, m_answer);

  DEBUG_FUNC(<< slave << '\t' << m_from.B << '\t' << m_to.B)
}

void FindRootsJob::handleResult(const int slave)
{
  DEBUG_FUNC(<< slave << '\t' << m_from.B << '\t' << m_to.B << '\t' << m_answer.size())

  ResonanceField::cleanupResonancyField(m_answer);
  for(int i = 0, c = m_answer.size(); i < c; ++i) {
    m_master->enqueueJob(new IntensityJob(m_master, m_answer.at(i)));
  }
}

void FindRootsJob::saveTo(binary_oarchive& archive) const
{
  archive << m_from << m_to;
}

MPIJob* FindRootsJob::constructFrom(boost::archive::binary_iarchive& archive, MPIMaster* master)
{
  BisectNode from, to;
  archive >> from >> to;
  return new FindRootsJob(master, from, to);
}

MPIJob::JobType FindRootsJob::type() const
{
  return FindRoots;
}

//BEGIN IntensityJob

IntensityJob::IntensityJob(MPIMaster* master, const fp B)
: MPIJob(master)
, m_B(B)
{
  DEBUG_FUNC(<< B)
}

IntensityJob::~IntensityJob()
{
  DEBUG_FUNC(<< m_B)
}

void IntensityJob::start()
{
  int slave = m_master->runCommand(this, CMD_INTENSITY,
                                   TAG_INTENSITY_INPUT, m_B,
                                   TAG_INTENSITY_RESULT, m_answer);

  DEBUG_FUNC(<< slave << m_B)
}

void IntensityJob::handleResult(const int slave)
{
  DEBUG_FUNC(<< slave << '\t' << m_B << '\t' << m_answer)
  m_master->intensityOutputFile() << m_B << '\t' << m_answer << endl;
}

void IntensityJob::saveTo(binary_oarchive& archive) const
{
  archive << m_B;
}

MPIJob* IntensityJob::constructFrom(boost::archive::binary_iarchive& archive, MPIMaster* master)
{
  fp B;
  archive >> B;
  return new IntensityJob(master, B);
}

MPIJob::JobType IntensityJob::type() const
{
  return Intensity;
}
