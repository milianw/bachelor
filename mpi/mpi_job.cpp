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

#define UNUSED(x) (void) x

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

string MPIJob::name() const
{
  stringstream stream;
  switch(type()) {
    case BisectStart:
      stream << "BisectStart";
      break;
    case Bisect:
      stream << "Bisect";
      break;
    case FindRoots:
      stream << "FindRoots";
      break;
    case Intensity:
      stream << "Intensity";
      break;
  }
  stream << '#' << m_id;
  return stream.str();
}

//BEGIN BisectStartJob
BisectStartJob::BisectStartJob(MPIMaster* master, const fp from, const fp to, const Orientation& orientation)
: MPIJob(master)
, m_from(from)
, m_to(to)
, m_orientation(orientation)
{
  DEBUG_FUNC(<< m_from << '\t' << m_to << '\t' << m_orientation.orientation.transpose() << '\t' << m_orientation.weight)
}

BisectStartJob::~BisectStartJob()
{
  DEBUG_FUNC(<< m_from << '\t' << m_to << '\t' << m_orientation.orientation.transpose() << '\t' << m_orientation.weight)
}

void BisectStartJob::start()
{
//   const int slavesToUse = m_master->availableSlaves();
  const int slavesToUse = 2;
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
                                     TAG_DIAGONALIZE_INPUT, DiagonalizeInput(B, m_orientation),
                                     TAG_DIAGONALIZE_RESULT, m_results.at(i));
    DEBUG_FUNC(<< slave << '\t' << B << '\t' << m_orientation.orientation.transpose() << '\t' << m_orientation.weight)
    m_slaveToResult[slave] = i;
    B += stepSize;
  }
}

void BisectStartJob::handleResult(const int slave)
{
  DEBUG_FUNC(<< slave)

  const unsigned int idx = m_slaveToResult.at(slave);

  // check previous result
  if (idx > 0 && m_results.at(idx - 1).E.size()) {
    // previous has also finished, we can bisect this segment further
    m_master->enqueueJob(new BisectJob(m_master, BisectInput(m_results.at(idx - 1), m_results.at(idx), m_orientation)));
  }
  // check next result
  if (idx + 1 < m_results.size() && m_results.at(idx + 1).E.size()) {
    // next has also finished, we can bisect this segment further
    m_master->enqueueJob(new BisectJob(m_master, BisectInput(m_results.at(idx), m_results.at(idx + 1), m_orientation)));
  }
}

void BisectStartJob::saveTo(binary_oarchive& archive) const
{
  archive << m_from << m_to << m_orientation;
}

MPIJob* BisectStartJob::constructFrom(boost::archive::binary_iarchive& archive, MPIMaster* master)
{
  fp from, to;
  Orientation orientation;
  archive >> from >> to >> orientation;
  return new BisectStartJob(master, from, to, orientation);
}

MPIJob::JobType BisectStartJob::type() const
{
  return BisectStart;
}

//BEGIN BisectJob

BisectJob::BisectJob(MPIMaster* master, const BisectInput& input)
: MPIJob(master)
, m_input(input)
{
  DEBUG_FUNC(<< m_input.from.B << '\t' << m_input.to.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight)
}

BisectJob::~BisectJob()
{
  DEBUG_FUNC(<< m_input.from.B << '\t' << m_input.to.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight)
}

void BisectJob::start()
{
  int slave = m_master->runCommand(this, CMD_BISECT,
                                   TAG_BISECT_INPUT, m_input,
                                   TAG_BISECT_RESULT, m_answer);

  DEBUG_FUNC(<< slave << m_input.from.B << '\t' << m_input.to.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight)
  UNUSED(slave);
}

void BisectJob::handleResult(const int slave)
{
  DEBUG_FUNC(<< slave << m_input.from.B << '\t' << m_input.to.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight)
  UNUSED(slave);
  switch (m_answer.status) {
    case BisectAnswer::Continue:
      m_master->enqueueJob(new BisectJob(m_master, BisectInput(m_input.from, m_answer.mid, m_input.orientation)));
      m_master->enqueueJob(new BisectJob(m_master, BisectInput(m_answer.mid, m_input.to, m_input.orientation)));
      break;
    case BisectAnswer::Resonant:
      m_master->enqueueJob(new FindRootsJob(m_master, BisectInput(m_input.from, m_answer.mid, m_input.orientation)));
      m_master->enqueueJob(new FindRootsJob(m_master, BisectInput(m_answer.mid, m_input.to, m_input.orientation)));
      break;
    case BisectAnswer::NotResonant:
      // nothing to do
      break;
  }
}

void BisectJob::saveTo(binary_oarchive& archive) const
{
  archive << m_input;
}

MPIJob* BisectJob::constructFrom(boost::archive::binary_iarchive& archive, MPIMaster* master)
{
  BisectInput input;
  archive >> input;
  return new BisectJob(master, input);
}

MPIJob::JobType BisectJob::type() const
{
  return Bisect;
}

//BEGIN FindRootsJob

FindRootsJob::FindRootsJob(MPIMaster* master, const BisectInput& input)
: MPIJob(master)
, m_input(input)
{
  DEBUG_FUNC(<< m_input.from.B << '\t' << m_input.to.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight)
}

FindRootsJob::~FindRootsJob()
{
  DEBUG_FUNC(<< m_input.from.B << '\t' << m_input.to.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight)
}

void FindRootsJob::start()
{
  int slave = m_master->runCommand(this, CMD_FINDROOTS,
                                   TAG_FINDROOTS_INPUT, m_input,
                                   TAG_FINDROOTS_RESULT, m_answer);

  DEBUG_FUNC(<< slave << '\t' << m_input.from.B << '\t' << m_input.to.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight)
  UNUSED(slave);
}

void FindRootsJob::handleResult(const int slave)
{
  DEBUG_FUNC(<< slave << '\t' << m_input.from.B << '\t' << m_input.to.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight << '\t' << m_answer.size())
  UNUSED(slave);

  for(int i = 0, c = m_answer.size(); i < c; ++i) {
    m_master->enqueueJob(new IntensityJob(m_master, IntensityInput(m_answer.at(i), m_input.orientation)));
  }
}

void FindRootsJob::saveTo(binary_oarchive& archive) const
{
  archive << m_input;
}

MPIJob* FindRootsJob::constructFrom(boost::archive::binary_iarchive& archive, MPIMaster* master)
{
  BisectInput input;
  archive >> input;
  return new FindRootsJob(master, input);
}

MPIJob::JobType FindRootsJob::type() const
{
  return FindRoots;
}

//BEGIN IntensityJob

IntensityJob::IntensityJob(MPIMaster* master, const IntensityInput& input)
: MPIJob(master)
, m_input(input)
{
  DEBUG_FUNC(<< m_input.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight)
}

IntensityJob::~IntensityJob()
{
  DEBUG_FUNC(<< m_input.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight)
}

void IntensityJob::start()
{
  int slave = m_master->runCommand(this, CMD_INTENSITY,
                                   TAG_INTENSITY_INPUT, m_input,
                                   TAG_INTENSITY_RESULT, m_answer);

  DEBUG_FUNC(<< slave << m_input.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight)
  UNUSED(slave);
}

void IntensityJob::handleResult(const int slave)
{
  DEBUG_FUNC(<< slave << m_input.B << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight << '\t' << m_answer)
  UNUSED(slave);
  m_master->intensityOutputFile() << m_input.B << '\t' << m_answer << '\t' << m_input.orientation.orientation.transpose() << '\t' << m_input.orientation.weight << endl;
}

void IntensityJob::saveTo(binary_oarchive& archive) const
{
  archive << m_input;
}

MPIJob* IntensityJob::constructFrom(boost::archive::binary_iarchive& archive, MPIMaster* master)
{
  IntensityInput input;
  archive >> input;
  return new IntensityJob(master, input);
}

MPIJob::JobType IntensityJob::type() const
{
  return Intensity;
}
