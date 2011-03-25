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

#include "helpers.h"

#include "experiment.h"
#include "orcaparser.h"
#include "nucleus.h"
#include "spins.h"

#include <string>

#include <boost/foreach.hpp>

using namespace std;

Experiment getExperiment(const string& orcaInput, int protons, int nitrogens)
{
  if (!orcaInput.empty()) {
    OrcaParser parser(orcaInput);
    Experiment exp(parser.nuclei());
    exp.gTensor = parser.electronGMatrix();
    return exp;
  } else {
    return Experiment::generateDummy(protons, nitrogens);
  }
}

string formatSize(long unsigned int size) {
  const char prefix[5] = {' ', 'K', 'M', 'G', 'T'};
  int i = 0;
  for(i; i < 4; ++i) {
    if (size > 1024) {
      size /= 1024;
    } else {
      break;
    }
  }
  stringstream stream;
  stream << size << prefix[i] << 'B';
  return stream.str();
}

void printExperiment(ostream& out, const Experiment& exp)
{
  out << "nuclei:------" << endl << endl;
  BOOST_FOREACH(const Nucleus& nucleus, exp.nuclei) {
    out << "  " << nucleus.name << ", J = " << (nucleus.twoJ / 2) << ", g = " << (nucleus.g) << endl;
    out << "  A = " << nucleus.A << endl;
  }
  out << "------/nuclei" << endl;
  out << "gTensor:\n" << exp.gTensor << endl
      << "B direction:\n" << exp.staticBFieldDirection << endl;
  out << "mwFreq: " << exp.mwFreqGHz << "GHz" << endl;
}

string identifierForExperiment(const Experiment& exp)
{
  stringstream stream;
  Spins s = exp.spinSystem();
  stream << s.spinHalfs << ':' << s.spinOnes << ':' << exp.mwFreqGHz;
  return stream.str();
}

string guessPeakMemConsumption(const Experiment& exp)
{
  return formatSize(2 * (sizeof(complex<fp>) * exp.dimension * exp.dimension));
}

