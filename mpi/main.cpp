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

#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "mpi_master.h"
#include "mpi_slave.h"
#include "mpi_iface.h"

#include "spinlib/experiment.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char* argv[]) {
  mpi::environment env(argc, argv);
  mpi::communicator world;

  // Declare the supported options.
  po::options_description desc("OPTIONS");
  desc.add_options()
    ("help,h", "show help message")
    ("protons,p", po::value<int>()->default_value(0), "number of protons in system")
    ("nitrogens,n", po::value<int>()->default_value(0), "number of nitrogens in system")
    ("from,f", po::value<fp>()->default_value(0), "minimum B range in Tesla")
    ("to,t", po::value<fp>()->default_value(1), "maximum B range in Tesla")
    ("mwFreq,m", po::value<fp>()->required(), "micro wave frequency in GHz")
    ("outputDir,o", po::value<string>()->required(), "path for writing intensity data to")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  try {
    po::notify(vm);
  } catch(boost::program_options::required_option e) {
    if (world.rank() == MASTER_RANK) {
      cerr << e.what() << endl;
    }
    return 1;
  }

  if (vm.count("help")) {
    cout << desc << endl;
    return 0;
  }

  Experiment exp(vm["protons"].as<int>(), vm["nitrogens"].as<int>());
  if (!exp.nNitrogens && !exp.nProtons) {
    cerr << "either protons or nitrogens must be set" << endl;
    return 1;
  }
  exp.mwFreqGHz = vm["mwFreq"].as<fp>();

  if (world.rank() == MASTER_RANK) {
    cout << "protons:" << exp.nProtons << endl
         << "mwFreq:" << exp.mwFreqGHz << " GHz" << endl;
    // master
    MPIMaster master(world, exp, vm["outputDir"].as<string>());

    master.startBisect(vm["from"].as<fp>(), vm["to"].as<fp>());
    cout << vm["outputDir"].as<string>();
  } else {
    MPISlave slave(world, exp);
    slave.work();
  }

  return 0;
}