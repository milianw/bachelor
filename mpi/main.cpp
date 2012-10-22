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
#include <vector>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "mpi_master.h"
#include "mpi_slave.h"
#include "mpi_iface.h"

#include "spinlib/experiment.h"
#include "spinlib/orcaparser.h"
#include "spinlib/helpers.h"

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main(int argc, char* argv[]) {
  mpi::environment env(argc, argv);
  mpi::communicator world;

  // Declare the supported options.
  po::options_description desc("OPTIONS");
  desc.add_options()
    ("help,h", "show help message")
    ("protons,p", po::value<int>()->default_value(0), "number of protons in system")
    ("nitrogens,n", po::value<int>()->default_value(0), "number of nitrogens in system")
    ("system,s", po::value<string>()->default_value(string()), "ORCA data file describing the system")
    ("cutoffcount,c", po::value<int>()->default_value(0), "consider only the N highest HFC couplings")
    ("from,f", po::value<fp>()->default_value(0), "minimum B range in Tesla")
    ("to,t", po::value<fp>()->default_value(0), "maximum B range in Tesla")
    ("mwFreq,m", po::value<fp>()->required(), "micro wave frequency in GHz")
    ("outputDir,o", po::value<string>()->default_value(string()), "path for writing intensity data to, defaults to current working directory")
    ("steps,i", po::value<int>()->default_value(-1), "number of B-steps, the default is auto-adapted segmentation")
    ("orientations,O", po::value<string>()->default_value(string()), "file containing orientation vectors, default is a single [0 0 1] orientation.")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);

  if (vm.count("help")) {
    cout << desc << endl << endl
         << "NOTE: -s is mutually exclusive with -p -n, but either one must be given" << endl
    ;
    return 0;
  }

  try {
    po::notify(vm);
  } catch(boost::program_options::required_option e) {
    if (world.rank() == MASTER_RANK) {
      cerr << e.what() << endl;
    }
    return 1;
  }

  Experiment exp = getExperiment(vm["system"].as<string>(), vm["cutoffcount"].as<int>(), vm["protons"].as<int>(), vm["nitrogens"].as<int>());
  if (exp.nuclei.empty()) {
    cerr << "either protons, nitrogens or system must be set" << endl;
    return 1;
  }
  exp.mwFreqGHz = vm["mwFreq"].as<fp>();

  if (world.rank() == MASTER_RANK) {
    fp from = vm["from"].as<fp>();
    fp to = vm["to"].as<fp>();
    if (!to) {
      guessBRange(exp, from, to);
    }

    const int steps = vm["steps"].as<int>();
    cout << "calculating intensity in range B = \t" << from << "T to " << to << "T" << endl;
    printExperiment(cout, exp);
    cout << "worker slaves:" << (world.size() - 1) << endl;
    cout << endl
         << "peak mem consumption (on slaves) at least:" << guessPeakMemConsumption(exp) << endl;

    // read orientation file
    string orientationsFile = vm["orientations"].as<string>();
    vector<Orientation> orientations;
    if (!orientationsFile.empty()) {
      ifstream stream(orientationsFile.c_str());
      if (!stream.is_open()) {
        cerr << "could not read orientations file: " << orientationsFile << endl;
        world.abort(4);
        return 4;
      }
      while(!stream.eof()) {
        Orientation orientation;
        stream >> orientation.orientation(0) >> orientation.orientation(1) >> orientation.orientation(2) >> orientation.weight;
        if (stream.good()) {
          orientations.push_back(orientation);
        }
      }
    }
    if (orientations.empty()) {
      orientations.push_back(Orientation((Vector3() << 0, 0, 1).finished(), 1));
    }
    cout << "number of orientations: " << orientations.size() << endl;

    // find output dir
    string outputDir = vm["outputDir"].as<string>();
    if (outputDir.empty()) {
      // fallback to current work dir
      outputDir = fs::current_path().string();
    }

    // append info about the experiment
    stringstream stream;
    stream << outputDir << "/" << from << '-' << to << ':' << steps << ':' << identifierForExperiment(exp) << ":mpi:" << orientations.size();
    outputDir = stream.str();

    // make sure the path exists
    fs::path oPath(outputDir);
    outputDir = fs::system_complete(oPath).string();
    if (!fs::exists(oPath) && !fs::create_directories(oPath)) {
      cerr << "could not create output path: " << oPath << endl;
      world.abort(2);
      return 2;
    } else if (!fs::is_directory(oPath)) {
      cerr << "file exists under name of output path: " << oPath << endl;
      world.abort(3);
      return 3;
    }

    MPIMaster master(world, exp, outputDir);
    master.calculateIntensity(from, to, steps, orientations);
  } else {
    MPISlave slave(world, exp);
    slave.work();
  }

  return 0;
}