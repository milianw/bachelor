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

#include <boost/mpi.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

namespace mpi = boost::mpi;
namespace po = boost::program_options;
using namespace std;

const int SERVER_RANK = 0;

enum Commands {
  CMD_CLOSE,
  CMD_BISECT,
//   CMD_INTENSITY
};

enum Tags {
  TAG_CMD,
  TAG_BISECT_RANGE,
  TAG_TEST
};

int main(int argc, char* argv[]) {
  mpi::environment env(argc, argv);
  mpi::communicator world;

  if (world.rank() == SERVER_RANK) {
    // server

    // Declare the supported options.
    po::options_description desc("OPTIONS");
    desc.add_options()
      ("help", "show help message")
      ("from", po::value<double>()->default_value(0), "minimum B range in Tesla")
      ("to", po::value<double>()->default_value(1), "maximum B range in Tesla")
      ("mwFreq", po::value<double>(), "micro wave frequency in GHz")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    bool returnEarly = false;
    int retVal = 0;
    if (vm.count("help")) {
      cout << desc << endl;
      retVal = 0;
      returnEarly = true;
    }
    if (!vm.count("mwFreq")) {
      cerr << "missng mwFreq argument" << endl;
      retVal = 1;
      returnEarly = true;
    }

    if (returnEarly) {
      for (int i = 1; i < world.size(); ++i) {
        world.isend(i, TAG_CMD, CMD_CLOSE);
      }
      return retVal;
    }

    const double B_min = vm["from"].as<double>();
    const double B_max = vm["to"].as<double>();
    const double mwFreqGHz = vm["mwFreq"].as<double>();

    const int nodes = world.size() - 1; // one node is the server
    // dispatch commands to nodes
    vector<mpi::request> requests(nodes);
    vector<string> responses(nodes);

    for (int i = 1; i < world.size(); ++i) {
      world.isend(i, TAG_CMD, CMD_BISECT);
      requests[i - 1] = world.irecv(i, TAG_TEST, responses.at(i - 1));
    }
    // accumulate responses
    mpi::wait_all(requests.begin(), requests.end());
    cout << "responses are:" << endl;
    for (int i = 1; i < world.size(); ++i) {
      cout << i << "\t" << responses.at(i - 1) << endl;
      world.isend(i, TAG_CMD, CMD_CLOSE);
    }
  } else {
    // node
    while(true) {
      int cmd;
      world.recv(SERVER_RANK, TAG_CMD, cmd);
      if (cmd == CMD_CLOSE) {
        break;
      } else if (cmd == CMD_BISECT) {
        world.isend(SERVER_RANK, TAG_TEST, env.processor_name() + " - " + boost::lexical_cast<string>(world.rank()));
      }
    }
  }

  return 0;
}