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

#include "spinhamiltonian.h"

int main(int argc, char* argv[])
{
  //cout << 2.023 * Bohrm << endl;
  //cout << NUC_MAGNETON << endl;
  //cout << g_1H * NUC_MAGNETON << endl;
  //cout << GAMMA_1H * hbar << endl;

  /////////////////////////////////////////////////////////////////////////
  // The input file should be structured as follows:                 //
  // natoms n                     (number of coupled nuclei)         //
  // field x y z                  (The static field direction)       //
  // g                            (g tensor)                         //
  //   xx xy xz                              //
  //   yx yy yz                              //
  //   zx zy zz                              //
  // N                            (A tensor)                 //
  //   xx xy xz                              //
  //   yx yy yz                              //
  //   zx zy zz                              //
  // H                            (A tensor)                 //
  //   .........                             //
  //   .........                             //
  //   .........                             //
  // And so on...                                //
  /////////////////////////////////////////////////////////////////////////
  

  //Parsing input file
  //ifstream inputfile;
  //inputfile.open(argv[1]);
  // if(!inputfile)
  //   {
  //     cerr << "This program expects a single input file\n";
  //     cerr << "Check the comments in the source code for details\n"
  //     return(1);
  //   }

  /*
  {
//       const double B = 0.362562;
//       const double B = 0.3;
      const double B = 0.3417757;
      SpinHamiltonian h(B);
      h.calculateTransitions();
      return 0;
  }
  */

  const int steps = 1024;
  const double B_min = 0.280;
  const double B_max = 0.400;
  const double mwFreq = 9.5; // in GHz
  const double B_stepSize = (B_max - B_min) / steps;
  double B = B_min;
  for(int i = 0; i < steps; ++i) {
    SpinHamiltonian(B).calculateIntensity(mwFreq);
    B += B_stepSize;
  }

  return 0;
}

