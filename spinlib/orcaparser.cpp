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

#include "orcaparser.h"

#include "constants.h"
#include "nucleus.h"

#include <fstream>
#include <iostream>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace Eigen;

OrcaParser::OrcaParser(const string& file)
{
  m_gMatrix = Matrix3::Identity() * Constants::g_E;
  parseFile(file);
}

OrcaParser::~OrcaParser()
{

}

Matrix3 OrcaParser::electronGMatrix() const
{
  return m_gMatrix;
}

vector<Nucleus> OrcaParser::nuclei() const
{
  return m_nuclei;
}

map<string, OrcaParser::AnglePrincipalPair> OrcaParser::eulerRotation() const
{
  return m_euler;
}

void OrcaParser::parseFile(const string& file)
{
  ifstream stream(file.data(), ios::in);
  if (!stream.is_open()) {
    throw "could not open input file";
  }

  string line;
  while(!stream.eof()) {
    getline(stream, line);
    if (boost::starts_with(line,"ELECTRONIC G-MATRIX")) {
      getline(stream, line); // ----...
      getline(stream, line); //
      getline(stream, line); // The g-matrix:
      for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
          stream >> m_gMatrix(i, j);
        }
      }
    } else if (boost::starts_with(line, " Nucleus")) {
      static boost::regex pattern("\\s*Nucleus\\s*([^: ]+)\\s*:"
                                  "\\s*A:ISTP=\\s*(\\d+)"
                                  "\\s+I=\\s*([0-9]+(?:\\.[0-9]+)?)"
                                  "\\s+P=\\s*([0-9]+(?:\\.[0-9]+)?)\\s+(.+)$");
      boost::smatch matches;
      if (!boost::regex_match(line, matches, pattern)) {
        cerr << "could not parse line:" << line << endl;
        throw "parser error";
      }

      string nucleus(matches.str(1));
      int isotope = boost::lexical_cast<int>(matches.str(2));
      fp spin = boost::lexical_cast<fp>(matches.str(3));
      fp g;
      switch(*(nucleus.end() - 1)) {
        case 'N':
          g = Constants::g_14N;
          break;
        case 'H':
          g = Constants::g_1H;
          break;
        default:
          cerr << "unhandled nucleus type, could not instantiate g-factor in line:" << line << endl;
          throw "parser error";
      }

      // Q = ...
      getline(stream, line);
      // ---...
      getline(stream, line);

      Matrix3 A = Matrix3::Zero();
      Vector3 Q = Vector3::Zero();
      bool haveA = false;
      bool haveQ = false;
      getline(stream, line);
      while(!stream.eof() && !boost::starts_with(line, " ---") && (!haveQ || !haveA)) {
        if(!haveA && boost::starts_with(line, " Raw HFC matrix (all values in MHz):")) {
          for(int i = 0; i < 3; ++i) {
            for(int j = 0; j < 3; ++j) {
              stream >> A(i, j);
            }
          }
          haveA = true;
        } else if (!haveQ && boost::starts_with(line, " Quadrupole tensor eigenvalues")) {
          string devNull;
          // e**2qQ = 
          stream >> devNull; stream >> devNull;
          fp alpha;
          stream >> alpha;
          // MHz
          stream >> devNull;
          // e**2qQ/(4I*(2I-1))=
          stream >> devNull;
          fp beta;
          stream >> beta;
          // MHz
          stream >> devNull;
          // eta = ;
          stream >> devNull; stream >> devNull;
          fp eta;
          stream >> eta;
          Q = (Vector3() << -1.0 + eta, -1.0 - eta, 2).finished() * beta;
          haveQ = true;
        }
        getline(stream, line);
      }
      m_nuclei.push_back(Nucleus(nucleus, spin * 2, isotope, A, g, Q));
    } else if (boost::contains(line, "Euler rotation of hyperfine tensor to g-tensor")) {
      getline(stream, line); // ---...
      getline(stream, line); //
      getline(stream, line); // ---...
      getline(stream, line); // Atom | ...
      getline(stream, line); //      | [degrees] ...
      getline(stream, line); // ---...
      while (!stream.eof()) {
        getline(stream, line);
        if (boost::starts_with(line, "---")) {
          break;
        }
        stringstream ls(line);
        string nucleus;
        ls >> nucleus;
        Vector3 angles;
        ls >> angles(0); ls >> angles(1); ls >> angles(2);
        Vector3 principals;
        ls >> principals(0); ls >> principals(1); ls >> principals(2);
        m_euler[nucleus] = AnglePrincipalPair(angles, principals);
      }
    }
  }
}
