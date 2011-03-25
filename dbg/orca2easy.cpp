#include <iostream>

#include <boost/foreach.hpp>

#include "spinlib/orcaparser.h"
#include "spinlib/nucleus.h"

using namespace std;
using namespace Eigen;

void printMatlab(const Matrix3c& m)
{
  for(int i = 0; i < 3; ++i) {
    if (i == 0) {
      cout << "[ ";
    } else {
      cout << "  ";
    }
    for(int j = 0; j < 3; ++j) {
      cout << m(i, j).real();
      if (m(i, j).imag()) {
        cout << " + " << m(i, j).imag() << 'i';
      }
      cout << ' ';
    }
    if (i < 2) {
      cout << endl;
    } else {
      cout << "]";
    }
  }
}

void printMatlab(const Vector3& v)
{
  cout << "[ ";
  for(int i = 0; i < 3; ++i) {
    cout << v(i);
    if (i + 1 < 3) {
      cout << ' ';
    }
  }
  cout << "]";
}

int main(int argc, char** argv) {
  if (argc < 2) {
    cerr << "orca FILE missing" << endl;
    return 1;
  }

  OrcaParser parser(argv[1]);

  cout << "Sys = struct();" << endl;
  cout << "Sys.g = "; printMatlab(parser.electronGMatrix()); cout << ";" << endl;
  const map<string, OrcaParser::AnglePrincipalPair>& euler = parser.eulerRotation();
  BOOST_FOREACH(const Nucleus& nuc, parser.nuclei()) {
    stringstream id;
    id << nuc.isotope << *(nuc.name.end() - 1);
    cout << "Sys = nucspinadd(Sys,";
    cout << '\'' << id.str() << '\'';
    cout << ", ";
    printMatlab(euler.at(nuc.name).first);
    cout << ", ";
    printMatlab(euler.at(nuc.name).second);
    cout << ");" << endl;
  }

  return 0;
}
