#include <iostream>

#include <boost/foreach.hpp>

#include "spinlib/orcaparser.h"
#include "spinlib/nucleus.h"
#include "spinlib/helpers.h"
#include "spinlib/experiment.h"

#include <boost/filesystem.hpp>

using namespace std;
using namespace Eigen;
using namespace boost::filesystem;

template<class T>
void printMatlab(const T& m)
{
  for(int i = 0; i < m.rows(); ++i) {
    if (i == 0) {
      cout << "[ ";
    } else {
      cout << "  ";
    }
    for(int j = 0; j < m.cols(); ++j) {
      cout << c_fp(m(i, j)).real();
      if (c_fp(m(i, j)).imag()) {
        cout << " + " << c_fp(m(i, j)).imag() << 'i';
      }
      cout << ' ';
    }
    if (i < m.rows() - 1) {
      cout << endl;
    } else {
      cout << "]";
    }
  }
}

int main(int argc, char** argv) {
  if (argc < 2) {
    cerr << "orca FILE missing" << endl;
    return 1;
  }

  int spinHalf = 0;
  int spinOne = 0;
  if (argc > 1) {
    spinHalf = atoi(argv[1]);
  }
  if (argc > 2) {
    spinOne = atoi(argv[2]);
  }
  path orcaFilePath(argv[1]);
  string orcaFile;
  if (exists(orcaFilePath)) {
    orcaFile = argv[1];
  }

  Experiment exp = getExperiment(orcaFile, spinHalf, spinOne);

  cout << "[Exp, Opt] = setupEasy();" << endl;
  cout << "Sys = struct();" << endl;
  cout << "Sys.g = "; printMatlab(exp.gTensor()); cout << ";" << endl;

  BOOST_FOREACH(const Nucleus& nuc, exp.nuclei) {
    stringstream id;
    id << nuc.isotope << *(nuc.name.end() - 1);
    cout << "A_full = ";
    printMatlab(nuc.A);
    cout << ";" << endl;
    cout << "[A, Apa] = euler(Sys.g, A_full, 1);" << endl;
    cout << "Sys = nucspinadd(Sys,";
    cout << '\'' << id.str() << '\'';
    cout << ", A, Apa);" << endl;
  }
  cout << "B_direction = ";
  printMatlab((exp.gTensorEigenVectors().col(2) / exp.gTensorEigenVectors().col(2).norm()).transpose());
  cout << ";" << endl;
  /*
  // Sys.Nucs = '1H,1H,...';
  stringstream nuclei;
  // Sys.A = [ [row1; row2; row3]; [...]; ... ];
  stringstream A;
  BOOST_FOREACH(const Nucleus& nuc, parser.nuclei()) {
    nuclei << nuc.isotope << *(nuc.name.end() - 1) << ',';
    // [5 0 0; 0 5 0; 0 0 5];

    A << '[';
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        A << nuc.A(i, j).real() << ' ';
      }
      A << "; ";
    }
    A << "];";
  }
  string nucs = nuclei.str();
  nucs.erase(nucs.size() - 1);
  cout << "Sys.Nucs = '" << nucs << "';" << endl;
  cout << "Sys.A = [ " << A.str() << " ];" << endl;
  */
  return 0;
}
