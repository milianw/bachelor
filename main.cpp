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

#include <QtCore/QCoreApplication>
#include <QtCore/QStringList>
#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtCore/QTextStream>
#include <QtCore/QRegExp>
#include <QtCore/QDateTime>
#include <QtCore/QVector>

#include <omp.h>

QTextStream qout(stdout);
QTextStream qerr(stderr);

void usage() {
  cout << "hs-N OPTIONS" << endl
       << endl
       /* TODO
       << " PARAMFILE:" << endl
       << "  a file that contains the parameters for the simulation" << endl
       << "  it should be structured as follows:" << endl
       << endl
       ///TODO: dynamic size of coupled nuclei? would require bitset alternative
//        << "  \tnatoms: n    # int, number of coupled nuclei" << endl
       << "  \tfield: x y z # 3x double vector, direction of static field" << endl
       << "  \tg:           # 3x3 double matrix, g tensor" << endl
       << "  \t  xx xy xz" << endl
       << "  \t  yx yy yz" << endl
       << "  \t  zx zy zz" << endl
       ///TODO: A tensor for H, N, ...
       << "  \tA:           # 3x3 double matrix, A tensor" << endl
       << "  \t  xx xy xz" << endl
       << "  \t  yx yy yz" << endl
       << "  \t  zx zy zz" << endl
       << endl
       */
       << " OPTIONS: (all are required, but -i and -p are mutually exclusive)" << endl
       << "  -n, --nprotons N  \tNumber of coupled nuclei" << endl
       << "  -i, --intensity B_MIN-B_MAX:STEPS:MW" << endl
       << "                    \tCalculate intensity for given MW over B_0 range" << endl
       << "         example: -i 0.2-0.4:1024:9.5 \t200 to 400mT, 1024 steps, 9.5GHz micro wave" << endl
       << "  -p, --peaks B     \tCalculate peaks for given B_0 (in T)" << endl
       << "         example: -p 0.2 \tpeaks at 0.2T" << endl
       << "  -o, --output DIR  \tOutput file path, each run creates a folder from timestamp" << endl
       << "                    \tand in there one file per thread gets created" << endl
       << "  -h, --help        \tDisplay help" << endl
       ;
}

#define ENSURE(cond, param) \
  if(!(cond)) { cerr << "ERROR: invalid parameter:" param << "\t" # cond << endl; return 1; }

int main(int argc, char* argv[])
{
  QCoreApplication app(argc, argv);

  enum Mode {
    Error,
    CalculateIntensity,
    CalculatePeaks,
  };

  Mode mode = Error;

  int nProtons = -1;
  QString outputPath;

  // peaks
  double B = 0;

  // intensity
  int steps = 0;
  double B_min = 0;
  double B_max = 0;
  double mwFreq = 0;
  QString intensityArg;

  { // argument parsing
    const QStringList args = app.arguments();
    QList<QString>::const_iterator it = args.constBegin() + 1;
    const QList<QString>::const_iterator end = args.constEnd();
    while (it != end) {
      const QString arg = *it;
      if (arg == "--intensity" || arg == "-i") {
        if ((it+1) != end) {
          ++it;
          intensityArg = *it;
          QRegExp pattern("(\\d+(?:\\.\\d+)?)-(\\d(?:\\.\\d+)?):(\\d+):(\\d+(?:\\.\\d+)?)", Qt::CaseSensitive, QRegExp::RegExp2);
          if (pattern.exactMatch(intensityArg)) {
            B_min = pattern.cap(1).toDouble();
            B_max = pattern.cap(2).toDouble();
            steps = pattern.cap(3).toInt();
            mwFreq = pattern.cap(4).toDouble();
          } else {
            mode = Error;
            break;
          }
          mode = CalculateIntensity;
        } else {
          mode = Error;
          break;
        }
      } else if (arg == "--peaks" || arg == "-p") {
        if ((it+1) != end) {
          ++it;
          B = it->toDouble();
          mode = CalculatePeaks;
        } else {
          mode = Error;
          break;
        }
      } else if (arg == "--nprotons" || arg == "-n") {
        if ((it+1) != end) {
          ++it;
          nProtons = it->toInt();
        } else {
          mode = Error;
          break;
        }
      } else if (arg == "--output" || arg == "-o") {
        if ((it+1) != end) {
          ++it;
          outputPath = *it;
        } else {
          mode = Error;
          break;
        }
      }
      ++it;
    }
  }

  ENSURE(nProtons >= 1 && nProtons <= 15, "nprotons")

  Experiment exp(nProtons);


  switch (mode) {
    case Error:
    {
      cerr << "CLI PARAM ERROR" << endl;
      usage();
      return 1;
    }
    case CalculateIntensity:
    {
      ENSURE(steps > 0, "intensity")
      ENSURE(B_min >= 0, "intensity")
      ENSURE(B_max > B_min, "intensity")
      ENSURE(mwFreq > 0, "intensity")
      const double B_stepSize = (B_max - B_min) / steps;
      cerr << "calculating intensity:" << endl
           << "mwFreq:\t" << mwFreq << "GHz" << endl
           << "B:\t" << B_min << "T to " << B_max << "T" << endl
           << "steps:\t" << steps << endl;
      cerr << "nProtons:\t" << exp.nProtons << endl
           << "aTensor:\n" << exp.aTensor << endl
           << "gTensor:\n" << exp.gTensor << endl
           << "B direction:\n" << exp.staticBFieldDirection << endl;
      cerr << "max OMP threads:\t" << omp_get_max_threads() << endl;

      QDir outputDir(outputPath);
      ENSURE(outputDir.exists(), "--output")
      QString base = intensityArg;
      ENSURE(outputDir.exists(base) || outputDir.mkdir(base), "--output");
      const QString outputFileTpl(outputDir.canonicalPath() + QDir::separator() + base + QDir::separator() + "%1");

      QVector<QFile*> outputFiles;
      QVector<QTextStream*> outputStreams;
      for(int i = 0; i < omp_get_max_threads(); ++i) {
        QFile* file = new QFile(outputFileTpl.arg(i));
        ENSURE(file->open(QIODevice::WriteOnly), "--output");
        outputFiles << file;
        QTextStream* out = new QTextStream(file);
        outputStreams << out;
      }

      #pragma omp parallel
      {
        #pragma omp for
        for(int i = 0; i < steps; ++i) {
          SpinHamiltonian(B_min + B_stepSize * i, exp).calculateIntensity(mwFreq, outputStreams.at(omp_get_thread_num()));
        }
      }

      qout << outputFileTpl.arg("") << endl;
      return 0;
    }
    case CalculatePeaks:
    {
      ENSURE(B > 0, "peaks")
      cerr << "calculating peaks:" << endl
           << "B:\t" << B << endl;
      cerr << "nProtons:\t" << exp.nProtons << endl
           << "aTensor:\n" << exp.aTensor << endl
           << "gTensor:\n" << exp.gTensor << endl
           << "B direction:\n" << exp.staticBFieldDirection << endl;

      SpinHamiltonian(B, exp).calculateTransitions();
      return 0;
    }
  }

  //Parsing input file
  //ifstream inputfile;
  //inputfile.open(argv[1]);
  // if(!inputfile)
  //   {
  //     cerr << "This program expects a single input file\n";
  //     cerr << "Check the comments in the source code for details\n"
  //     return(1);
  //   }

  return 0;
}

