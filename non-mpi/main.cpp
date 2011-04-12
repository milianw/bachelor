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

#include <QtCore/QCoreApplication>
#include <QtCore/QStringList>
#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtCore/QTextStream>

QTextStream& operator<<(QTextStream& out, long double d)
{
  return out << static_cast<double>(d);
}

#include <QtCore/QDebug>

QDebug& operator<<(QDebug& out, long double d)
{
  return out.operator<<(static_cast<double>(d));
}

#include <QtCore/QRegExp>
#include <QtCore/QDateTime>
#include <QtCore/QVector>
#include <QtCore/QStack>
#include <QtCore/QHash>

#include <omp.h>

#include "spinlib/resonancefield.h"
#include "spinlib/spinhamiltonian.h"
#include "spinlib/experiment.h"
#include "spinlib/helpers.h"

QTextStream qout(stdout);
QTextStream qerr(stderr);

using namespace std;

void usage() {
  cout << "hs-N OPTIONS" << endl
       << endl
       << " OPTIONS: (all are required, but -i and -p are mutually exclusive)" << endl
       << "  -n, --nprotons N  \tNumber of 1H nuclei" << endl
       << "  -N, --nitrogens N \tNumber of 14N nuclei" << endl
       << "  -i, --intensity B_MIN-B_MAX:STEPS:MW" << endl
       << "                    \tCalculate intensity for given MW over B_0 range" << endl
       << "         example: -i 0.2-0.4:1024:9.5 \t200 to 400mT, 1024 steps, 9.5GHz micro wave" << endl
       << "  -p, --peaks B     \tCalculate peaks for given B_0 (in T)" << endl
       << "         example: -p 0.2 \tpeaks at 0.2T" << endl
       << "  -o, --output DIR  \tOutput file path, each run creates a folder from timestamp" << endl
       << "                    \tand in there one file per thread gets created" << endl
       << "  -s, --system FILE \tOrca input file to use for system setup" << endl
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

  int nProtons = 0;
  int nNitrogens = 0;
  QString outputPath;
  QString orcaInput;

  // peaks
  fp B = 0;

  // intensity
  int steps = 0;
  fp B_min = 0;
  fp B_max = 0;
  fp mwFreq = 0;
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
          QRegExp pattern("(\\d+(?:\\.\\d+)?)-(\\d(?:\\.\\d+)?):(\\d+|auto):(\\d+(?:\\.\\d+)?)", Qt::CaseSensitive, QRegExp::RegExp2);
          if (pattern.exactMatch(intensityArg)) {
            B_min = pattern.cap(1).toDouble();
            B_max = pattern.cap(2).toDouble();
            steps = pattern.cap(3) == "auto" ? -1 : pattern.cap(3).toInt();
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
      } else if (arg == "--nitrogens" || arg == "-N") {
        if ((it+1) != end) {
          ++it;
          nNitrogens = it->toInt();
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
      } else if (arg == "--system" || arg == "-s") {
        if ((it+1) != end) {
          ++it;
          orcaInput = *it;
        } else {
          mode = Error;
          break;
        }
      }
      ++it;
    }
  }

  Experiment exp = getExperiment(orcaInput.toStdString(), nProtons, nNitrogens);

  ENSURE(!exp.nuclei.empty(), "nitrogens/protons/system")

  switch (mode) {
    case Error:
    {
      cerr << "CLI PARAM ERROR" << endl;
      usage();
      return 1;
    }
    case CalculateIntensity:
    {
      ENSURE(mwFreq > 0, "intensity")
      exp.mwFreqGHz = mwFreq;
      if (!B_max) {
        guessBRange(exp, B_min, B_max);
      }
      ENSURE(B_min >= 0, "intensity")
      ENSURE(B_max > B_min, "intensity")
      ENSURE(steps > 0 || steps == -1, "intensity")
      const double B_stepSize = (B_max - B_min) / steps;
      cout << "calculating intensity:" << endl
           << "B range:\t" << B_min << "T to " << B_max << "T" << endl
           << "steps:\t" << steps << endl;
      printExperiment(cout, exp);
      cerr << "max OMP threads:\t" << omp_get_max_threads() << endl;

      cerr << endl
           << "peak mem consumption at least:" << guessPeakMemConsumption(exp) << endl;

      QDir outputDir(outputPath);
      ENSURE(outputDir.exists(), "--output")
      QString base = QString("%1-%2:%3:%4").arg(B_min).arg(B_max).arg(steps).arg(QString::fromStdString(identifierForExperiment(exp)));
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

      if (steps != -1) { // dirty version that simply devides the max-min range into equal sized steps
        #pragma omp parallel
        {
          #pragma omp for
          for(int i = 0; i < steps; ++i) {
            const fp B = B_min + B_stepSize * i;
            *outputStreams.at(omp_get_thread_num()) << B << '\t' << SpinHamiltonian(B, exp).calculateIntensity() << endl;
          }
        }
      } else {
        // wip implementation of "S. Stoll, A. Schweiger / Chemical Physics Letters 380 (2003) 464 - 470"
        const vector<fp> resonanceField = ResonanceField(exp).calculate(B_min, B_max);
        #pragma omp parallel
        {
          #pragma omp for
          for(unsigned int i = 0; i < resonanceField.size(); ++i) {
            const fp B = resonanceField.at(i);
            *outputStreams.at(omp_get_thread_num()) << B << '\t' << SpinHamiltonian(B, exp).calculateIntensity() << endl;
          }
        }
      }
      qout << outputFileTpl.arg("") << endl;
      return 0;
    }
    case CalculatePeaks:
    {
      ENSURE(B > 0, "peaks")
      cout << "calculating peaks:" << endl
           << "B:\t" << B << endl;
      printExperiment(cout, exp);

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

