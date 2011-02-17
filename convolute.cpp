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
#include <vector>

#include <QtCore/QCoreApplication>
#include <QtCore/QStringList>
#include <QtCore/qtextstream.h>
#include <QtCore/QDir>
#include <QtCore/QMap>
#include <QtCore/QDebug>
#include <QtCore/QProcessEnvironment>

#include "types.h"

using namespace std;

QTextStream qout(stdout);
QTextStream qerr(stderr);

const int DEFAULT_STEPS = 100;
const fp DEFAULT_WIDTH = 0.001;
const fp WIDTHS_TO_ZERO = 4.0; // = exp(-N^2 /2), with N = 4 this is approx. 0.033%
const bool DEFAULT_DERIV = false;

void usage()
{
  qerr << "convolute DATA_DIR [steps=" << DEFAULT_STEPS << "] [width=" << DEFAULT_WIDTH << "] [deriv=" << DEFAULT_DERIV << "]" << endl;
}

class Gaussian {
public:
  Gaussian(double center, double intensity, double width)
  : m_center(center)
  , m_intensity(intensity)
  , m_width(width)
  {
  }

  bool applies(const double x) const
  {
    return x >= m_center - WIDTHS_TO_ZERO * m_width && x <= m_center + WIDTHS_TO_ZERO * m_width;
  }

  fp y(const double x) const
  {
    if (!applies(x)) {
      return 0;
    } else {
      return m_intensity * exp( - pow(x - m_center, 2) / (2.0 * pow(m_width, 2)) );
    }
  }

  fp yDeriv(const double x) const
  {
    if (!applies(x)) {
      return 0;
    } else {
      return - m_intensity * (x - m_center) / pow(m_width, 2) * exp( - pow(x - m_center, 2) / (2.0 * pow(m_width, 2)) );
    }
  }

private:
  const double m_center;
  const double m_intensity;
  const double m_width;
};

int main(int argc, char* argv[]) {
  QCoreApplication app(argc, argv);
  const QProcessEnvironment env = QProcessEnvironment::systemEnvironment();

  if (app.arguments().count() < 2) {
    usage();
    return 1;
  }

  QDir dir(app.arguments().at(1));
  if (!dir.exists() || !dir.isReadable()) {
    qerr << "could not open data dir for reading:" << dir.path() << endl << endl;
    usage();
    return 2;
  }

  int steps = DEFAULT_STEPS;
  QString stepsStr;
  if (app.arguments().count() >= 3) {
    stepsStr = app.arguments().at(2);
  } else if (env.contains("CONVOLUTE_STEPS")) {
    stepsStr = env.value("CONVOLUTE_STEPS");
  }
  if (!stepsStr.isEmpty()) {
    bool ok = true;
    steps = stepsStr.toInt(&ok);
    if (!ok) {
      qerr << "invalid steps argument" << endl << endl;
      usage();
      return 3;
    }
  }

  fp width = DEFAULT_WIDTH;
  QString widthStr;
  if (app.arguments().count() >= 4) {
    widthStr = app.arguments().at(3);
  } else if (env.contains("CONVOLUTE_WIDTH")) {
    widthStr = env.value("CONVOLUTE_WIDTH");
  }
  if (!widthStr.isEmpty()) {
    bool ok = true;
    width = widthStr.toFloat(&ok);
    if (!ok) {
      qerr << "invalid width argument" << endl << endl;
      usage();
      return 3;
    }
  }

  bool deriv = DEFAULT_DERIV;
  QString derivStr;
  if (app.arguments().count() >= 4) {
    derivStr = app.arguments().at(3);
  } else if (env.contains("CONVOLUTE_DERIV")) {
    derivStr = env.value("CONVOLUTE_DERIV");
  }
  if (!derivStr.isEmpty()) {
    bool ok = true;
    deriv = derivStr.toFloat(&ok);
    if (!ok) {
      qerr << "invalid deriv argument" << endl << endl;
      usage();
      return 3;
    }
  }

  // step 1: read data
  QMap<fp, Gaussian*> data;

  foreach(const QString& dataFile, dir.entryList(QDir::Files)) {
    QFile file(dir.absoluteFilePath(dataFile));
    if (!file.open(QIODevice::ReadOnly)) {
      qerr << "could not open data file for reading:" << dataFile << endl;
      return 3;
    }
    while(!file.atEnd()) {
      QStringList tuple = QString::fromLocal8Bit(file.readLine()).split(QLatin1Char('\t'));
      fp B = tuple.first().toDouble();
      fp I = tuple.last().toDouble();
      data[B] = new Gaussian(B, I, width);
    }
  }

  // step 2: sum each peak
  fp lastMax = -1;
  const fp stepSize = 2.0 * width * WIDTHS_TO_ZERO / steps;

  cout << (data.begin().key() - width * WIDTHS_TO_ZERO - stepSize) << '\t' << 0 << endl;

  foreach(const fp center, data.keys()) {
    if (lastMax != -1 && center - width * WIDTHS_TO_ZERO > lastMax + stepSize) {
      cout << (lastMax + stepSize) << '\t' << 0 << endl;
      cout << (center - width * WIDTHS_TO_ZERO - stepSize) << '\t' << 0 << endl;
    }
    for(fp x = max(lastMax, center - width * WIDTHS_TO_ZERO); x <= center + width * WIDTHS_TO_ZERO; x += stepSize) {
      fp y = 0;
      foreach(Gaussian* g, data) {
        if (deriv) {
          y += g->yDeriv(x);
        } else {
          y += g->y(x);
        }
      }
      cout << x << '\t' << y << endl;
    }
    lastMax = max(lastMax, center + width * WIDTHS_TO_ZERO);
  }

  cout << (lastMax + stepSize) << '\t' << 0 << endl;

  qDeleteAll(data.values());
  data.clear();

  return 0;
}
