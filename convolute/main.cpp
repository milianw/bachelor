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
#include <fstream>

#include <QtCore/QCoreApplication>
#include <QtCore/QStringList>
#include <QtCore/qtextstream.h>
#include <QtCore/QDir>
#include <QtCore/QMap>
#include <QtCore/QDebug>
#include <QtCore/QProcessEnvironment>

#include "spinlib/fptype.h"

using namespace std;

QTextStream qout(stdout);
QTextStream qerr(stderr);

const int DEFAULT_STEPS = 100;
const fp DEFAULT_WIDTH = 0.001;
const fp WIDTHS_TO_ZERO = 4.0; // = exp(-N^2 /2), with N = 4 this is approx. 0.033%
const bool DEFAULT_DERIV = false;
const bool DEFAULT_NORMALIZE = true;

void usage()
{
  qerr << "convolute DATA_DIR [steps=" << DEFAULT_STEPS << "] [width=" << DEFAULT_WIDTH << "] [deriv=" << DEFAULT_DERIV << "] [norm=" << DEFAULT_NORMALIZE << "]" << endl;
}

class Gaussian {
public:
  Gaussian(double center, double intensity, double width, int orientation)
  : m_center(center)
  , m_intensity(intensity)
  , m_width(width)
  , m_orientation(orientation)
  {
  }

  bool applies(const double x) const
  {
    return x >= m_center - WIDTHS_TO_ZERO * m_width && x <= m_center + WIDTHS_TO_ZERO * m_width;
  }

  fp y(const double x) const
  {
    return m_intensity * exp( - pow(x - m_center, 2) / (2.0 * pow(m_width, 2)) );
  }

  fp yDeriv(const double x) const
  {
    return - m_intensity * (x - m_center) / pow(m_width, 2) * exp( - pow(x - m_center, 2) / (2.0 * pow(m_width, 2)) );
  }

  void normalizeIntensity(const double maxValue)
  {
    m_intensity /= maxValue;
  }

  double center() const
  {
    return m_center;
  }

  double height() const
  {
    return m_intensity;
  }

  void addIntensity(fp intensity)
  {
    m_intensity += intensity;
  }

  int orientation() const
  {
    return m_orientation;
  }

private:
  double m_center;
  double m_intensity;
  double m_width;

  int m_orientation;
};

struct Orientation {
  fp x;
  fp y;
  fp z;
  fp w;

  bool operator==(const Orientation& o) const
  {
    return o.x == x && o.y == y && o.z == z && o.w == w;
  }
};

int main(int argc, char* argv[]) {
  QCoreApplication app(argc, argv);
  const QProcessEnvironment env = QProcessEnvironment::systemEnvironment();

  if (app.arguments().count() < 2) {
    usage();
    return 1;
  }

  QFileInfo info(app.arguments().at(1));
  if (!info.exists() || !info.isReadable()) {
    qerr << "data not readable or does not exist:" << info.filePath();
    usage();
    return 2;
  }

  QStringList dataFiles;

  if (info.isDir()) {
    QDir dir(info.absoluteFilePath());
    foreach(const QString& file, dir.entryList(QDir::Files)) {
      dataFiles << dir.absoluteFilePath(file);
    }
  } else {
    dataFiles << info.absoluteFilePath();
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

  bool normalize = DEFAULT_NORMALIZE;
  QString normalizeStr;
  if (app.arguments().count() >= 5) {
    normalizeStr = app.arguments().at(4);
  } else if (env.contains("CONVOLUTE_NORMALIZE")) {
    normalizeStr = env.value("CONVOLUTE_NORMALIZE");
  }
  if (!normalizeStr.isEmpty()) {
    bool ok = true;
    normalize = normalizeStr.toFloat(&ok);
    if (!ok) {
      qerr << "invalid normalize argument" << endl << endl;
      usage();
      return 3;
    }
  }

  // step 1: read data
  QMultiMap<fp, Gaussian*> data;
  fp maxI = 0;

  QList<Orientation> orientations;

  foreach(const QString& dataFile, dataFiles) {
    ifstream file(dataFile.toLocal8Bit());
    if (!file.is_open()) {
      qerr << "could not open data file for reading:" << dataFile << endl;
      return 3;
    }
    while(!file.eof()) {
      fp B;
      fp I;
      Orientation o;
      file >> B >> I >> o.x >> o.y >> o.z >> o.w;
      if (file.bad()) {
        file.clear();
        continue;
      }
      I *= o.w;
      if (!I) {
        continue;
      }
      if (I > maxI) {
        maxI = I;
      }
      int oIdx = orientations.indexOf(o);
      if (oIdx == -1) {
        oIdx = orientations.size();
        orientations << o;
      }
      data.insert(B, new Gaussian(B, I, width, oIdx));
    }
  }

  if (normalize) {
    QMultiMap< fp, Gaussian* >::iterator it = data.begin();
    while(it != data.end()) {
      it.value()->normalizeIntensity(maxI);
      // get rid of data far below the peak threshold
      if (it.value()->height() < 1E-6) {
        delete it.value();
        it = data.erase(it);
      } else {
        ++it;
      }
    }
  }

  // output raw data marked by third row == 1
  foreach(Gaussian* g, data) {
      cout << g->center() << '\t' << g->height() << "\t1" << endl;
  }


  // step 2: sum each peak
  fp lastMax = -1;
  const fp stepSize = 2.0 * width * WIDTHS_TO_ZERO / steps;

  cout << (data.begin().key() - width * WIDTHS_TO_ZERO - 50 * stepSize) << '\t' << 0 << '\t' << 0 << endl;
  cout << (data.begin().key() - width * WIDTHS_TO_ZERO - stepSize) << '\t' << 0 << '\t' << 0 << endl;

  QVector<fp> yPerOrientation;
  yPerOrientation.resize(orientations.size());

  // normalize faktor for summation
  fp sumNormFaktor = 0;
  foreach(const fp center, data.uniqueKeys()) {
    yPerOrientation.fill(0);
    foreach(Gaussian* g, data) {
      if (!g->applies(center)) {
        continue;
      }
      fp& y_now = yPerOrientation[g->orientation()];
      fp y = deriv ? g->yDeriv(center) : g->y(center);
      y_now = qMax(y_now, y);
    }
    fp y_tot = 0;
    for(int i = 0; i < orientations.size(); ++i) {
      y_tot += yPerOrientation.at(i);
    }
    sumNormFaktor = qMax(y_tot, sumNormFaktor);
  }

  foreach(const fp center, data.uniqueKeys()) {
    if (lastMax != -1 && center - width * WIDTHS_TO_ZERO > lastMax + stepSize) {
      cout << (lastMax + stepSize) << '\t' << 0 << '\t' << 0 << endl;
      cout << (center - width * WIDTHS_TO_ZERO - stepSize) << '\t' << 0 << '\t' << 0 << endl;
    }

    for(fp x = max(lastMax, center - width * WIDTHS_TO_ZERO); x <= center + width * WIDTHS_TO_ZERO; x += stepSize) {
      yPerOrientation.fill(0);
      foreach(Gaussian* g, data) {
        if (!g->applies(x)) {
          continue;
        }
        fp& y_now = yPerOrientation[g->orientation()];
        fp y = deriv ? g->yDeriv(x) : g->y(x);
        y_now = qMax(y_now, y);
      }
      fp y_tot = 0;
      for(int i = 0; i < orientations.size(); ++i) {
        y_tot += yPerOrientation.at(i);
      }
      cout << x << '\t' << (y_tot / sumNormFaktor) << '\t' << 0 << endl;
    }

    lastMax = max(lastMax, center + width * WIDTHS_TO_ZERO);
  }

  cout << (lastMax + stepSize) << '\t' << 0 << '\t' << 0 << endl;
  cout << (lastMax + 50 * stepSize) << '\t' << 0 << '\t' << 0 << endl;

  qDeleteAll(data.values());
  data.clear();

  return 0;
}
