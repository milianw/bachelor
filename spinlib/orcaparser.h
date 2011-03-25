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

#ifndef MW_BACHELOR_ORCAPARSER_H
#define MW_BACHELOR_ORCAPARSER_H

#include "eigentypes.h"

#include <vector>

class Nucleus;

/**
 * Parser for Orca data files
 */
class OrcaParser {
public:
  explicit OrcaParser(const std::string& file);
  ~OrcaParser();

  ///TODO: complex data allowed?
  Matrix3c electronGMatrix() const;

  std::vector<Nucleus> nuclei() const;

private:
  void parseFile(const std::string& file);

  ///TODO: complex data?
  Matrix3c m_gMatrix;
  std::vector<Nucleus> m_nuclei;
};

#endif // MW_BACHELOR_ORCAPARSER_H
