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

#ifndef MW_BACHELOR_MPITYPES_H
#define MW_BACHELOR_MPITYPES_H

#include <boost/mpi.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/serialization/collection_size_type.hpp>
#include <boost/serialization/nvp.hpp>

#include "spinlib/eigentypes.h"
#include "spinlib/resonancefield.h"

namespace mpi = boost::mpi;

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive& ar, VectorX& vec, unsigned int version)
{
  int size;
  if ( !Archive::is_loading::value ) {
      size = vec.size();
  }
  ar & size;
  if ( Archive::is_loading::value ) {
      vec.resize(size);
  }
  ar & make_array(vec.data(), size);
}

template<class Archive>
void serialize(Archive& ar, BisectNode& node, unsigned int version)
{
  ar & node.B;
  ar & node.E;
  ar & node.E_deriv;
}

template<class Archive>
void serialize(Archive& ar, BisectAnswer& answer, const unsigned int version)
{
  ar & answer.status;
  ar & answer.from;
  ar & answer.to;
  ar & answer.mid;
}

} // namespace serialization
} // namespace boost

class BisectInput {
public:
  BisectInput()
  { }

  BisectInput(const BisectNode& _from, const BisectNode& _to)
  : from(_from), to(_to)
  { }

  BisectNode from;
  BisectNode to;

private:
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar & from;
    ar & to;
  }
};


class IntensityAnswer {
public:
  IntensityAnswer()
  { }

  IntensityAnswer(const fp _B, const fp _intensity)
  : B(_B)
  , intensity(_intensity)
  { }

  fp B;
  fp intensity;

private:
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar & B;
    ar & intensity;
  }
};

#endif // MW_BACHELOR_MPITYPES_H
