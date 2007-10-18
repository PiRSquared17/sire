/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#ifndef SIREDB_ATOMPARAMETER_HPP
#define SIREDB_ATOMPARAMETER_HPP

#include "SireMol/atominfo.h"

#include "SireStream/datastream.h"

SIRE_BEGIN_HEADER

namespace SireDB
{
template<class Param>
class AtomParameter;
}

template<class Param>
QDataStream& operator<<(QDataStream&, const SireDB::AtomParameter<Param>&);
template<class Param>
QDataStream& operator>>(QDataStream&, SireDB::AtomParameter<Param>&);

namespace SireDB
{

using SireMol::AtomInfo;

/** This is a convienience class that provides the parameter for
    a specific atom. This class combines an AtomInfo, which provides
    all of the metadata about an atom, with the Parameter class itself.

    Sire separates the metadata for an Atom from its coordinates.
    SireMol::AtomInfo contains the metadata, while SireMaths::Vector
    contains the coordinates. Making this separation allows other
    data types to be combined with the meta-data. This class is
    used to combine forcefield parameter data types with the
    metadata, e.g. AtomCharge is an alias for AtomParameter<ChargeParameter>,
    which represents the combination of the atom's metadata with
    its charge. Similarly, AtomLJ is an alias for
    AtomParameter<LJParameter>, which represents the combination
    of the atom's metadata with its LJ parameter.

    \param Param  The type of parameter stored in this object

    @author Christopher Woods
*/
template<class Param>
class SIREDB_EXPORT AtomParameter : public AtomInfo, public Param
{

friend QDataStream& ::operator<<<>(QDataStream&, const AtomParameter<Param>&);
friend QDataStream& ::operator>><>(QDataStream&, AtomParameter<Param>&);

public:
    AtomParameter();
    AtomParameter(const AtomInfo &atominfo);
    AtomParameter(const AtomInfo &atominfo, const Param &param);

    AtomParameter(const AtomParameter<Param> &other);

    ~AtomParameter();

    AtomParameter<Param>& operator=(const AtomParameter &other);
    AtomParameter<Param>& operator=(const AtomInfo &other);
    AtomParameter<Param>& operator=(const Param &other);

    bool operator==(const AtomParameter<Param> &other) const;
    bool operator!=(const AtomParameter<Param> &other) const;
};

#ifndef SKIP_TEMPLATE_DEFINITIONS

/** Null constructor */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
AtomParameter<Param>::AtomParameter() : AtomInfo(), Param()
{}

/** Construct a default-constructed parameter for the specified
    atom */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
AtomParameter<Param>::AtomParameter(const AtomInfo &atominfo)
                     : AtomInfo(atominfo), Param()
{}

/** Construct an AtomParameter for the passed atom, setting it equal
    to the passed parameter */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
AtomParameter<Param>::AtomParameter(const AtomInfo &atominfo,
                                    const Param &param)
                     : AtomInfo(atominfo), Param(param)
{}

/** Copy constructor */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
AtomParameter<Param>::AtomParameter(const AtomParameter<Param> &other)
                     : AtomInfo(other), Param(other)
{}

/** Destructor */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
AtomParameter<Param>::~AtomParameter()
{}

/** Assignment operator */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
AtomParameter<Param>& AtomParameter<Param>::operator=(const AtomParameter &other)
{
    AtomInfo::operator=(other);
    Param::operator=(other);
    return *this;
}

/** Set the AtomInfo part of this object equal to 'other' */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
AtomParameter<Param>& AtomParameter<Param>::operator=(const AtomInfo &other)
{
    AtomInfo::operator=(other);
    return *this;
}

/** Set the parameter part of this object equal to 'other' */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
AtomParameter<Param>& AtomParameter<Param>::operator=(const Param &other)
{
    Param::operator=(other);
    return *this;
}

/** Comparison operator */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
bool AtomParameter<Param>::operator==(const AtomParameter<Param> &other) const
{
    return Param::operator==(other) and AtomInfo::operator==(other);
}

/** Comparison operator */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
bool AtomParameter<Param>::operator!=(const AtomParameter<Param> &other) const
{
    return Param::operator!=(other) or AtomInfo::operator!=(other);
}

static const MagicID atomparameter_magic = getMagic("SireDB::AtomParameter");

#endif // SKIP_TEMPLATE_DEFINITIONS

}

#ifndef SKIP_TEMPLATE_DEFINITIONS

/** Serialise to a binary datastream */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
QDataStream &operator<<(QDataStream &ds, const SireDB::AtomParameter<Param> &param)
{
    SireStream::writeHeader(ds, SireDB::atomparameter_magic, 1)
                        << static_cast<const SireMol::AtomInfo&>(param)
                        << static_cast<const Param&>(param);

    return ds;
}

/** Deserialise from a binary datastream */
template<class Param>
SIRE_OUTOFLINE_TEMPLATE
QDataStream &operator>>(QDataStream &ds, SireDB::AtomParameter<Param> &param)
{
    SireStream::VersionID v = SireStream::readHeader(ds, SireDB::atomparameter_magic,
                                                     "SireDB::AtomParameter");
    if (v == 1)
    {
        ds >> static_cast<SireMol::AtomInfo&>(param)
           >> static_cast<Param&>(param);
    }
    else
        throw SireStream::version_error(v, "1", "SireDB::AtomParameter", CODELOC);

    return ds;
}

#endif // SKIP_TEMPLATE_DEFINITIONS

SIRE_END_HEADER

#endif
