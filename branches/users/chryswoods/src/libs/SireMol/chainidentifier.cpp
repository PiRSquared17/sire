/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "chainidentifier.h"
#include "moleculeinfodata.h"

using namespace SireMol;
using namespace SireID;

/** Null constructor */
ChainIdentifier::ChainIdentifier() : ChainID()
{}

/** Construct from the passed ChainID */
ChainIdentifier::ChainIdentifier(const ChainID &chainid)
               : ChainID()
{
    if (chainid.isA<ChainIdentifier>())
        d = chainid.asA<ChainIdentifier>().d;
    else if (not chainid.isNull())
        d.reset( chainid.clone() );
}

/** Copy constructor */
ChainIdentifier::ChainIdentifier(const ChainIdentifier &other)
               : ChainID(other), d(other.d)
{}

/** Destructor */
ChainIdentifier::~ChainIdentifier()
{}

/** Is this selection null? */
bool ChainIdentifier::isNull() const
{
    return d.get() == 0;
}

/** Return a hash of this identifier */
uint ChainIdentifier::hash() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->hash();
}
            
/** Return a string representatio of this ID */
QString ChainIdentifier::toString() const
{
    if (d.get() == 0)
        return "null";
    else
        return d->toString();
}

/** Return the base type of this ID */
const ChainID& ChainIdentifier::base() const
{
    if (d.get() == 0)
        return *this;
    else
        return *d;
}

/** Copy assignment operator */
ChainIdentifier& ChainIdentifier::operator=(const ChainIdentifier &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
ChainIdentifier& ChainIdentifier::operator=(const ChainID &other)
{
    if (other.isA<ChainIdentifier>())
        d = other.asA<ChainIdentifier>().d;
    else if (other.isNull())
        d.reset();
    else
        d.reset(other.clone());
    
    return *this;
}

/** Comparison operator */
bool ChainIdentifier::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ChainIdentifier>(*this, other);
}

/** Comparison operator */
bool ChainIdentifier::operator==(const ChainIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() == other.d.get();
    else
        return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool ChainIdentifier::operator!=(const ChainIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() != other.d.get();
    else
        return d != other.d and *d != *(other.d);
}

/** Comparison operator */
bool ChainIdentifier::operator==(const ChainID &other) const
{
    if (d.get() == 0)
        return other.isNull();
    else if (other.isA<ChainIdentifier>())
        return this->operator==(other.asA<ChainIdentifier>());
    else
        return d->operator==(other);
}

/** Comparison operator */
bool ChainIdentifier::operator!=(const ChainID &other) const
{
    if (d.get() == 0)
        return not other.isNull();
    else if (other.isA<ChainIdentifier>())
        return this->operator!=(other.asA<ChainIdentifier>());
    else
        return d->operator!=(other);
}

/** Map this ID to the list of indicies of atoms that match this ID

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QList<ChainIdx> ChainIdentifier::map(const MoleculeInfoData &molinfo) const
{
    if (d.get() == 0)
        return molinfo.getChains();
    else
        return d->map(molinfo);
}

