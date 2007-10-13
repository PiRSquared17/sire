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

#include "residentifier.h"
#include "moleculeinfodata.h"

using namespace SireMol;
using namespace SireID;

/** Null constructor */
ResIdentifier::ResIdentifier() : ResID()
{}

/** Construct from the passed ResID */
ResIdentifier::ResIdentifier(const ResID &resid)
              : ResID()
{
    if (resid.isA<ResIdentifier>())
        d = resid.asA<ResIdentifier>().d;
    else if (not resid.isNull())
        d.reset( resid.clone() );
}

/** Copy constructor */
ResIdentifier::ResIdentifier(const ResIdentifier &other)
              : ResID(other), d(other.d)
{}

/** Destructor */
ResIdentifier::~ResIdentifier()
{}

/** Is this selection null? */
bool ResIdentifier::isNull() const
{
    return d.get() == 0;
}

/** Return a hash of this identifier */
uint ResIdentifier::hash() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->hash();
}
            
/** Return a string representatio of this ID */
QString ResIdentifier::toString() const
{
    if (d.get() == 0)
        return "null";
    else
        return d->toString();
}

/** Return the base type of this ID */
const ResID& ResIdentifier::base() const
{
    if (d.get() == 0)
        return *this;
    else
        return *d;
}

/** Copy assignment operator */
ResIdentifier& ResIdentifier::operator=(const ResIdentifier &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
ResIdentifier& ResIdentifier::operator=(const ResID &other)
{
    if (other.isA<ResIdentifier>())
        d = other.asA<ResIdentifier>().d;
    else if (other.isNull())
        d.reset();
    else
        d.reset(other.clone());
    
    return *this;
}

/** Comparison operator */
bool ResIdentifier::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ResIdentifier>(*this, other);
}

/** Comparison operator */
bool ResIdentifier::operator==(const ResIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() == other.d.get();
    else
        return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool ResIdentifier::operator!=(const ResIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() != other.d.get();
    else
        return d != other.d and *d != *(other.d);
}

/** Comparison operator */
bool ResIdentifier::operator==(const ResID &other) const
{
    if (d.get() == 0)
        return other.isNull();
    else if (other.isA<ResIdentifier>())
        return this->operator==(other.asA<ResIdentifier>());
    else
        return d->operator==(other);
}

/** Comparison operator */
bool ResIdentifier::operator!=(const ResID &other) const
{
    if (d.get() == 0)
        return not other.isNull();
    else if (other.isA<ResIdentifier>())
        return this->operator!=(other.asA<ResIdentifier>());
    else
        return d->operator!=(other);
}

/** Map this ID to the list of indicies of atoms that match this ID

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QList<ResIdx> ResIdentifier::map(const MoleculeInfoData &molinfo) const
{
    if (d.get() == 0)
        return molinfo.getResidues();
    else
        return d->map(molinfo);
}
