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

#ifndef SIREMOL_ATOMATOMID_H
#define SIREMOL_ATOMATOMID_H

#include "atomidentifier.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

/** This class holds a pair of AtomIDs, which are used together
    to identify atom(s)
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomAtomID : public AtomID
{
public:
    AtomAtomID();
    
    AtomAtomID(const AtomID &id0, const AtomID &id1);
    
    AtomAtomID(const AtomAtomID &other);
    
    ~AtomAtomID();
    
    static const char* typeName()
    {
        return QMetaType::typeName( qMetaTypeId<AtomAtomID>() );
    }
    
    const char* what() const
    {
        return AtomAtomID::typeName();
    }

    AtomAtomID* clone() const
    {
        return new AtomAtomID(*this);
    }

    QString toString() const;

    bool operator==(const AtomAtomID &other) const
    {
        return (atomid0 == other.atomid0 and atomid1 == other.atomid1) or
               (atomid0 == other.atomid1 and atomid1 == other.atomid0);
    }

    bool operator!=(const AtomAtomID &other) const
    {
        return not this->operator==(other);
    }

    bool operator==(const SireID::ID &other) const
    {
        return SireID::ID::compare<AtomAtomID>(*this, other);
    }

    uint hash() const
    {
        return atomid0.hash() + atomid1.hash();
    }

    bool isNull() const
    {
        return atomid0.isNull() and atomid1.isNull();
    }

    QList<AtomIdx> map(const MolInfo &molinfo) const;

private:
    /** The pair of IDs that are combined */
    AtomIdentifier atomid0, atomid1;
};

}

Q_DECLARE_METATYPE(SireMol::AtomAtomID);

SIRE_EXPOSE_CLASS( SireMol::AtomAtomID )

SIRE_END_HEADER

#endif
