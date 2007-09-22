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

#ifndef SIREMOL_ATOMID_H
#define SIREMOL_ATOMID_H

#include "SireID/id.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

class MoleculeInfo;
class ResID;
class SegID;
class CGID;

class CGAtomIdx;

class SIREMOL_EXPORT AtomIDPart
{
public:
    AtomIDPart()
    {}
    
    virtual ~AtomIDPart()
    {}
    
    AtomIDCombo operator+(const AtomIDPart &other) const;
    
    virtual void update(AtomIDCombo &atomid) const=0;
};

/** This is the base class of all identifiers that are used 
    to identify an atom

    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomID : public SireID::ID, public AtomIDPart
{
public:
    AtomID() : SireID::ID()
    {}

    AtomID(const AtomID &other) : SireID::ID(other)
    {}

    ~AtomID()
    {}

    using AtomIDPart::operator+;

    virtual AtomID* clone() const=0;

    /** Map this ID back to the indicies of the matching atoms in the molecule, 
        using the passed MoleculeInfo to do the mapping */
    virtual QSet<AtomIdx> map(const MoleculeInfo &molinfo) const=0;
};

}

SIRE_END_HEADER

#endif
