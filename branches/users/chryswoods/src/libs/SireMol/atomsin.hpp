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

#ifndef SIREMOL_ATOMSIN_HPP
#define SIREMOL_ATOMSIN_HPP

#include "SireID/index.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

class MoleculeInfoData;

/** This helper class is used to provide the '.atoms()' functionality
    of the group ID classes. This allows the class to an atom, or
    range of atoms by index from the group that has been identified.
    
    @author Christopher Woods
*/
template<class GROUP>
class AtomsIn : public AtomID
{

public:
    AtomsIn();
    
    AtomsIn(const GROUP &id);
    AtomsIn(const GROUP &id, qint32 i);
    AtomsIn(const GROUP &id, qint32 i, qint32 j);
    
    AtomsIn(const AtomsIn<GROUP> &other);
    
    ~AtomsIn();
    
    static const char* typeName()
    {
        return QMetaType::typeName( qMetaTypeId< AtomsIn<GROUP> >() );
    }
    
    const char* what() const
    {
        return AtomsIn<GROUP>::typeName();
    }
    
    AtomsIn<GROUP>* clone() const
    {
        return new AtomsIn<GROUP>(*this);
    }
    
    AtomsIn<GROUP>& operator=(const AtomsIn<GROUP> &other);

    bool operator==(const AtomsIn<GROUP> &other) const;
    bool operator==(const SireID::ID &other) const;
    
    bool operator!=(const AtomsIn<GROUP> &other) const
    {
        return not this->operator==(other);
    }
    
    bool operator!=(const SireID::ID &other) const
    {
        return not this->operator==(other);
    }
    
    bool isNull() const
    {
        return groupid.isNull() and strt == 0 and end == -1;
    }
    
    uint hash() const
    {
        return groupid.hash() + strt + end;
    }
    
    QString toString() const;
    
    QList<AtomIdx> map(const MoleculeInfoData &molinfo) const;

private:
    /** The ID of the group that contains the atoms */
    typename GROUP::Identifier groupid;

    /** The indicies of the range of atoms that this specifies */
    SireID::Index strt,end;
};

/** Null constructor */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::AtomsIn()
               : AtomID(), strt(0), end(-1)
{}

/** Construct to get all of the atoms in the passed group */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::AtomsIn(const GROUP &id)
               : AtomID(), groupid(id), strt(0), end(-1)
{}

/** Construct for a specified atom */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::AtomsIn(const GROUP &id, qint32 i)
            : AtomID(), groupid(id), strt(i), end(i)
{}

/** Construct for a range of atoms */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::AtomsIn(const GROUP &id, qint32 i, qint32 j)
            : AtomID(), groupid(id), strt(i), end(j)
{}

/** Copy constructor */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::AtomsIn(const AtomsIn<GROUP> &other)
            : AtomID(other), groupid(other.groupid), 
              strt(other.strt), end(other.end)
{}

/** Destructor */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::~AtomsIn()
{}

/** Return a string representation of this ID */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
QString AtomsIn<GROUP>::toString() const
{
    if (strt == 0 and end == -1)
        return QString("(%1).atoms()").arg(groupid.toString());
    else if (strt == end)
        return QString("(%1).atom(%2)").arg(groupid.toString()).arg(strt);
    else
        return QString("(%1).atoms(%2,%3)")
                    .arg(groupid.toString()).arg(strt).arg(end);
}

/** Copy assignment operator */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>& AtomsIn<GROUP>::operator=(const AtomsIn<GROUP> &other)
{
    if (this != &other)
    {
        AtomID::operator=(other);
        groupid = other.groupid;
        strt = other.strt;
        end = other.end;
    }

    return *this;
}

/** Comparison operator */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool AtomsIn<GROUP>::operator==(const AtomsIn<GROUP> &other) const
{
    return strt == other.strt and end == other.end and
           groupid == other.groupid;
}

/** Comparison operator */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool AtomsIn<GROUP>::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare< AtomsIn<GROUP> >(*this, other);
}

}

SIRE_END_HEADER

#endif
