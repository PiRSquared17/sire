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

#ifndef SIREMOL_CHAINNAME_H
#define SIREMOL_CHAINNAME_H

#include "SireID/name.h"

#include "chainid.h"
#include "molinfo.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class ChainName;
}

XMLStream& operator<<(XMLStream&, const SireMol::ChainName&);
XMLStream& operator>>(XMLStream&, SireMol::ChainName&);

namespace SireMol
{

/** This class holds the name of an atom. This can be used
    to identify an atom within a residue.
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT ChainName : public SireID::Name, public ChainID
{

friend XMLStream& ::operator<<(XMLStream&, const ChainName&);
friend XMLStream& ::operator>>(XMLStream&, ChainName&);

public:
    ChainName() : SireID::Name(), ChainID()
    {}
    
    explicit ChainName(const QString &name) : SireID::Name(name), ChainID()
    {}
    
    ChainName(const ChainName &other) : SireID::Name(other), ChainID(other)
    {}
    
    ~ChainName()
    {}
    
    static const char* typeName()
    {
        return "SireMol::ChainName";
    }
    
    const char* what() const
    {
        return ChainName::typeName();
    }
    
    ChainName* clone() const
    {
        return new ChainName(*this);
    }
    
    bool isNull() const
    {
        return SireID::Name::isNull();
    }
    
    uint hash() const
    {
        return qHash(_name);
    }
    
    QString toString() const
    {
        return QString("ChainName('%1')").arg(_name);
    }
    
    ChainName& operator=(const ChainName &other)
    {
        SireID::Name::operator=(other);
        ChainID::operator=(other);
        return *this;
    }
    
    bool operator==(const SireID::ID &other) const
    {
        return SireID::ID::compare<ChainName>(*this, other);
    }
    
    bool operator==(const ChainName &other) const
    {
        return _name == other._name;
    }
    
    bool operator!=(const ChainName &other) const
    {
        return _name != other._name;
    }

    QList<ChainIdx> map(const MolInfo &molinfo) const
    {
        return molinfo.map(*this);
    }
};

}

Q_DECLARE_METATYPE(SireMol::ChainName);

SIRE_EXPOSE_CLASS( SireMol::ChainName )

SIRE_END_HEADER

#endif

