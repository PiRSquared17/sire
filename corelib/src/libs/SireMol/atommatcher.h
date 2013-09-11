/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREMOL_ATOMMATCHER_H
#define SIREMOL_ATOMMATCHER_H

#include <QHash>

#include "SireBase/property.h"

namespace SireMol
{
class AtomMatcher;
class AtomIdxMatcher;
class AtomNameMatcher;
}

QDataStream& operator<<(QDataStream&, const SireMol::AtomMatcher&);
QDataStream& operator>>(QDataStream&, SireMol::AtomMatcher&);

QDataStream& operator<<(QDataStream&, const SireMol::AtomIdxMatcher&);
QDataStream& operator>>(QDataStream&, SireMol::AtomIdxMatcher&);

QDataStream& operator<<(QDataStream&, const SireMol::AtomNameMatcher&);
QDataStream& operator>>(QDataStream&, SireMol::AtomNameMatcher&);

namespace SireMol
{

class AtomIdx;
class AtomName;
class MoleculeView;
class MoleculeInfoData;

/** Virtual base class of all of the functions used to match
    atoms in one molecule layout with atoms in another layout
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomMatcher : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const AtomMatcher&);
friend QDataStream& ::operator>>(QDataStream&, AtomMatcher&);

public:
    AtomMatcher();
    AtomMatcher(const AtomMatcher &other);
    
    virtual ~AtomMatcher();
    
    static const char* typeName()
    {
        return "SireMol::AtomMatcher";
    }
    
    virtual AtomMatcher* clone() const=0;
    
    virtual bool unchangedAtomOrder(const MoleculeInfoData &molinfo0,
                                    const MoleculeInfoData &molinfo1) const=0;
    
    virtual QHash<AtomIdx,AtomIdx> match(const MoleculeInfoData &molinfo0,
                                         const MoleculeInfoData &molinfo1) const=0;
    
    QHash<AtomIdx,AtomIdx> match(const MoleculeView &molview0, 
                                 const MoleculeView &molview1) const;
    
    virtual QString toString() const;
};

/** This is a simple atom matcher that matches the atoms based
    on their index in the molecule - e.g. it matches the first
    atom in molinfo0 to the first atom in molinfo1, the second
    atom in molinfo0 to the second atom in molinfo1, and the
    nth atom in molinfo0 to the nth atom in molinfo1
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomIdxMatcher
         : public SireBase::ConcreteProperty<AtomIdxMatcher,AtomMatcher>
{

friend QDataStream& ::operator<<(QDataStream&, const AtomIdxMatcher&);
friend QDataStream& ::operator>>(QDataStream&, AtomIdxMatcher&);

public:
    AtomIdxMatcher();
    AtomIdxMatcher(const AtomIdxMatcher&);
    
    ~AtomIdxMatcher();
    
    static const char* typeName();
    
    const char* what() const
    {
        return AtomIdxMatcher::typeName();
    }

    AtomIdxMatcher& operator=(const AtomIdxMatcher &other);
    
    bool operator==(const AtomIdxMatcher &other) const;
    bool operator!=(const AtomIdxMatcher &other) const;

    bool unchangedAtomOrder(const MoleculeInfoData &molinfo0,
                            const MoleculeInfoData &molinfo1) const;
    
    QHash<AtomIdx,AtomIdx> match(const MoleculeInfoData &molinfo0,
                                 const MoleculeInfoData &molinfo1) const;
};

/** This is a simple atom matcher that matches the atoms based
    on their names, so the atom called "CA1" in molinfo0 will
    be matched to the atom called "CA1" in molinfo1
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomNameMatcher
         : public SireBase::ConcreteProperty<AtomNameMatcher,AtomMatcher>
{

friend QDataStream& ::operator<<(QDataStream&, const AtomNameMatcher&);
friend QDataStream& ::operator>>(QDataStream&, AtomNameMatcher&);

public:
    AtomNameMatcher();
    AtomNameMatcher(const AtomNameMatcher&);
    
    ~AtomNameMatcher();
    
    static const char* typeName();
    
    const char* what() const
    {
        return AtomNameMatcher::typeName();
    }

    AtomNameMatcher& operator=(const AtomNameMatcher &other);
    
    bool operator==(const AtomNameMatcher &other) const;
    bool operator!=(const AtomNameMatcher &other) const;
    
    bool unchangedAtomOrder(const MoleculeInfoData &molinfo0,
                            const MoleculeInfoData &molinfo1) const;
    
    QHash<AtomIdx,AtomIdx> match(const MoleculeInfoData &molinfo0,
                                 const MoleculeInfoData &molinfo1) const;
};

typedef SireBase::PropPtr<AtomMatcher> AtomMatcherPtr;

}

Q_DECLARE_METATYPE( SireMol::AtomIdxMatcher )
Q_DECLARE_METATYPE( SireMol::AtomNameMatcher )

SIRE_EXPOSE_CLASS( SireMol::AtomMatcher )
SIRE_EXPOSE_CLASS( SireMol::AtomIdxMatcher )
SIRE_EXPOSE_CLASS( SireMol::AtomNameMatcher )

SIRE_EXPOSE_PROPERTY( SireMol::AtomMatcherPtr, SireMol::AtomMatcher )

#endif
