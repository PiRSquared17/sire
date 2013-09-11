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

#include "atommatcher.h"
#include "atomidx.h"
#include "atomname.h"
#include "moleculeinfodata.h"
#include "moleculeview.h"

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

/////////
///////// Implmentation of AtomMatcher
/////////

static const RegisterMetaType<AtomMatcher> r_atommatcher( MAGIC_ONLY,
                                                          "SireMol::AtomMatcher" );

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const AtomMatcher &matcher)
{
    writeHeader(ds, r_atommatcher, 1);
    
    ds << static_cast<const Property&>(matcher);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, AtomMatcher &matcher)
{
    VersionID v = readHeader(ds, r_atommatcher);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(matcher);
    }
    else
        throw version_error(v, "1", r_atommatcher, CODELOC);
        
    return ds;
}

/** Constructor */
AtomMatcher::AtomMatcher() : Property()
{}

/** Copy constructor */
AtomMatcher::AtomMatcher(const AtomMatcher &other) : Property(other)
{}

/** Destructor */
AtomMatcher::~AtomMatcher()
{}

/** Return a string representation of this matcher */
QString AtomMatcher::toString() const
{
    return QString("%s()").arg(this->what());
}
 
/** Return the match for the atoms in 'molview0' to the atoms in 'molview1' */   
QHash<AtomIdx,AtomIdx> AtomMatcher::match(const MoleculeView &molview0, 
                                          const MoleculeView &molview1) const
{
    return this->match( molview0.data().info(), molview1.data().info() );
}

/////////
///////// Implmentation of AtomIdxMatcher
/////////

static const RegisterMetaType<AtomIdxMatcher> r_idxmatcher;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const AtomIdxMatcher &idxmatcher)
{
    writeHeader(ds, r_idxmatcher, 1);
    ds << static_cast<const AtomMatcher&>(idxmatcher);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, AtomIdxMatcher &idxmatcher)
{
    VersionID v = readHeader(ds, r_idxmatcher);
    
    if (v == 1)
    {
        ds >> static_cast<AtomMatcher&>(idxmatcher);
    }
    else
        throw version_error(v, "1", r_idxmatcher, CODELOC);

    return ds;
}

/** Constructor */
AtomIdxMatcher::AtomIdxMatcher() : ConcreteProperty<AtomIdxMatcher,AtomMatcher>()
{}

/** Copy constructor */
AtomIdxMatcher::AtomIdxMatcher(const AtomIdxMatcher &other)
               : ConcreteProperty<AtomIdxMatcher,AtomMatcher>(other)
{}

/** Destructor */
AtomIdxMatcher::~AtomIdxMatcher()
{}

/** Copy assignment operator */
AtomIdxMatcher& AtomIdxMatcher::operator=(const AtomIdxMatcher &other)
{
    return *this;
}

/** Comparison operator */
bool AtomIdxMatcher::operator==(const AtomIdxMatcher &other) const
{
    return true;
}

/** Comparison operator */
bool AtomIdxMatcher::operator!=(const AtomIdxMatcher &other) const
{
    return false;
}

/** Return whether the application of this matcher will not change
    the order of the atoms (i.e. all the atoms in 'molinfo1' are
    in 'molinfo0', and they are all in the same order as they
    are in 'molinfo0') */
bool AtomIdxMatcher::unchangedAtomOrder(const MoleculeInfoData &molinfo0,
                                        const MoleculeInfoData &molinfo1) const
{
    return molinfo0.nAtoms() == molinfo1.nAtoms();
}

/** Match the atoms in 'molinfo1' to the atoms in 'molinfo0' - this
    returns the AtomIdxs of the atoms in 'molinfo1' that are in 
    'molinfo0', indexed by the AtomIdx of the atom in 'molinfo0'.
    
     This skips atoms in 'molinfo1' that are not in 'molinfo0'
*/
QHash<AtomIdx,AtomIdx> AtomIdxMatcher::match(const MoleculeInfoData &molinfo0,
                                             const MoleculeInfoData &molinfo1) const
{
    int nats = qMin( molinfo0.nAtoms(), molinfo1.nAtoms() );
    
    QHash<AtomIdx,AtomIdx> map;
    map.reserve(nats);
    
    for (int i=0; i<nats; ++i)
    {
        map.insert( AtomIdx(i), AtomIdx(i) );
    }
    
    return map;
}

const char* AtomIdxMatcher::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomIdxMatcher>() );
}

/////////
///////// Implmentation of AtomNameMatcher
/////////

static const RegisterMetaType<AtomNameMatcher> r_namematcher;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, 
                                       const AtomNameMatcher &namematcher)
{
    writeHeader(ds, r_namematcher, 1);
    ds << static_cast<const AtomMatcher&>(namematcher);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, AtomNameMatcher &namematcher)
{
    VersionID v = readHeader(ds, r_namematcher);
    
    if (v == 1)
    {
        ds >> static_cast<AtomMatcher&>(namematcher);
    }
    else
        throw version_error(v, "1", r_namematcher, CODELOC);

    return ds;
}

/** Constructor */
AtomNameMatcher::AtomNameMatcher() : ConcreteProperty<AtomNameMatcher,AtomMatcher>()
{}

/** Copy constructor */
AtomNameMatcher::AtomNameMatcher(const AtomNameMatcher &other)
               : ConcreteProperty<AtomNameMatcher,AtomMatcher>(other)
{}

/** Destructor */
AtomNameMatcher::~AtomNameMatcher()
{}

/** Copy assignment operator */
AtomNameMatcher& AtomNameMatcher::operator=(const AtomNameMatcher &other)
{
    return *this;
}

/** Comparison operator */
bool AtomNameMatcher::operator==(const AtomNameMatcher &other) const
{
    return true;
}

/** Comparison operator */
bool AtomNameMatcher::operator!=(const AtomNameMatcher &other) const
{
    return false;
}

/** Return whether the application of this matcher will not change
    the order of the atoms (i.e. all the atoms in 'molinfo1' are
    in 'molinfo0', and they are all in the same order as they
    are in 'molinfo0') */
bool AtomNameMatcher::unchangedAtomOrder(const MoleculeInfoData &molinfo0,
                                         const MoleculeInfoData &molinfo1) const
{
    if (molinfo0.nAtoms() != molinfo1.nAtoms())
        return false;
        
    int nats = molinfo1.nAtoms();
    
    for (AtomIdx i(0); i<nats; ++i)
    {
        try
        {
            if (molinfo0.atomIdx( molinfo1.name(i) ) != i )
                return false;
        }
        catch(...)
        {
            return false;
        }
    }
    
    return true;
}

/** Match the atoms in 'molinfo1' to the atoms in 'molinfo0' - this
    returns the AtomIdxs of the atoms in 'molinfo1' that are in 
    'molinfo0', indexed by the AtomIdx of the atom in 'molinfo0'.
    
     This skips atoms in 'molinfo1' that are not in 'molinfo0'
*/
QHash<AtomIdx,AtomIdx> AtomNameMatcher::match(const MoleculeInfoData &molinfo0,
                                              const MoleculeInfoData &molinfo1) const
{
    QHash<AtomIdx,AtomIdx> map;
    
    int nats = molinfo1.nAtoms();
    map.reserve(nats);
    
    for (AtomIdx i(0); i<nats; ++i)
    {
        try
        {
            map.insert( molinfo0.atomIdx( molinfo1.name(i) ), i );
        }
        catch(...){}
    }

    return map;
}

const char* AtomNameMatcher::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomNameMatcher>() );
}
