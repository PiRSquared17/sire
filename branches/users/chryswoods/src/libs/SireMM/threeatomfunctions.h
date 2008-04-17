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

#ifndef SIREMM_THREEATOMFUNCTIONS_H
#define SIREMM_THREEATOMFUNCTIONS_H

#include "atomfunctions.h"

#include "SireMol/cgatomidx.h"
#include "SireMol/atomidx.h"
#include "SireMol/atomid.h"

#include "SireMol/angleid.h"

namespace SireMM
{
class ThreeAtomFunction;
class ThreeAtomFunctions;
}

QDataStream& operator<<(QDataStream&, const SireMM::ThreeAtomFunction&);
QDataStream& operator>>(QDataStream&, SireMM::ThreeAtomFunction&);

QDataStream& operator<<(QDataStream&, const SireMM::ThreeAtomFunctions&);
QDataStream& operator>>(QDataStream&, SireMM::ThreeAtomFunctions&);

namespace SireMM
{

using SireMol::CGAtomIdx;
using SireMol::AtomIdx;
using SireMol::AtomID;
using SireMol::AngleID;

/** This class holds a function that acts using the 
    coordinate information of just three atoms */
class SIREMM_EXPORT ThreeAtomFunction : public AtomFunction
{

friend QDataStream& ::operator<<(QDataStream&, const ThreeAtomFunction&);
friend QDataStream& ::operator>>(QDataStream&, ThreeAtomFunction&);

public:
    ThreeAtomFunction();
    ThreeAtomFunction(const CGAtomIdx &atom0, const CGAtomIdx &atom1,
                      const CGAtomIdx &atom2,
                      const SireCAS::Expression &function);
                  
    ThreeAtomFunction(const ThreeAtomFunction &other);
    
    ~ThreeAtomFunction();

    ThreeAtomFunction& operator=(const ThreeAtomFunction &other);
    
    bool operator==(const ThreeAtomFunction &other) const;
    bool operator!=(const ThreeAtomFunction &other) const;
    
    const CGAtomIdx& atom0() const;
    const CGAtomIdx& atom1() const;
    const CGAtomIdx& atom2() const;

private:
    /** The indicies of the three atoms */
    CGAtomIdx atm0, atm1, atm2;
};

namespace detail
{

class IDTriple
{
public:
    IDTriple(quint32 atom0=0, quint32 atom1=0, quint32 atom2=0);
    IDTriple(const IDTriple &other);
    
    ~IDTriple();
    
    IDTriple& operator=(const IDTriple &other);
    
    bool operator==(const IDTriple &other) const;
    bool operator!=(const IDTriple &other) const;
    
    quint32 atom0;
    quint32 atom1;
    quint32 atom2;
};

}

/** This class holds the set of ThreeAtomFunction potentials that
    act between the atoms in a molecule
    
    @author Christopher Woods
*/
class SIREMM_EXPORT ThreeAtomFunctions
        : public SireBase::ConcreteProperty<ThreeAtomFunctions,AtomFunctions>
{

friend QDataStream& ::operator<<(QDataStream&, const ThreeAtomFunctions&);
friend QDataStream& ::operator>>(QDataStream&, ThreeAtomFunctions&);

public:
    ThreeAtomFunctions();
    
    ThreeAtomFunctions(const MoleculeData &moldata);
    
    ThreeAtomFunctions(const ThreeAtomFunctions &other);
    
    ~ThreeAtomFunctions();
    
    static const char* typeName()
    {
        return QMetaType::typeName( qMetaTypeId<ThreeAtomFunctions>() );
    }
    
    ThreeAtomFunctions& operator=(const ThreeAtomFunctions &other);
    
    bool operator==(const ThreeAtomFunctions &other) const;
    bool operator!=(const ThreeAtomFunctions &other) const;

    void set(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2,
             const Expression &expression);
             
    void set(const AtomID &atom0, const AtomID &atom1,
             const AtomID &atom2,
             const Expression &expression);
             
    void set(const AngleID &angleid, const Expression &expression);

    void clear(AtomIdx atom);
    void clear(const AtomID &atom);

    void clear(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2);
    void clear(const AtomID &atom0, const AtomID &atom1,
               const AtomID &atom2);
    void clear(const AngleID &angleid);
    
    void clear();

    void substitute(const Identities &identities);

    Expression potential(AtomIdx atom0, AtomIdx atom1,
                         AtomIdx atom2) const;
    Expression potential(const AtomID &atom0, const AtomID &atom1,
                         const AtomID &atom2) const;
    Expression potential(const AngleID &angleid) const;
    
    Expression force(AtomIdx atom0, AtomIdx atom1, 
                     AtomIdx atom2, const Symbol &symbol) const;
    Expression force(const AtomID &atom0, const AtomID &atom1,
                     const AtomID &atom2, const Symbol &symbol) const;
    Expression force(const AngleID &angleid, const Symbol &symbol) const;

    QVector<ThreeAtomFunction> potentials() const;
    QVector<ThreeAtomFunction> forces(const Symbol &symbol) const;
    
private:
    void removeSymbols(QSet<Symbol> symbols);

    /** All of the potential functions, identified by the atom triple
        that contains that function */
    QHash<detail::IDTriple,Expression> potentials_by_atoms;
};

//////
////// Inline functions of ThreeAtomFunction
//////

/** Return the first atom of the triple */
inline const CGAtomIdx& ThreeAtomFunction::atom0() const
{
    return atm0;
}

/** Return the second atom of the triple */
inline const CGAtomIdx& ThreeAtomFunction::atom1() const
{
    return atm1;
}

/** Return the third atom of the triple */
inline const CGAtomIdx& ThreeAtomFunction::atom2() const
{
    return atm2;
}

}

Q_DECLARE_METATYPE( SireMM::ThreeAtomFunctions );

#endif
