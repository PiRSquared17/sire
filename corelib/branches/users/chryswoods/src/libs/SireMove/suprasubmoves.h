/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIREMOVE_SUPRASUBMOVES_H
#define SIREMOVE_SUPRASUBMOVES_H

#include "suprasubmove.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class SupraSubMoves;
class SameSupraSubMoves;
}

QDataStream& operator<<(QDataStream&, const SireMove::SupraSubMoves&);
QDataStream& operator>>(QDataStream&, SireMove::SupraSubMoves&);

QDataStream& operator<<(QDataStream&, const SireMove::SameSupraSubMoves&);
QDataStream& operator>>(QDataStream&, SireMove::SameSupraSubMoves&);

namespace SireMove
{

/** This is the base class of holders of collections of SupraSubMove objects

    @author Christopher Woods
*/
class SIREMOVE_EXPORT SupraSubMoves : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const SupraSubMoves&);
friend QDataStream& ::operator>>(QDataStream&, SupraSubMoves&);

public:
    SupraSubMoves();
    
    SupraSubMoves(const SupraSubMoves &other);
    
    virtual ~SupraSubMoves();
    
    virtual SupraSubMoves* clone() const=0;
    
    static const char* typeName()
    {
        return "SireMove::SupraSubMoves";
    }
    
    const SupraSubMove& operator[](int i) const;
    
    int count() const;
    int size() const;
    int nSubMoveTypes() const;
    
    virtual QString toString() const=0;
    
    virtual void move(SupraSubSystem &system, int nsubmoves,
                      bool record_substats) const=0;
                      
    virtual void clearStatistics()=0;
    
    virtual QList<SupraSubMovePtr> subMoves() const=0;
    
    static const SameSupraSubMoves& null();
};

/** This is a simple SupraSubMoves class that just repeats the 
    same SupraSubMove multiple times
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT SameSupraSubMoves
        : public SireBase::ConcreteProperty<SameSupraSubMoves,SupraSubMoves>
{

QDataStream& ::operator<<(QDataStream&, const SameSupraSubMoves&);
QDataStream& ::operator>>(QDataStream&, SameSupraSubMoves&);

public:
    SameSupraSubMoves();
    SameSupraSubMoves(const SupraSubMove &move);
    
    SameSupraSubMoves(const SameSupraSubMoves &other);
    
    ~SameSupraSubMoves();
    
    SameSupraSubMoves& operator=(const SameSupraSubMoves &other);
    
    bool operator==(const SameSupraSubMoves &other) const;
    bool operator!=(const SameSupraSubMoves &other) const;
    
    static const char* typeName()
    {
        return QMetaType::typeName( qMetaTypeId<SameSupraSubMoves>() );
    }
    
    QString toString() const;
    
    void move(SupraSubSystem &system, int nsubmoves,
              bool record_substats) const

    void clearStatistics();
    
    QList<SupraSubMove> subMoves() const;

private:
    /** The move that will be repeatedly applied to the SupraSubSystem */
    SupraSubMove submove;
};

typedef SireBase::PropPtr<SupraSubMoves> SupraSubMovesPtr;

}

Q_DECLARE_METATYPE( SireMove::SameSupraSubMoves )

SIRE_EXPOSE_CLASS( SireMove::SupraSubMoves )
SIRE_EXPOSE_CLASS( SireMove::SameSupraSubMoves )

SIRE_EXPOSE_PROPERTY( SireMove::SupraSubMovesPtr, SireMove::SupraSubMoves )

SIRE_END_HEADER

#endif