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

#ifndef SIREID_NAME_H
#define SIREID_NAME_H

#include <QHash>
#include <QString>

#include "id.h"

#include "Siren/stream.h"

SIRE_BEGIN_HEADER

namespace SireID
{

enum CaseSensitivity{  CaseInsensitive = 0,
                       CaseSensitive   = 1 };

/** This is the base class of all Name ID objects. A Name is used
    to provide an object with a human-readable name that can be
    used to identify the object, e.g. identifying atoms by their
    name within a residue, or identifying forcefields by their name.
    
    A name does not have to uniquely identify an object, though it 
    helps, as generally only the first object with a specified name
    is returned if there are in fact multiple objects with the same
    name.

    Derive from this class with code like;
    
    class AtomName : public Implements<AtomName,Name>,
                     public Interfaces<AtomName,AtomID>
    
    @author Christopher Woods
*/
class SIREID_EXPORT Name : public Siren::Extends<Name,ID>
{
public:
    Name(const QString &name = QString::null,
         CaseSensitivity = CaseSensitive);
    
    Name(const Name &other);

    ~Name();

    uint hashCode() const;
    void stream(Siren::Stream &s);

    bool isNull() const;
    
    operator QString() const;

    bool isEmpty() const;
    
    const QString& value() const;

    bool isCaseSensitive() const;

protected:
    Name& operator=(const Name &other);
    
    bool operator==(const Name &other) const;
    bool operator!=(const Name &other) const;
    
    /** The actual name */
    QString _name;
    
    /** Should this name be case sensitive or not? */
    bool case_sensitive;
};

}

SIRE_EXPOSE_CLASS( SireID::Name )

SIRE_END_HEADER

#endif
