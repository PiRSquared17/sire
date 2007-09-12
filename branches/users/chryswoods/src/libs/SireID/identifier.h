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

#ifndef SIREID_IDENTIFIER_H
#define SIREID_IDENTIFIER_H

#include "id.h"

SIRE_BEGIN_HEADER

namespace SireID
{
class Identifier;
}

namespace SireID
{

/** This is a generic holder that is capable of holding any
    class that is derived from 'ID'. This is useful if a generic
    function wants to return a generic ID object.
    
    @author Christopher Woods
*/
class SIREID_EXPORT Identifier
{
public:
    Identifier();
    Identifier(const ID &idobj);
    
    Identifier(const Identifier &other);
    
    ~Identifier();
    
    Identifier& operator=(const Identifier &other);
    Identifier& operator=(const ID &idobj);
    
    bool operator==(const Identifier &other) const;
    bool operator==(const ID &idobj) const;
    
    bool operator!=(const Identifier &other) const;
    bool operator!=(const ID &idobj) const;
    
    uint hash() const;

    QString toString() const;
    
    const char* what() const;
    
    template<class T>
    bool isA() const
    {
        return idptr->isA<T>();
    }
    
    template<class T>
    const T& asA() const
    {
        return idptr->asA<T>();
    }
    
    template<class T>
    T& asA()
    {
        return idptr->asA<T>();
    }
    
private:
    /** Pointer to the actual ID object */
    ID *idptr;
}; 

}

SIRE_END_HEADER

#endif
