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

#ifndef SIREN_OBJREF_H
#define SIREN_OBJREF_H

#include "object.h"
#include "detail/sharedpolypointer.hpp"

SIREN_BEGIN_HEADER

namespace Siren
{
class ObjRef;
class Stream;

template<class T> class ObjPtr;
}

namespace Siren
{

namespace detail
{
template<class T> struct StreamHelper;
}

Object* extractPointer(ObjRef &objref);

/** This is a light-weight reference to an object. This class is 
    used when a new object is returned from a function. This class
    will automatically cast itself back into a real const Object& 
    reference
*/
class SIREN_EXPORT ObjRef
{
public:
    ObjRef();
    ObjRef(const Object &obj);
    ObjRef(Object *obj);
    
    template<class T>
    ObjRef(const ObjPtr<T> &ptr);
    
    ObjRef(const ObjRef &other);
    
    ~ObjRef();
    
    ObjRef& operator=(const ObjRef &other);
    ObjRef& operator=(const Object &other);
    
    bool operator==(const ObjRef &other) const;
    bool operator!=(const ObjRef &other) const;

    bool operator==(const Object &other) const;
    bool operator!=(const Object &other) const;

    Class getClass() const;

    static QString typeName();

    QString what() const;

    ObjRef clone() const;
    
    bool unique() const;
    bool isNone() const;
    
    QString toString() const;

    bool test() const;
    bool test(Logger &logger) const;

    void copy(const Object &other);

    bool equals(const Object &other) const;
    
    void save(Stream &s) const;
    void load(Stream &s);
    
    void stream(Stream &s);
    
    uint hashCode() const;

    template<class T>
    bool isA() const;

    template<class T>
    const T& asA() const;
    
    const Object& operator*() const;
    
    operator const Object&() const;
    
private:
    friend class detail::StreamHelper<ObjRef>;

    friend Object* extractPointer(ObjRef &objref);

    /** The shared pointer to the actual object */
    detail::SharedPolyPointer<Object> d;
};

#ifndef SIREN_SKIP_INLINE_FUNCTIONS

template<class T>
SIREN_OUTOFLINE_TEMPLATE
ObjRef::ObjRef( const ObjPtr<T> &ptr )
{
    if (not ptr.isNull())
        this->operator=( *ptr );
}

/** Return whether or not this object is an object of type 'T' */
template<class T>
SIREN_INLINE_TEMPLATE
bool ObjRef::isA() const
{
    return d->isA<T>();
}

/** Cast this object to type 'T'. This raises a Siren::invalid_cast
    exception if this cast is not possible 
    
    \throw Siren::invalid_cast
*/
template<class T>
SIREN_INLINE_TEMPLATE
const T& ObjRef::asA() const
{
    return d->asA<T>();
}

#endif // SIREN_SKIP_INLINE_FUNCTIONS

}

SIREN_EXPOSE_CLASS( Siren::ObjRef )

SIREN_END_HEADER

#endif