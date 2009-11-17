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

#ifndef SIREN_PRIMITIVE_H
#define SIREN_PRIMITIVE_H

#include "objref.h"
#include "stream.h"

SIREN_BEGIN_HEADER

namespace Siren
{
class String;
class Number;
}

namespace Siren
{

class Logger;

namespace detail
{

template<class T>
struct PrimitiveHelper
{
    static QString typeName();
    
    static QString toString(const T &object);
    
    static HASH_CODE hashCode(const T &object);
    
    static void stream(T &object, Stream &s);
    
    static bool test(const T &object, Logger &logger);
};

} // end of namespace detail

/** This is the base class of Primitive types. This just
    provides the functionality to allow automatic
    conversion of Object types. Primitives are simple
    classes that are too small or performance/memory 
    intensive to be made virtual or have the costs
    associated with being part of the Object hierarchy.
    
    This class provides the functionality to allow 
    Primitive objects to be automatically converted
    to PrimitiveObject<T> types so that they can
    be used as Object types when necessary.
    
    @author Christopher Woods
*/
template<class Derived>
class SIREN_EXPORT Primitive
{
public:
    Primitive();
    ~Primitive();
    
    operator ObjRef() const;
    
    void save(Stream &s) const;
    void load(Stream &s);
};

/** This is a simple wrapper that provides a Siren::Object
    interface to Primitive objects 
*/
template<class T>
class SIREN_EXPORT PrimitiveObject : public Implements<PrimitiveObject<T>,Object>
{
public:
    PrimitiveObject();
    PrimitiveObject(const T &primitive);
    PrimitiveObject(const PrimitiveObject<T> &other);
    
    ~PrimitiveObject();
    
    PrimitiveObject& operator=(const PrimitiveObject<T> &other);
    PrimitiveObject& operator=(const T &primitive);
    
    bool operator==(const PrimitiveObject<T> &other) const;
    bool operator!=(const PrimitiveObject<T> &other) const;
    
    bool operator==(const T &primitive) const;
    bool operator!=(const T &primitive) const;
    
    QString toString() const;
    
    HASH_CODE hashCode() const;

    void stream(Stream &s);
    
    bool test(Logger &logger) const;
    
    operator const T&() const;

private:
    /** The actual primitive object */
    T primitive_object;
};

/** This is a simple wrapper class that signifies that
    the passed primitive is a String. This is used with
    functions that expect a generic Object reference, and
    allows the user to write;
    
    object.func( obj0, String("Hello") );
    
    void object::func( const Object &obj0, const Object &obj1 )
    {
        QString s = primitive_cast<String>(obj1);
    }
*/
class SIREN_EXPORT String : public Primitive<String>
{
public:
    String();
    String(const QString &text);
    
    String(const String &other);
    
    ~String();
    
    String& operator=(const String &other);
    
    bool operator==(const String &other) const;
    bool operator!=(const String &other) const;
    
    operator const QString&() const;
    
    static QString typeName();
    
    QString toString() const;
    
    HASH_CODE hashCode() const;

    void stream(Stream &s);
    
    bool test() const;
    bool test(Logger &logger) const;
    
private:
    /** The actual text in the string */
    QString text;
};

/** This is a simple wrapper class that signifies that 
    the passed primitive is a Number. This is used with
    functions that expect a generic Object reference, and
    allows the user to write;
    
    object.func( obj0, Number(42) );
    
    void object::func( const Object &obj0, const Object &obj1 )
    {
        Number n = primitive_cast<Number>(obj1);
        
        int val = primitive_cast<Number>(obj1).toInteger();
    }
*/
class SIREN_EXPORT Number : public Primitive<Number>
{
public:
    Number();
    Number(qint64 value);
    Number(double value);
    
    Number(const Number &other);
    
    ~Number();
    
    Number& operator=(const Number &other);
    
    bool operator==(const Number &other) const;
    bool operator!=(const Number &other) const;
    
    static QString typeName();
    
    QString toString() const;
    
    HASH_CODE hashCode() const;
    
    bool test() const;
    bool test(Logger &logger) const;
    
    static Number fromInteger(qint64 value);
    static Number fromFloat(double value);

    qint64 toInteger() const;
    double toFloat() const;

    void stream(Stream &s);
    
    operator qint64() const;
    operator double() const;

private:
    union
    {
        qint64 integer_value;
        double float_value;
    
    } number_data;
    
    /** Whether or not this is an integer value */
    bool is_integer;
};

#ifndef SIREN_SKIP_INLINE_FUNCTIONS

/////////
///////// Implementation of PrimitiveHelper
/////////

namespace detail
{

template<class T>
SIREN_OUTOFLINE_TEMPLATE
QString PrimitiveHelper<T>::typeName()
{
    return T::typeName();
}
    
template<class T>
SIREN_OUTOFLINE_TEMPLATE
QString PrimitiveHelper<T>::toString(const T &object)
{
    return object.toString();
}
    
template<class T>
SIREN_OUTOFLINE_TEMPLATE
HASH_CODE PrimitiveHelper<T>::hashCode(const T &object)
{
    return object.hashCode();
}
    
template<class T>
SIREN_OUTOFLINE_TEMPLATE
void PrimitiveHelper<T>::stream(T &object, Stream &s)
{
    t.stream(s);
}
    
template<class T>
SIREN_OUTOFLINE_TEMPLATE
bool PrimitiveHelper<T>::test(const T &object, Logger &logger)
{
    return object.test(logger);
}

} // end of namespace detail

/////////
///////// Implementation of Primitive
/////////

/** Constructor */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
Primitive<T>::Primitive()
{}

/** Destructor */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
Primitive<T>::~Primitive()
{}
  
/** Allow automatic casting to an Object reference */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
Primitive<T>::operator ObjRef() const
{
    return ObjRef( PrimitiveObject<T>(static_cast<const T&>(*this)) );
}

/** Save this object to the passed stream */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
void Primitive<T>::save(Stream &s) const
{
    s.assertIsSaving();
    T *t_ptr = const_cast<T*>( static_cast<const T*>(this) );
    t_ptr->stream(s);
}

/** Load this object from the passed stream */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
void Primitive<T>::load(Stream &s)
{
    s.assertIsLoading();
    T *t_ptr = static_cast<T*>(this);
    t_ptr->stream(s);
}

/////////
///////// Implementation of PrimitiveObject
/////////

/** Default constructor */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
PrimitiveObject<T>::PrimitiveObject() : Implements<PrimitiveObject<T>,Object>()
{}

/** Construct as a copy of 'primitive' */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
PrimitiveObject<T>::PrimitiveObject(const T &primitive)
                   : Implements<PrimitiveObject<T>,Object>(), primitive_object(primitive)
{}

/** Copy constructor */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
PrimitiveObject<T>::PrimitiveObject(const PrimitiveObject<T> &other)
                   : Implements<PrimitiveObject<T>,Object>(other),
                     primitive_object(other.primitive_object)
{}

/** Destructor */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
PrimitiveObject<T>::~PrimitiveObject()
{}

/** Copy assignment operator */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
PrimitiveObject<T>& PrimitiveObject<T>::operator=(const PrimitiveObject<T> &other)
{
    if (this != &other)
    {
        primitive_object = other.primitive_object;
        Object::operator=(other);
    }
    
    return *this;
}

/** Copy assignment operator */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
PrimitiveObject<T>& PrimitiveObject<T>::operator=(const T &primitive)
{
    primitive_object = primitive;
    return *this;
}

/** Comparison operator */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
bool PrimitiveObject<T>::operator==(const PrimitiveObject<T> &other) const
{
    return primitive_object == other.primitive_object;
}

/** Comparison operator */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
bool PrimitiveObject<T>::operator!=(const PrimitiveObject<T> &other) const
{
    return primitive_object != other.primitive_object;
}

/** Comparison operator */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
bool PrimitiveObject<T>::operator==(const T &primitive) const
{
    return primitive_object == primitive;
}

/** Comparison operator */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
bool PrimitiveObject<T>::operator!=(const T &primitive) const
{
    return primitive_object != primitive;
}

/** Return a string representation of this object */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
QString PrimitiveObject<T>::toString() const
{
    return detail::PrimitiveHelper<T>::toString(primitive_object);
}

/** Return a hash code for this object */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
HASH_CODE PrimitiveObject<T>::hashCode() const
{
    return detail::PrimitiveHelper<T>::hashCode(primitive_object);
}

/** Test this object */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
bool PrimitiveObject<T>::test(Logger &logger) const
{
    return detail::PrimitiveHelper<T>::test(primitive_object, logger);
}

/** Stream this object */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
void PrimitiveObject<T>::stream(Stream &s)
{
    Schema schema = s.schema<T>(1);

    if (schema.version() == 1)
    {
        detail::PrimitiveHelper<T>::stream(primitive_object, s);
    }
    else
        schema.throwUnsupportedVersion("1");

    super().stream(s);
}

/** Allow automatic casting back to the primitive */
template<class T>
SIREN_OUTOFLINE_TEMPLATE
PrimitiveObject<T>::operator const T&() const
{
    return primitive_object;
}

/** Use a 'primitive_cast' to convert an Object to the Primitive 
    of type 'T' 
    
    \throw Siren::invalid_cast
*/
template<class T>
const T& primitive_cast( const Object &object )
{
    return object.asA< PrimitiveObject<T> >();
}

#endif

typedef PrimitiveObject<String> StringObject;
typedef PrimitiveObject<Number> NumberObject;

}

Q_DECLARE_METATYPE( Siren::StringObject )
Q_DECLARE_METATYPE( Siren::NumberObject )

SIREN_EXPOSE_CLASS( Siren::String )
SIREN_EXPOSE_CLASS( Siren::Number )

SIREN_END_HEADER

#endif