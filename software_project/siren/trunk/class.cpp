/********************************************\
  *
  *  Siren - C++ metaobject library
  *
  *  Copyright (C) 2011  Christopher Woods
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

#include "Siren/class.h"
#include "Siren/none.h"
#include "Siren/obj.h"
#include "Siren/string.h"
#include "Siren/bytearray.h"
#include "Siren/siren.hpp"

using namespace Siren;

REGISTER_SIREN_CLASS( Siren::Class )

/** Return a full list of the registered Siren classes in this program */
StringList Class::registeredTypes()
{
    return detail::ClassRegistry::registeredClasses();
}

/** Null constructor */
Class::Class() : super(), d(0)
{}

/** Construct the class handle for the passed object */
Class::Class(const Object &object) : super(), d(0)
{
    Class::copy_object(object.getClass());
}

/** Construct the class handle for the passed class metadata */
Class::Class(const detail::ClassData *data) : super(), d(data)
{}

/** Construct the class handle for the class with specified name */
Class::Class(const char *type_name) : super(), d(0)
{
    Class::copy_object( detail::ClassRegistry::getClass(type_name) );
}

/** Construct the class handle for the class with specified name */
Class::Class(const String &type_name) : super(), d(0)
{
    Class::copy_object( detail::ClassRegistry::getClass(type_name) );
}

/** Copy constructor */
Class::Class(const Class &other) : super(other), d(other.d)
{} 

/** Destructor */
Class::~Class()
{}

/** Copy assignment operator */
void Class::copy_object(const Class &other)
{
    d = other.d;
    super::copy_object(other);
}

/** Comparison operator */
bool Class::compare_object(const Class &other) const
{
    return d == other.d and super::compare_object(other);
}

/** Return a string representation of this class */
String Class::toString() const
{
    StringList ifaces = this->interfaces();
    
    if (this->hasSuper())
    {
        if (ifaces.isEmpty())
            return String::tr("Class{ \"%1\", Base = \"%2\" }")
                    .arg(this->name(), this->superClass().name());
        else
            return String::tr("Class{ \"%1\", Base = \"%2\", Interfaces = ( %3 ) }")
                    .arg(this->name(), this->superClass().name(),
                         this->interfaces().join(", "));
    }
    else
    {
        if (ifaces.isEmpty())
            return String::tr("Class{ \"%1\" }").arg(this->name());
        else
            return String::tr("Class{ \"%1\", Interfaces = ( %3 ) }")
                    .arg(this->name(), this->interfaces().join(", "));
    }
}

/** Check that the class is not null. This will ensure that
    a null Class describes the Siren::None class */
void Class::checkNotNull() const
{
    if (d == 0)
    {
        None none;
    
        const_cast<Class*>(this)->operator=(none.getClass());
    }
}

/** Return the name of the class */
String Class::name() const
{
    checkNotNull();
    return d->typeName();
}

/** Return whether or not the class has a super class (should do!) */
bool Class::hasSuper() const
{
    checkNotNull();
    return d->hasSuper();
}

/** Return the class object that represents the superclass of this class */
Class Class::superClass() const
{
    checkNotNull();
    return Class( d->baseTypeName() );
}

/** Return whether or not this class inherits from the type "type_name",
    or implements the interface with that type */
bool Class::implements(const String &type_name) const
{
    if (this->name() == type_name)
        return true;
        
    else if (this->interfaces().contains(type_name))
        return true;

    else if (this->hasSuper())
        return this->superClass().implements(type_name);
    
    else
        return false;
}

/** Return whether or not this class inherits from the type "type_name",
    or implements the interface with that type */
bool Class::implements(const char* type_name) const
{
    return this->implements( String(type_name) );
}

/** Return the list of interfaces supported by this class */
StringList Class::interfaces() const
{
    return d->interfaces();
}

/** Return whether or not this class is a concrete type */
bool Class::isConcrete() const
{
    checkNotNull();
    return d->isConcrete();
}

/** Return whether or not objects of this type can be created */
bool Class::canCreate() const
{
    return this->isConcrete();
}

/** Return a newly created object of this class */
Obj Class::create() const
{
    if (this->isConcrete())
        return Obj( d->createObject() );
    else
    {
        sirenDebug() << "CANNOT PRODUCE OBJECT OF VIRTUAL TYPE" << this->name();
        return Obj();
    }
}

/** Return whether or not this type can be cast to the class "type_name" */
bool Class::canCast(const String &type_name) const
{
    return this->implements(type_name);
}

/** Return whether or not this type can be cast to the class "type_name" */
bool Class::canCast(const char *type_name) const
{
    return this->implements(type_name);
}

/** Assert that is is possible to cast this type to the class "type_name" */
void Class::assertCanCast(const String &type_name) const
{
    if (not this->canCast(type_name))
        sirenDebug() << "Cannot cast from type" << this->name() << "to" << type_name;
}

/** Assert that it is possible to cast this type to the class "type_name" */
void Class::assertCanCast(const char *type_name) const
{
    assertCanCast( String(type_name) );
}