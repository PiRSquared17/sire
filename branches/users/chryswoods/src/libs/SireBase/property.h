/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREBASE_PROPERTY_H
#define SIREBASE_PROPERTY_H

#include "sharedpolypointer.hpp"

#include <QVariant>

SIRE_BEGIN_HEADER

namespace SireBase
{
class PropertyBase;
class Property;
}

QDataStream& operator<<(QDataStream&, const SireBase::PropertyBase&);
QDataStream& operator>>(QDataStream&, SireBase::PropertyBase&);

QDataStream& operator<<(QDataStream&, const SireBase::Property&);
QDataStream& operator>>(QDataStream&, SireBase::Property&);

SireStream::XMLStream& operator<<(SireStream::XMLStream&, const SireBase::Property&);
SireStream::XMLStream& operator>>(SireStream::XMLStream&, SireBase::Property&);

SireStream::XMLStream& operator<<(SireStream::XMLStream&,
                                  const SireBase::PropertyBase&);
SireStream::XMLStream& operator>>(SireStream::XMLStream&,
                                  SireBase::PropertyBase&);


namespace SireBase
{

class Property;

/** This is the base class of all properties that may be attached to a
    molecule. Properties are used to assign extra information to a molecule,
    which may then be carried by the molecule throughout its passage
    through the simulation. Examples of properties may include the file
    from which the molecule was read, the charge parameters on the atoms,
    the PDB code etc.

    Properties form a polymorphic hierarchy which are implicitly shared
    via SireBase::SharedPolyPointer smart pointers.

    @author Christopher Woods
*/
class SIREBASE_EXPORT PropertyBase : public QSharedData
{

friend QDataStream& ::operator<<(QDataStream&, const PropertyBase&);
friend QDataStream& ::operator>>(QDataStream&, PropertyBase&);

public:
    PropertyBase();

    PropertyBase(const PropertyBase &other);

    virtual ~PropertyBase();

    virtual PropertyBase& operator=(const PropertyBase &other);

    PropertyBase& operator=(const Property &other);

    virtual bool operator==(const PropertyBase &other) const;
    virtual bool operator!=(const PropertyBase &other) const;

    bool operator==(const Property &other) const;
    bool operator!=(const Property &other) const;

    virtual PropertyBase* clone() const=0;
    virtual PropertyBase& copy(const PropertyBase &other)=0;
    virtual bool compare(const PropertyBase &other) const=0;

    virtual void save(QDataStream &ds) const=0;
    virtual void load(QDataStream &ds)=0;

    virtual void save(SireStream::XMLStream &xs) const=0;
    virtual void load(SireStream::XMLStream &xs);

    virtual const char* what() const=0;

    static Property null_property();

    static const char* typeName()
    {
        return "SireBase::PropertyBase";
    }

    template<class T>
    bool isA() const
    {
        return dynamic_cast<const T*>(this) != 0;
    }

    template<class T>
    const T& asA() const
    {
        return dynamic_cast<const T&>(*this);
    }
};

/** This is the second-to-top class of all Properties. Any instantiatable
    properties must inherit from this template class so that all of the
    virtual functions are supplied correctly, e.g. if you have a
    class called Derived, that inherits from Base, and Base derives
    from PropertyBase, then you need to inherit Derived from
    ConcreteProperty<Derived,Base>

    @author Christopher Woods
*/
template<class Derived, class Base>
ConcreteProperty : public Base
{
public:
    ConcreteProperty() : Base()
    {}

    template<class T0>
    ConcreteProperty(const T0 &t0) : Base(t0)
    {}

    template<class T0, class T1>
    ConcreteProperty(const T0 &t0, const T1 &t1) : Base(t0,t1)
    {}

    template<class T0, class T1, class T2>
    ConcreteProperty(const T0 &t0, const T1 &t1,
                     const T2 &t2) : Base(t0,t1,t2)
    {}

    ~ConcreteProperty()
    {}

    ConcreteProperty<Derived,Base>& operator=(const Base &other)
    {
        const Derived* other_t = dynamic_cast<const Derived*>(&other);

        if (!other_t)
            Base::throwInvalidCast(other);

        return static_cast<Derived*>(this)->operator=(*other_t);
    }

    bool operator==(const Base &other) const
    {
        const Derived* other_t = dynamic_cast<const Derived*>(&other);

        if (other_t)
            return static_cast<const Derived*>(this)->operator==(*other_t);
        else
            return false;
    }

    bool operator!=(const Base &other) const
    {
        const Derived* other_t = dynamic_cast<const Derived*>(&other);

        if (other_t)
            return static_cast<const Derived*>(this)->operator!=(*other_t);
        else
            return true;
    }

    const char* what() const
    {
        return Derived::typeName();
    }

    ConcreteProperty<Derived,Base>* clone() const
    {
        return new Derived( static_cast<const T&>(*this) );
    }

    ConcreteProperty<Derived,Base>& copy(const Base &other)
    {
        return ConcreteProperty<Derived,Base>::operator=(other);
    }

    bool compare(const Base &other) const
    {
        return ConcreteProperty<Derived,Base>::operator==(other);
    }

    void save(QDataStream &ds) const
    {
        ds << static_cast<const T&>(*this);
    }

    void load(QDataStream &ds)
    {
        ds >> static_cast<T&>(*this);
    }

protected:
    ConcreteProperty<Derived,Base>&
    operator=(const ConcreteProperty<Derived,Base> &other)
    {
        Base::operator=(other);
        return *this;
    }

    bool operator==(const ConcreteProperty<Derived,Base> &other) const
    {
        return Base::operator==(other);
    }

    bool operator!=(const ConcreteProperty<Derived,Base> &other) const
    {
        return Base::operator!=(other);
    }
};

/** This is a simple property that holds any value as a QVariant. This
    is designed to be used for metadata that doesn't need any tight
    checking (e.g. the author of the molecule file, the source of
    the coordinates, the 'header' lines etc.)

    @author Christopher Woods
*/
class SIREBASE_EXPORT VariantProperty
          : public ConcreteProperty<VariantProperty,PropertyBase>, public QVariant
{
public:
    VariantProperty();

    VariantProperty(const QVariant &value);

    VariantProperty(const Property &other);

    VariantProperty(const VariantProperty &other);

    ~VariantProperty();

    using PropertyBase::operator=();
    using PropertyBase::operator==();
    using PropertyBase::operator!=();

    VariantProperty& operator=(const QVariant &value);
    VariantProperty& operator=(const VariantProperty &other);

    bool operator==(const VariantProperty &other) const;
    bool operator!=(const VariantProperty &other) const;

    static const char* typeName()
    {
        return "SireBase::VariantProperty";
    }
};

/** This is the visible holder class for PropertyBase. This is just
    a thin wrapper around SireBase::SharedPolyPointer<PropertyBase>
    that has extra tests that ensure that this is never null.

    @author Christopher Woods
*/
class SIREBASE_EXPORT Property
{

friend QDataStream& ::operator<<(QDataStream&, const Property&);
friend QDataStream& ::operator>>(QDataStream&, Property&);

public:
    Property();

    Property(const PropertyBase &property);
    Property(PropertyBase *property);

    Property(const QVariant &value);

    Property(const Property &other);

    ~Property();

    Property& operator=(const Property &other);

    bool operator==(const Property &other) const;
    bool operator!=(const Property &other) const;

    void save(QDataStream &ds) const;
    void load(QDataStream &ds);

    const PropertyBase& base() const
    {
        return *d;
    }

    const char* what() const
    {
        return d->what();
    }

    template<class T>
    bool isA() const
    {
        return d->isA<T>();
    }

    template<class T>
    const T& asA() const
    {
        return d->asA<T>();
    }

private:
    /** Shared pointer to the actual property */
    SharedPolyPointer<PropertyBase> d;
};

}

Q_DECLARE_METATYPE(SireBase::Property);
Q_DECLARE_METATYPE(SireBase::VariantProperty);

SIRE_END_HEADER

#endif
