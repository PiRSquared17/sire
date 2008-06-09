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

#ifndef SIREMOL_WEIGHTFUNCTION_H
#define SIREMOL_WEIGHTFUNCTION_H

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class WeightFunction;
class WeightFunc;

class RelFromMass;
class RelFromNumber;
class AbsFromMass;
class AbsFromNumber;
}

QDataStream& operator<<(QDataStream&, const SireMol::WeightFunction&);
QDataStream& operator>>(QDataStream&, SireMol::WeightFunction&);

QDataStream& operator<<(QDataStream&, const SireMol::WeightFunc&);
QDataStream& operator>>(QDataStream&, SireMol::WeightFunc&);

QDataStream& operator<<(QDataStream&, const SireMol::AbsFromMass&);
QDataStream& operator>>(QDataStream&, SireMol::AbsFromMass&);

QDataStream& operator<<(QDataStream&, const SireMol::RelFromMass&);
QDataStream& operator>>(QDataStream&, SireMol::RelFromMass&);

QDataStream& operator<<(QDataStream&, const SireMol::AbsFromNumber&);
QDataStream& operator>>(QDataStream&, SireMol::AbsFromNumber&);

QDataStream& operator<<(QDataStream&, const SireMol::RelFromNumber&);
QDataStream& operator>>(QDataStream&, SireMol::RelFromNumber&);

namespace SireMol
{

class AtomSelection;
class MoleculeData;
class MoleculeView;

using SireBase::PropertyMap;

/** This is the base class of all weight functions. A weight function
    is a simple function that takes two groups in a molecule, and
    returns the relative weight of those two groups (0 == 100% group A,
    0.5 == 50% group A, 50% group B, 1 == 100% groupB)

    @author Christopher Woods
*/
class SIREMOL_EXPORT WeightFunc : public SireBase::PropertyBase
{

friend QDataStream& ::operator<<(QDataStream&, const WeightFunc&);
friend QDataStream& ::operator>>(QDataStream&, WeightFunc&);

public:
    WeightFunc();

    WeightFunc(const WeightFunc &other);

    virtual ~WeightFunc();

    static const char* typeName()
    {
        return "SireMol::WeightFunc";
    }

    virtual WeightFunc* clone() const=0;

    /** Return the relative weight of group0 and group1 in the
        molecule whose data is in 'moldata', using the supplied
        map to find the required properties

        \throw SireError::incompatible_error
        \throw SireBase::missing_property
    */
    virtual double operator()(const MoleculeData &moldata,
                              const AtomSelection &group0,
                              const AtomSelection &group1,
                              const PropertyMap &map = PropertyMap()) const=0;

    /** Return the relative weight of two molecule views (view0
        and view1), using map0 to find the required properties
        from view0, and map1 to find the required properties from
        view1.

        \throw SireBase::missing_property
    */
    virtual double operator()(const MoleculeView &view0,
                              const PropertyMap &map0,
                              const MoleculeView &view1,
                              const PropertyMap &map1) const=0;

    /** Return the relative weight of the two molecule views
        (view0 and view1) using the supplied map to find the
        required properties from both views

        \throw SireBase::missing_property
    */
    double operator()(const MoleculeView &view0,
                      const MoleculeView &view1,
                      const PropertyMap &map = PropertyMap()) const
    {
        return this->operator()(view0, map, view1, map);
    }
};

/** This is the wrapper class that holds a generic WeightFunc
    Property

    @author Christopher Woods
*/
class SIREMOL_EXPORT WeightFunction : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const WeightFunction&);
friend QDataStream& ::operator>>(QDataStream&, WeightFunction&);

public:
    WeightFunction();
    WeightFunction(const SireBase::PropertyBase &property);
    WeightFunction(const WeightFunc &molgroup);

    WeightFunction(const WeightFunction &other);

    ~WeightFunction();

    virtual WeightFunction& operator=(const SireBase::PropertyBase &property);
    virtual WeightFunction& operator=(const WeightFunc &other);

    double operator()(const MoleculeData &moldata,
                      const AtomSelection &group0,
                      const AtomSelection &group1,
                      const PropertyMap &map = PropertyMap()) const;

    double operator()(const MoleculeView &view0,
                      const PropertyMap &map0,
                      const MoleculeView &view1,
                      const PropertyMap &map1) const;

    double operator()(const MoleculeView &view0,
                      const MoleculeView &view1,
                      const PropertyMap &map = PropertyMap()) const;

    const WeightFunc* operator->() const;
    const WeightFunc& operator*() const;

    const WeightFunc& read() const;
    WeightFunc& edit();

    const WeightFunc* data() const;
    const WeightFunc* constData() const;

    WeightFunc* data();

    operator const WeightFunc&() const;

    static const WeightFunction& shared_null();
};

/** This class calculates the weight by assigning all of the weight to
    the group with the largest number of atoms.

    @author Christopher Woods
*/
class SIREMOL_EXPORT AbsFromNumber
            : public SireBase::ConcreteProperty<AbsFromNumber,WeightFunc>
{

friend QDataStream& ::operator<<(QDataStream&, const AbsFromNumber&);
friend QDataStream& ::operator>>(QDataStream&, AbsFromNumber&);

public:
    AbsFromNumber();
    AbsFromNumber(const AbsFromNumber &other);

    ~AbsFromNumber();

    static const char* typeName()
    {
        return "SireMol::AbsFromNumber";
    }

    using SireBase::PropertyBase::operator=;
    using SireBase::PropertyBase::operator==;
    using SireBase::PropertyBase::operator!=;

    AbsFromNumber& operator=(const AbsFromNumber&)
    {
        return *this;
    }

    bool operator==(const AbsFromNumber&) const
    {
        return true;
    }

    bool operator!=(const AbsFromNumber&) const
    {
        return false;
    }

    AbsFromNumber* clone() const
    {
        return new AbsFromNumber(*this);
    }

    double operator()(const MoleculeData &moldata,
                      const AtomSelection &group0,
                      const AtomSelection &group1,
                      const PropertyMap &map = PropertyMap()) const;

    double operator()(const MoleculeView &view0,
                      const PropertyMap &map0,
                      const MoleculeView &view1,
                      const PropertyMap &map1) const;

private:
    double weight(int nats0, int nats1) const;
};

/** This class calculates the weight by assigning the weight based on the
    ratio of the number of atoms in the two groups.

    @author Christopher Woods
*/
class SIREMOL_EXPORT RelFromNumber
            : public SireBase::ConcreteProperty<RelFromNumber,WeightFunc>
{

friend QDataStream& ::operator<<(QDataStream&, const RelFromNumber&);
friend QDataStream& ::operator>>(QDataStream&, RelFromNumber&);

public:
    RelFromNumber();
    RelFromNumber(const RelFromNumber &other);

    ~RelFromNumber();

    static const char* typeName()
    {
        return "SireMol::RelFromNumber";
    }

    using SireBase::PropertyBase::operator=;
    using SireBase::PropertyBase::operator==;
    using SireBase::PropertyBase::operator!=;

    RelFromNumber& operator=(const RelFromNumber&)
    {
        return *this;
    }

    bool operator==(const RelFromNumber&) const
    {
        return true;
    }

    bool operator!=(const RelFromNumber&) const
    {
        return false;
    }

    RelFromNumber* clone() const
    {
        return new RelFromNumber(*this);
    }

    double operator()(const MoleculeData &moldata,
                      const AtomSelection &group0,
                      const AtomSelection &group1,
                      const PropertyMap &map = PropertyMap()) const;

    double operator()(const MoleculeView &view0,
                      const PropertyMap &map0,
                      const MoleculeView &view1,
                      const PropertyMap &map1) const;

private:
    double weight(int nats0, int nats1) const;
};

/** This class calculates the weight by assigning all of the weight to
    the group with the most mass

    @author Christopher Woods
*/
class SIREMOL_EXPORT AbsFromMass
            : public SireBase::ConcreteProperty<AbsFromMass,WeightFunc>
{

friend QDataStream& ::operator<<(QDataStream&, const AbsFromMass&);
friend QDataStream& ::operator>>(QDataStream&, AbsFromMass&);

public:
    AbsFromMass();
    AbsFromMass(const AbsFromMass &other);

    ~AbsFromMass();

    static const char* typeName()
    {
        return "SireMol::AbsFromMass";
    }

    using SireBase::PropertyBase::operator=;
    using SireBase::PropertyBase::operator==;
    using SireBase::PropertyBase::operator!=;

    AbsFromMass& operator=(const AbsFromMass&)
    {
        return *this;
    }

    bool operator==(const AbsFromMass&) const
    {
        return true;
    }

    bool operator!=(const AbsFromMass&) const
    {
        return false;
    }

    AbsFromMass* clone() const
    {
        return new AbsFromMass(*this);
    }

    double operator()(const MoleculeData &moldata,
                      const AtomSelection &group0,
                      const AtomSelection &group1,
                      const PropertyMap &map = PropertyMap()) const;

    double operator()(const MoleculeView &view0,
                      const PropertyMap &map0,
                      const MoleculeView &view1,
                      const PropertyMap &map1) const;

private:
    double weight(double mass0, double mass1) const;
};

/** This class calculates the weight by assigning the weight based on the
    ratio of the masses of the two groups.

    @author Christopher Woods
*/
class SIREMOL_EXPORT RelFromMass
            : public SireBase::ConcreteProperty<RelFromMass,WeightFunc>
{

friend QDataStream& ::operator<<(QDataStream&, const RelFromMass&);
friend QDataStream& ::operator>>(QDataStream&, RelFromMass&);

public:
    RelFromMass();
    RelFromMass(const RelFromMass &other);

    ~RelFromMass();

    static const char* typeName()
    {
        return "SireMol::RelFromMass";
    }

    using SireBase::PropertyBase::operator=;
    using SireBase::PropertyBase::operator==;
    using SireBase::PropertyBase::operator!=;

    RelFromMass& operator=(const RelFromMass&)
    {
        return *this;
    }

    bool operator==(const RelFromMass&) const
    {
        return true;
    }

    bool operator!=(const RelFromMass&) const
    {
        return false;
    }

    RelFromMass* clone() const
    {
        return new RelFromMass(*this);
    }

    double operator()(const MoleculeData &moldata,
                      const AtomSelection &group0,
                      const AtomSelection &group1,
                      const PropertyMap &map = PropertyMap()) const;

    double operator()(const MoleculeView &view0,
                      const PropertyMap &map0,
                      const MoleculeView &view1,
                      const PropertyMap &map1) const;

private:
    double weight(double mass0, double mass1) const;
};

}

Q_DECLARE_METATYPE(SireMol::WeightFunction);
Q_DECLARE_METATYPE(SireMol::RelFromMass);
Q_DECLARE_METATYPE(SireMol::RelFromNumber);
Q_DECLARE_METATYPE(SireMol::AbsFromMass);
Q_DECLARE_METATYPE(SireMol::AbsFromNumber);

SIRE_EXPOSE_CLASS( SireMol::WeightFunc )
SIRE_EXPOSE_CLASS( SireMol::RelFromMass )
SIRE_EXPOSE_CLASS( SireMol::RelFromNumber )
SIRE_EXPOSE_CLASS( SireMol::AbsFromMass )
SIRE_EXPOSE_CLASS( SireMol::AbsFromNumber )

SIRE_EXPOSE_PROPERTY( SireMol::WeightFunction, SireMol::WeightFunc )

SIRE_END_HEADER

#endif
