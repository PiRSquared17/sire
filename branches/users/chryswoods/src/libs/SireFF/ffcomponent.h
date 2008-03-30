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

#ifndef SIREFF_FFCOMPONENT_H
#define SIREFF_FFCOMPONENT_H

#include "SireCAS/symbol.h"

#include "SireFF/ffname.h"

#include "SireUnits/dimensions.h"

#include <QLatin1String>

SIRE_BEGIN_HEADER

namespace SireFF
{

class FF;

class TotalComponent;

/** This is the base class of all Symbols that represent forcefield
    components.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFComponent : public SireCAS::Symbol
{
public:
    virtual ~FFComponent();
    
    static const char* typeName()
    {
        return "SireFF::FFComponent";
    }
    
    FFName forceFieldName() const;

    QString componentName() const;

    virtual const FFComponent& total() const=0;

protected:
    FFComponent(const FFName &ffname, const QLatin1String &name);
    FFComponent(const SireCAS::Symbol &symbol, const QLatin1String &name);
    
    FFComponent(const FFComponent &other);
    
    void setEnergy(FF &ff, const Symbol &symbol, double value) const;
    void changeEnergy(FF &ff, const Symbol &symbol, double delta) const;
    
    static QString symbolName(const FFName &ffname, const QLatin1String &name)
    {
        return QString("E_{%1}^{%2}").arg(ffname).arg(name);
    }
};

/** This class holds the energy for the component of type T

    @author Christopher Woods
*/
template<class T>
class SIREFF_EXPORT ComponentEnergy
{
public:
    typedef T Components;
    
    ComponentEnergy() : nrg(0)
    {}
    
    explicit ComponentEnergy(double inrg) : nrg(inrg)
    {}
    
    ComponentEnergy(const ComponentEnergy<T> &other) : nrg(other.nrg)
    {}
    
    ~ComponentEnergy()
    {}
    
    ComponentEnergy<T>& operator+=(const ComponentEnergy<T> &other)
    {
        nrg += other.nrg;
        return *this;
    }
    
    ComponentEnergy<T>& operator-=(const ComponentEnergy<T> &other)
    {
        nrg -= other.nrg;
        return *this;
    }
    
    ComponentEnergy<T> operator+(const ComponentEnergy<T> &other) const
    {
        return ComponentEnergy<T>(nrg + other.nrg);
    }
    
    ComponentEnergy<T> operator-(const ComponentEnergy<T> &other) const
    {
        return ComponentEnergy<T>(nrg - other.nrg);
    }
    
    Components components() const
    {
        return Components();
    }
    
    double component(const T&) const
    {
        return nrg;
    }
    
    double total() const
    {
        return nrg;
    }
    
    operator double() const
    {
        return nrg;
    }
    
    operator SireUnits::Dimension::Energy() const
    {
        return SireUnits::Dimension::Energy(nrg);
    }
    
private:
    /** The component of the energy */
    double nrg;
};

} // end of namespace SireFF

SIRE_END_HEADER

#endif
