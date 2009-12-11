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

#ifndef SIRECAS_COMPLEXVALUES_H
#define SIRECAS_COMPLEXVALUES_H

#include <QHash>

#include "values.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{

/** This class holds a set of Symbols and their associated complex values. 
    This is used when numerically evaluating an equation using complex maths.

    @author Christopher Woods
*/
class SIRECAS_EXPORT ComplexValues 
            : public Siren::Implements<ComplexValues,Siren::Object>
{
public:
    ComplexValues();

    ComplexValues(const QList<SymbolComplex> &values);
    ComplexValues(const QHash<Symbol,SireMaths::Complex> &values);

    ComplexValues(const SymbolComplex &symval0);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4, const SymbolComplex &symval5);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4, const SymbolComplex &symval5,
                  const SymbolComplex &symval6);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4, const SymbolComplex &symval5,
                  const SymbolComplex &symval6, const SymbolComplex &symval7);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4, const SymbolComplex &symval5,
                  const SymbolComplex &symval6, const SymbolComplex &symval7,
                  const SymbolComplex &symval8);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4, const SymbolComplex &symval5,
                  const SymbolComplex &symval6, const SymbolComplex &symval7,
                  const SymbolComplex &symval8, const SymbolComplex &symval9);

    ComplexValues(const Values &other);

    ComplexValues(const ComplexValues &other);

    ~ComplexValues();

    ComplexValues& operator=(const ComplexValues &other);
    
    bool operator==(const ComplexValues &other) const;
    bool operator!=(const ComplexValues &other) const;
    
    QString toString() const;
    uint hashCode() const;
    void stream(Siren::Stream &s);

    void set(const Symbol &symbol, const SireMaths::Complex &value);

    void add(const SymbolComplex &symval0);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4, const SymbolComplex &symval5);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4, const SymbolComplex &symval5,
             const SymbolComplex &symval6);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4, const SymbolComplex &symval5,
             const SymbolComplex &symval6, const SymbolComplex &symval7);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4, const SymbolComplex &symval5,
             const SymbolComplex &symval6, const SymbolComplex &symval7,
             const SymbolComplex &symval8);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4, const SymbolComplex &symval5,
             const SymbolComplex &symval6, const SymbolComplex &symval7,
             const SymbolComplex &symval8, const SymbolComplex &symval9);

    SireMaths::Complex value(const Symbol &sym) const;

    const QHash<SymbolID,SireMaths::Complex>& values() const;

private:
    /** Hash mapping Symbol IDs to actual numerical values */
    QHash<SymbolID,SireMaths::Complex> vals;
};

}

Q_DECLARE_METATYPE(SireCAS::ComplexValues)

SIRE_EXPOSE_CLASS( SireCAS::ComplexValues )

SIRE_END_HEADER

#endif
