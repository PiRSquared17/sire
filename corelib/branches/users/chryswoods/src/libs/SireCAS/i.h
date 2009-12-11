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

#ifndef SIRECAS_I_H
#define SIRECAS_I_H

#include "constant.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{

/** This is the complex number, i

    @author Christopher Woods
*/
class SIRECAS_EXPORT I : public Siren::Implements<I,Constant>
{
public:
    I();
    I(const I &other);

    ~I();

    I& operator=(const I &other);
    
    bool operator==(const I &other) const;
    bool operator!=(const I &other) const;

    uint hashCode() const;
    QString toString() const;
    void stream(Siren::Stream &s);

    double evaluate(const Values&) const;
    SireMaths::Complex evaluate(const ComplexValues&) const;

    Expression conjugate() const;

    bool isComplex() const;
};

}

Q_DECLARE_METATYPE(SireCAS::I)

SIRE_EXPOSE_CLASS( SireCAS::I )

SIRE_END_HEADER

#endif
