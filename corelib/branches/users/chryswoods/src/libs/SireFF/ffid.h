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

#ifndef SIREFF_FFID_H
#define SIREFF_FFID_H

#include <QList>

#include "SireID/id.h"

SIRE_BEGIN_HEADER

namespace SireFF
{

class FFIdx;
class FFIdentifier;
class FFName;

class ForceFields;

/** The base class of all ForceField identifiers

    @author Christopher Woods
*/
class SIREFF_EXPORT FFID : public SireID::ID
{
public:
    typedef FFIdx Index;
    typedef FFIdentifier Identifier;

    FFID();
    FFID(const FFID &other);

    virtual ~FFID();
    
    static const char* typeName()
    {
        return "SireFF::FFID";
    }
    
    virtual FFID* clone() const=0;
    
    virtual QList<FFIdx> map(const ForceFields &ffields) const=0;
};

}

SIRE_EXPOSE_CLASS( SireFF::FFID )

SIRE_END_HEADER

#endif
