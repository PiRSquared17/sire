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

#ifndef SIREMOL_RESIDENTIFIER_H
#define SIREMOL_RESIDENTIFIER_H

#include "SireID/identifier.h"

#include "resid.h"
#include "resnum.h"

namespace SireMol
{

class ResIdentifier : public SireID::Identifier_T_<ResIdentifier,ResID>
{
public:
    ResIdentifier() : SireID::Identifier_T_<ResIdentifier,ResID>()
    {}
    
    ResIdentifier(const ResID &resid) 
            : SireID::Identifier_T_<ResIdentifier,ResID>(resid)
    {}
    
    ResIdentifier(const ResIdentifier &other) 
            : SireID::Identifier_T_<ResIdentifier,ResID>(other)
    {}
    
    ~ResIdentifier()
    {}
    
    ResNum map(const MoleculeInfo &molinfo) const
    {
        this->assertNotNull();
        return this->asA<ResID>().map(molinfo);
    }
};

}

QDataStream& operator<<(QDataStream&, const SireMol::ResIdentifier&);
QDataStream& operator>>(QDataStream&, SireMol::ResIdentifier&);

XMLStream& operator<<(XMLStream&, const SireMol::ResIdentifier&);
XMLStream& operator>>(XMLStream&, SireMol::ResIdentifier&);

uint qHash(const SireMol::ResIdentifier &resid);

Q_DECLARE_METATYPE(SireMol::ResIdentifier);

#endif


