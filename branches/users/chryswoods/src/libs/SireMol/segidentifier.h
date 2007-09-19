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

#ifndef SIREMOL_SEGIDENTIFIER_H
#define SIREMOL_SEGIDENTIFIER_H

#include "SireID/identifier.h"

#include "segid.h"
#include "segname.h"

namespace SireMol
{

class SegIdentifier : public SireID::Identifier_T_<SegIdentifier,SegID>
{
public:
    SegIdentifier() : SireID::Identifier_T_<SegIdentifier,SegID>()
    {}
    
    SegIdentifier(const SegID &segid) 
            : SireID::Identifier_T_<SegIdentifier,SegID>(segid)
    {}
    
    SegIdentifier(const SegIdentifier &other) 
            : SireID::Identifier_T_<SegIdentifier,SegID>(other)
    {}
    
    ~SegIdentifier()
    {}
    
    SegName map(const MoleculeInfo &molinfo) const
    {
        this->assertNotNull();
        return this->asA<SegID>().map(molinfo);
    }
};

}

QDataStream& operator<<(QDataStream&, const SireMol::SegIdentifier&);
QDataStream& operator>>(QDataStream&, SireMol::SegIdentifier&);

XMLStream& operator<<(XMLStream&, const SireMol::SegIdentifier&);
XMLStream& operator>>(XMLStream&, SireMol::SegIdentifier&);

uint qHash(const SireMol::SegIdentifier &segid);

Q_DECLARE_METATYPE(SireMol::SegIdentifier);

#endif

