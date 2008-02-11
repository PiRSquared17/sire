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

#include "chainresid.h"

#include <QSet>

#include "SireMol/errors.h"

using namespace SireMol;

ChainResID::ChainResID() : ResID()
{}

ChainResID::ChainResID(const ChainID &chain, const ResID &res)
           : ResID(), chainid(chain), resid(res)
{}

ChainResID::ChainResID(const ChainResID &other)
           : ResID(other),
             chainid(other.chainid), resid(other.resid)
{}

ChainResID::~ChainResID()
{}

QString ChainResID::toString() const
{
    return QObject::tr("%1 and %2").arg(chainid.toString(), resid.toString());
}

QList<ResIdx> ChainResID::map(const MolInfo &molinfo) const
{
    if (this->isNull())
        return molinfo.getResidues();
    else if (resid.isNull())
        return molinfo.getResiduesIn(chainid);
    else if (resid.isNull())
        return resid.map(molinfo);
    
    QList<ResIdx> residxs = 
                MolInfo::intersection(resid.map(molinfo),
                                      molinfo.getResiduesIn(chainid) );
                                             
    if (residxs.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "There are no residues that match %1.")
                .arg(this->toString()), CODELOC );
            
    return residxs;
}

ChainResID SireMol::operator+(const ChainID &chainid, const ResID &resid)
{
    return ChainResID(chainid, resid);
}

ChainResID SireMol::operator+(const ResID &resid, const ChainID &chainid)
{
    return ChainResID(chainid, resid);
}
