/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#ifndef SIREIO_DX_H
#define SIREIO_DX_H

#include "SireVol/gridinfo.h"

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireIO
{
class DX;
}

QDataStream& operator<<(QDataStream&, const SireIO::DX&);
QDataStream& operator>>(QDataStream&, SireIO::DX&);

namespace SireMol
{
class VolumeMap;
}

namespace SireIO
{

using SireVol::GridInfo;
using SireBase::PropertyMap;
using SireMol::VolumeMap;

/** This IO class is used to read and write grid data in DX format.
    Note that this has only been written to support basic, cubic/rectangular
    grids, and doesn't support the whole file format!
    
    @author Christopher Woods
*/
class SIREIO_EXPORT DX : public SireBase::ConcreteProperty<DX,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const DX&);
friend QDataStream& ::operator>>(QDataStream&, DX&);

public:
    DX();
    DX(const DX &other);
    
    ~DX();
    
    DX& operator=(const DX &other);
    
    bool operator==(const DX &other) const;
    bool operator!=(const DX &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;
    
    void write(const GridInfo &gridinfo, const QVector<float> &values,
               const QString &filename, const PropertyMap &map = PropertyMap()) const;

    void write(const VolumeMap &volumemap, const QString &filename,
               const PropertyMap &map = PropertyMap()) const;
};

}

Q_DECLARE_METATYPE( SireIO::DX )

SIRE_EXPOSE_CLASS( SireIO::DX )

SIRE_END_HEADER

#endif
