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

#include "SireIO/dx.h"

#include "SireMol/volumemap.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include <QTextStream>
#include <QFile>

using namespace SireIO;
using namespace SireVol;
using namespace SireBase;
using namespace SireStream;

using boost::tuple;

static const RegisterMetaType<DX> r_dx;

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const DX &dx)
{
    writeHeader(ds, r_dx, 1);
    
    ds << static_cast<const Property&>(dx);
    
    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, DX &dx)
{
    VersionID v = readHeader(ds, r_dx);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(dx);
    }
    else
        throw version_error(v, "1", r_dx, CODELOC);
    
    return ds;
}

/** Constructor */
DX::DX() : ConcreteProperty<DX,Property>()
{}

/** Copy constructor */
DX::DX(const DX &other) : ConcreteProperty<DX,Property>(other)
{}

/** Destructor */
DX::~DX()
{}

/** Copy assignment operator */
DX& DX::operator=(const DX &other)
{
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool DX::operator==(const DX &other) const
{
    return Property::operator==(other);
}

/** Comparison operator */
bool DX::operator!=(const DX &other) const
{
    return not operator==(other);
}

const char* DX::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DX>() );
}

const char* DX::what() const
{
    return DX::typeName();
}

QString DX::toString() const
{
    return QObject::tr("DX()");
}

/** Write the passed grid data to the DX formatted file called 'filename' */
void DX::write(const GridInfo &gridinfo, const QVector<float> &values,
               const QString &filename, const PropertyMap &map) const
{
    if (values.isEmpty() or gridinfo.isEmpty())
        return;
    
    if (values.count() != gridinfo.nPoints())
        throw SireError::incompatible_error( QObject::tr(
                "Cannot write the grid data to file %1 as there is a disagreement "
                "in the number of points, %2 vs. %3.")
                    .arg(filename).arg(gridinfo.toString()).arg(values.count()), CODELOC );
    
    QFile f(filename);
    
    if (not f.open(QIODevice::WriteOnly))
        throw SireError::file_error(f, CODELOC);
    
    QTextStream ts(&f);
    
    ts << QString("# File written by SireIO::DX\n");
    ts << QString("object 1 class gridpositions counts %1 %2 %3\n")
                .arg(gridinfo.dimX()).arg(gridinfo.dimY()).arg(gridinfo.dimZ());
    ts << QString("origin %1 %2 %3\n")
                .arg(gridinfo.dimensions().minCoords().x())
                .arg(gridinfo.dimensions().minCoords().y())
                .arg(gridinfo.dimensions().minCoords().z());
    ts << QString("delta %1 0 0\n").arg(gridinfo.spacing().value());
    ts << QString("delta 0 %1 0\n").arg(gridinfo.spacing().value());
    ts << QString("delta 0 0 %1\n").arg(gridinfo.spacing().value());
    
    ts << QString("object 2 class gridconnections counts %1 %2 %3\n")
                .arg(gridinfo.dimX()).arg(gridinfo.dimY()).arg(gridinfo.dimZ());

    ts << QString("object 3 class array type double rank 0 items %1 data follows\n")
                .arg(gridinfo.nPoints());
    
    for (int i=1; i<=values.count(); ++i)
    {
        ts << values.constData()[i-1];
        
        if (i % 3 == 0)
            ts << "\n";
        else
            ts << " ";
    }
    
    if (values.count() % 3 != 0)
        ts << "\n";
    
    ts << QString("\nobject \"griddata\" class field\n\n");

    f.close();
}

/** Write the passed volume map to the file 'filename' */
void DX::write(const VolumeMap &volmap, const QString &filename, const PropertyMap &map) const
{
    this->write(volmap.gridInfo(), volmap.occupancy(), filename, map);
}
