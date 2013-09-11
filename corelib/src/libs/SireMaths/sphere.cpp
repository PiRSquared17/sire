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

#include "sphere.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireMaths;

static const RegisterMetaType<Sphere> r_sphere;

/** Serialise to a binary data stream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const Sphere &sphere)
{
    writeHeader(ds, r_sphere, 1) << sphere._center << sphere._radius;

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, Sphere &sphere)
{
    VersionID v = readHeader(ds, r_sphere);

    if (v == 1)
    {
        ds >> sphere._center >> sphere._radius;
    }
    else
        throw version_error(v, "1", r_sphere, CODELOC);

    return ds;
}

/** Create a default sphere, centered on the origin with zero radius */
Sphere::Sphere() : _center(0.0), _radius(0.0)
{}

/** Create a sphere with radius 'radius', centered on the origin */
Sphere::Sphere(const double &radius) : _center(0.0), _radius(radius)
{
    if (radius < 0.0)
        _radius = -radius;
}

/** Create a sphere centered at 'position' and with radius 'radius' */
Sphere::Sphere(const Vector &position, const double &radius)
       : _center(position), _radius(radius)
{
    if (radius < 0.0)
        _radius = -radius;
}

/** Copy constructor */
Sphere::Sphere(const Sphere &other)
       : _center(other._center), _radius(other._radius)
{}

/** Destructor */
Sphere::~Sphere()
{}

/** Comparison operator */
bool Sphere::operator==(const Sphere &other) const
{
    return _radius == other._radius and _center == other._center;
}

/** Comparison operator */
bool Sphere::operator!=(const Sphere &other) const
{
    return _radius != other._radius or _center != other._center;
}

/** Set the position of the center of this sphere */
void Sphere::setPosition(const Vector &position)
{
    _center = position;
}

/** Set the position of the center of this sphere */
void Sphere::setCenter(const Vector &center)
{
    _center = center;
}

/** Set the radius of this sphere */
void Sphere::setRadius(double radius)
{
    _radius = radius;

    if (radius < 0.0)
        _radius = 0.0;
}

const char* Sphere::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Sphere>() );
}
