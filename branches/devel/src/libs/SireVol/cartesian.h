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

#ifndef SIREVOL_CARTESIAN_H
#define SIREVOL_CARTESIAN_H

#include "space.h"

SIRE_BEGIN_HEADER

namespace SireVol
{
class Cartesian;
}

QDataStream& operator<<(QDataStream&, const SireVol::Cartesian&);
QDataStream& operator>>(QDataStream&, SireVol::Cartesian&);

namespace SireVol
{

/**
This class overloads SimVolume to provide an infinite Cartesian (3-dimensional, orthoganol dimensions) volume. This corresponds to a traditional gas-phase or no-boundary system.

@author Christopher Woods
*/
class SIREVOL_EXPORT Cartesian : public SpaceBase
{

friend QDataStream& ::operator<<(QDataStream&, const Cartesian&);
friend QDataStream& ::operator>>(QDataStream&, Cartesian&);

public:
    Cartesian();
    Cartesian(const Cartesian &other);

    ~Cartesian();

    static const char* typeName()
    {
        return "SireVol::Cartesian";
    }

    const char* what() const
    {
        return Cartesian::typeName();
    }

    SpaceBase* clone() const
    {
        return new Cartesian(*this);
    }

    double calcDist(const CoordGroup &group, DistMatrix &mat) const;
    double calcDist2(const CoordGroup &group, DistMatrix &mat) const;
    double calcInvDist(const CoordGroup &group, DistMatrix &mat) const;
    double calcInvDist2(const CoordGroup &group, DistMatrix &mat) const;

    double calcDist(const CoordGroup &group1, const CoordGroup &group2,
                    DistMatrix &mat) const;

    double calcDist2(const CoordGroup &group1, const CoordGroup &group2,
                     DistMatrix &mat) const;

    double calcInvDist(const CoordGroup &group1, const CoordGroup &group2,
                       DistMatrix &mat) const;

    double calcInvDist2(const CoordGroup &group1, const CoordGroup &group2,
                        DistMatrix &mat) const;

    bool beyond(double dist, const CoordGroup &group0,
                const CoordGroup &group1) const;

    double minimumDistance(const CoordGroup &group0, const CoordGroup &group1) const;

    double minimumDistance(const CoordGroup &group) const;

    CoordGroup moveToCenterBox(const CoordGroup &group) const;

    QVector<CoordGroup> moveToCenterBox(const QVector<CoordGroup> &groups) const;

    CoordGroup getMinimumImage(const CoordGroup &group, const Vector &center) const;

    QVector<CoordGroup> getMinimumImage(const QVector<CoordGroup> &groups,
                                        const Vector &center) const;

    QList< boost::tuple<double,CoordGroup> >
               getCopiesWithin(const CoordGroup &group,
                               const CoordGroup &center, double dist) const;
};

}

Q_DECLARE_METATYPE(SireVol::Cartesian)

SIRE_END_HEADER

#endif
