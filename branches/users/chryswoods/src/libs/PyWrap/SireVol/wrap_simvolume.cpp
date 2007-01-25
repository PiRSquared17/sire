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

/**
  * This file contains the boost::python wrapping of Space and PairMatrix
  *
  * @author Christopher Woods
  *
  */

#include <Python.h>
#include <boost/python.hpp>

#include "SireMol/qhash_siremol.h"

#include "SireVol/space.h"
#include "SireVol/cartesian.h"
#include "SireVol/periodicbox.h"

#include "SireVol/coordgroup.h"

using namespace boost::python;

namespace SireVol
{

void
SIREVOL_EXPORT
export_Space()
{
    class_<SpaceBase, boost::noncopyable>("SpaceBase", no_init)

        .def( "what", (const char* (SpaceBase::*)() const)
                          &SpaceBase::what )

        .def( "calcDist", (double (SpaceBase::*)(const CoordGroup&,
                                             DistMatrix&) const)
                          &SpaceBase::calcDist )

        .def( "calcDist2", (double (SpaceBase::*)(const CoordGroup&,
                                              DistMatrix&) const)
                          &SpaceBase::calcDist2 )
        .def( "calcInvDist", (double (SpaceBase::*)(const CoordGroup&,
                                                DistMatrix&) const)
                          &SpaceBase::calcInvDist )

        .def( "calcInvDist2", (double (SpaceBase::*)(const CoordGroup&,
                                                 DistMatrix&) const)
                          &SpaceBase::calcInvDist2 )

        .def( "calcDist", (double (SpaceBase::*)(const CoordGroup&,
                                             const CoordGroup&,
                                             DistMatrix&) const)
                          &SpaceBase::calcDist )

        .def( "calcDist2", (double (SpaceBase::*)(const CoordGroup&,
                                              const CoordGroup&,
                                              DistMatrix&) const)
                          &SpaceBase::calcDist2 )

        .def( "calcInvDist", (double (SpaceBase::*)(const CoordGroup&,
                                                const CoordGroup&,
                                                DistMatrix&) const)
                          &SpaceBase::calcInvDist )

        .def( "calcInvDist2", (double (SpaceBase::*)(const CoordGroup&,
                                                 const CoordGroup&,
                                                 DistMatrix&) const)
                          &SpaceBase::calcInvDist2 )

        .def( "beyond", (bool (SpaceBase::*)(double,
                                         const CoordGroup&,
                                         const CoordGroup&) const)
                          &SpaceBase::beyond )

        .def( "minimumDistance", (double (SpaceBase::*)(const CoordGroup&,
                                                        const CoordGroup&) const)
                          &SpaceBase::minimumDistance )

        .def( "minimumDistance", (double (SpaceBase::*)(const CoordGroup&) const)
                          &SpaceBase::minimumDistance )

        .def( "moveToCenterBox", (CoordGroup (SpaceBase::*)(const CoordGroup&) const)
                          &SpaceBase::moveToCenterBox )

        .def( "moveToCenterBox", (QVector<CoordGroup> (SpaceBase::*)(const QVector<CoordGroup>&) const)
                          &SpaceBase::moveToCenterBox )

        .def( "getMinimumImage", (CoordGroup (SpaceBase::*)(const CoordGroup&,
                                                            const Vector&) const)
                          &SpaceBase::getMinimumImage )

        .def( "getMinimumImage", (QVector<CoordGroup> (SpaceBase::*)(const QVector<CoordGroup>&,
                                                                     const Vector&) const)
                          &SpaceBase::getMinimumImage )

        .def( "getCopiesWithin", (QList< boost::tuple<double,CoordGroup> >
                                  (SpaceBase::*)(const CoordGroup&, const CoordGroup&,
                                                 double) const)
                          &SpaceBase::getCopiesWithin )
    ;

    class_<Space>( "Space", init<>() )

        .def( init<const Space&>() )
        .def( init<const SpaceBase&>() )

        .def( "what", (const char* (Space::*)() const)
                          &Space::what )

        .def( "calcDist", (double (Space::*)(const CoordGroup&,
                                             DistMatrix&) const)
                          &Space::calcDist )

        .def( "calcDist2", (double (Space::*)(const CoordGroup&,
                                              DistMatrix&) const)
                          &Space::calcDist2 )
        .def( "calcInvDist", (double (Space::*)(const CoordGroup&,
                                                DistMatrix&) const)
                          &Space::calcInvDist )

        .def( "calcInvDist2", (double (Space::*)(const CoordGroup&,
                                                 DistMatrix&) const)
                          &Space::calcInvDist2 )

        .def( "calcDist", (double (Space::*)(const CoordGroup&,
                                             const CoordGroup&,
                                             DistMatrix&) const)
                          &Space::calcDist )

        .def( "calcDist2", (double (Space::*)(const CoordGroup&,
                                              const CoordGroup&,
                                              DistMatrix&) const)
                          &Space::calcDist2 )

        .def( "calcInvDist", (double (Space::*)(const CoordGroup&,
                                                const CoordGroup&,
                                                DistMatrix&) const)
                          &Space::calcInvDist )

        .def( "calcInvDist2", (double (Space::*)(const CoordGroup&,
                                                 const CoordGroup&,
                                                 DistMatrix&) const)
                          &Space::calcInvDist2 )

        .def( "beyond", (bool (Space::*)(double,
                                         const CoordGroup&,
                                         const CoordGroup&) const)
                          &Space::beyond )

        .def( "minimumDistance", (double (Space::*)(const CoordGroup&,
                                                    const CoordGroup&) const)
                          &Space::minimumDistance )

        .def( "minimumDistance", (double (Space::*)(const CoordGroup&) const)
                          &Space::minimumDistance )

        .def( "moveToCenterBox", (CoordGroup (Space::*)(const CoordGroup&) const)
                          &Space::moveToCenterBox )

        .def( "moveToCenterBox", (QVector<CoordGroup> (Space::*)(const QVector<CoordGroup>&) const)
                          &Space::moveToCenterBox )

        .def( "getMinimumImage", (CoordGroup (Space::*)(const CoordGroup&,
                                                        const Vector&) const)
                          &Space::getMinimumImage )

        .def( "getMinimumImage", (QVector<CoordGroup> (Space::*)(const QVector<CoordGroup>&,
                                                                 const Vector&) const)
                          &Space::getMinimumImage )

        .def( "getCopiesWithin", (QList< boost::tuple<double,CoordGroup> >
                                  (Space::*)(const CoordGroup&, const CoordGroup&,
                                             double) const)
                          &Space::getCopiesWithin )
    ;


    class_<Cartesian, bases<SpaceBase> >("Cartesian", init<>())
        .def( init<const Cartesian&>() )
    ;

    class_<PeriodicBox, bases<SpaceBase> >("PeriodicBox", init<>())

        .def(init<const Vector&, const Vector&>())
        .def( init<const PeriodicBox&>() )

        .def( "setDimension", (void (PeriodicBox::*)(const Vector&,
                                                     const Vector&))
                          &PeriodicBox::setDimension )

        .def( "minCoords", (const Vector& (PeriodicBox::*)() const)
                          &PeriodicBox::minCoords,
                          return_value_policy<copy_const_reference>() )

        .def( "maxCoords", (const Vector& (PeriodicBox::*)() const)
                          &PeriodicBox::maxCoords,
                          return_value_policy<copy_const_reference>() )

        .def( "center", (Vector (PeriodicBox::*)() const)
                          &PeriodicBox::center )

        .def( "contains", (bool (PeriodicBox::*)(const Vector&) const)
                          &PeriodicBox::contains )
    ;
}

}
