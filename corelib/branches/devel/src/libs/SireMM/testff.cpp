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

#include "testff.h"

#include "SireVol/cartesian.h"

#include "SireUnits/units.h"

#include <QDebug>
#include <QElapsedTimer>

#include "tostring.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireMaths;
using namespace SireUnits;
using namespace SireVol;

TestFF::TestFF() : cljfunc(new CLJVacShiftAriFunction(15*angstrom, 15*angstrom))
{}

TestFF::TestFF(const TestFF &other)
       : atoms0(other.atoms0), atoms1(other.atoms1), cljfunc(other.cljfunc)
{}

TestFF::~TestFF()
{}

TestFF& TestFF::operator=(const TestFF &other)
{
    atoms0 = other.atoms0;
    atoms1 = other.atoms1;
    cljfunc = other.cljfunc;
    return *this;
}

void TestFF::add(const Molecules &molecules)
{
    atoms0 = CLJAtoms(molecules);
    cljboxes0 = CLJBoxes(atoms0);
}

void TestFF::addFixedAtoms(const Molecules &molecules)
{
    atoms1 = CLJAtoms(molecules);
    cljboxes1 = CLJBoxes(atoms1);
    
    Cartesian space;
    
    QList<CLJBoxDistance> dists = CLJBoxes::getDistances(space, cljboxes1, 15*angstrom);
    
    qDebug() << "N-pairs" << dists.count();
    
    dists = CLJBoxes::getDistances(space, cljboxes0, cljboxes1, 15*angstrom);
    
    qDebug() << "N-pairs" << dists.count();
}

void TestFF::setCutoff(Length coul_cutoff, Length lj_cutoff)
{
    cljfunc.reset( new CLJVacShiftAriFunction(coul_cutoff, lj_cutoff) );
}

void TestFF::calculateEnergy()
{
    QElapsedTimer t;

    double cnrg;
    double ljnrg;

    qDebug() << "inter group energy";

    t.start();

    (*cljfunc)(atoms0, atoms1, cnrg, ljnrg);

    quint64 ns = t.nsecsElapsed();

    qDebug() << "TestFF" << (cnrg+ljnrg) << cnrg << ljnrg << "took" << (0.000001*ns) << "ms";

    qDebug() << "\nintra group energy";
    
    t.start();
    
    (*cljfunc)(atoms1, cnrg, ljnrg);
    
    ns = t.nsecsElapsed();
    
    qDebug() << "TestFF" << (cnrg+ljnrg) << cnrg << ljnrg << "took" << (0.000001*ns) << "ms";
}
