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

#ifndef SIREFF_FFTHREADWORKER_H
#define SIREFF_FFTHREADWORKER_H

#include "ffworker.h"
#include "SireCluster/threadworker.h"

SIRE_BEGIN_HEADER

namespace SireFF
{

class FFCalculatorBase;

/** This class provides a worker that can calculate the energy and
    forces of a ForceField in a background thread (via the passed
    FFCalculator)

    @author Christopher Woods
*/
class SIREFF_EXPORT FFThreadWorker : public FFLocalWorker,
                                     public SireCluster::ThreadWorker
{
public:
    FFThreadWorker(FFCalculatorBase *ffcalculator);

    ~FFThreadWorker();

    ForceField forcefield() const;

protected:
    void _pvt_recalculateEnergy();
    void _pvt_recalculateEnergyFG();

    void _pvt_waitUntilReady();

    double _pvt_getEnergies(Values &nrg_components);

private:
    void calculate();
};

}

SIRE_END_HEADER

#endif
