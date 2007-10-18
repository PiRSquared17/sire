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

#include "localsimsystem.h"
#include "checkpoint.h"

using namespace SireSystem;

/** Construct a simulation that will simulate a copy of a
    system at the checkpoint 'checkpoint' */
LocalSimSystem::LocalSimSystem(const CheckPoint &checkpoint)
               : SimSystem()
{
    this->rollBack(checkpoint);
}

/** Destructor */
LocalSimSystem::~LocalSimSystem()
{}

/** Roll this simulation back to a previous checkpoint */
void LocalSimSystem::rollBack(const CheckPoint &checkpoint)
{
    local_sysdata = checkpoint.info();
    local_ffields = checkpoint.forceFields();
    local_monitors = checkpoint.monitors();

    SimSystem::setSystem(local_sysdata,
                         local_ffields,
                         local_monitors);
}
