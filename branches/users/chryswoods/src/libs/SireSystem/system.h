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

#ifndef SIRESYSTEM_SYSTEM_H
#define SIRESYSTEM_SYSTEM_H

#include <QString>

#include "systemdata.h"
#include "systemmonitors.h"

#include "SireFF/forcefields.h"

namespace SireSystem
{
class System;
}

QDataStream& operator<<(QDataStream&, const SireSystem::System&);
QDataStream& operator>>(QDataStream&, SireSystem::System&);

namespace SireSystem
{

class Move;
class Moves;

class CheckPoint;

using SireCAS::Function;

using SireFF::ForceField;
using SireFF::ForceFields;

using SireFF::FFExpression;
using SireFF::FFComponent;

/** This class holds all of the data necessary to specify a single
    system of molecules (including how to calculate their energy).

    This is the class/interface visible to the user for creating,
    editing and using simulation systems. This class is also used
    as the checkpoint object for saving the complete state
    of a simulation.

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT System
{

friend QDataStream& ::operator<<(QDataStream&, const System&);
friend QDataStream& ::operator>>(QDataStream&, System&);

friend class SimSystem;  //friend so can mess with this System ;-)

public:
    System(const MoleculeGroups &groups = MoleculeGroups(),
           const ForceFields &forcefields = ForceFields(),
           const SystemMonitors &monitors = SystemMonitors());

    System(const QString &name,
           const MoleculeGroups &groups = MoleculeGroups(),
           const ForceFields &forcefields = ForceFields(),
           const SystemMonitors &monitors = SystemMonitors());

    System(const MoleculeGroups &groups,
           const SystemMonitors &monitors);

    System(const QString &name,
           const MoleculeGroups &groups,
           const SystemMonitors &monitors);

    System(const ForceFields &forcefields,
           const SystemMonitors &monitors = SystemMonitors());

    System(const QString &name,
           const ForceFields &forcefields,
           const SystemMonitors &monitors = SystemMonitors());

    System(const SystemMonitors &monitors);
    System(const QString &name, const SystemMonitors &monitors);

    System(const CheckPoint &checkpoint);

    System(QuerySystem &simsystem);

    System(const System &other);

    ~System();

    System& operator=(const System &other);
    System& operator=(const CheckPoint &checkpoint);
    System& operator=(QuerySystem &simsystem);

    bool operator==(const System &other) const;
    bool operator!=(const System &other) const;

    const ForceFields& forceFields() const;
    const SystemData& info() const;
    const SystemMonitors& monitors() const;

    void setSpace(const Space &space);

    Moves run(const Move &move, quint32 nmoves=1);
    Moves run(const Moves &moves);
    Moves run(const Moves &moves, quint32 nmoves);

    SystemID ID() const;
    const Version& version() const;

    void prepareForSimulation();

    CheckPoint checkPoint();

private:
    System(const SystemData &sysdata,
           const ForceFields &ffields);

    /** The object holding information about this system */
    SystemData sysdata;

    /** The forcefields that are used to calculate
        energies / forces of molecules */
    ForceFields ffields;

    /** The monitors that are used to monitor properties of
        this system */
    SystemMonitors sysmonitors;

    /** Whether or not the system is in a consistent state */
    bool is_consistent;
};

/** Return the ID number of the system */
inline SystemID System::ID() const
{
    return sysdata.ID();
}

/** Return the version number of the system */
inline const Version& System::version() const
{
    return sysdata.version();
}

/** Return the forcefields contained in this system */
inline const ForceFields& System::forceFields() const
{
    return ffields;
}

/** Return the information about this system */
inline const SystemData& System::info() const
{
    return sysdata;
}

/** Return the monitors that are used to monitor properties
    of this system */
inline const SystemMonitors& System::monitors() const
{
    return sysmonitors;
}

}

Q_DECLARE_METATYPE(SireSystem::System);

#endif
