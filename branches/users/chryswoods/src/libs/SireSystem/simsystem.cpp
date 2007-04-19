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

#include "simsystem.h"
#include "system.h"
#include "moves.h"
#include "move.h"

#include "SireFF/forcefield.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/moleculegroups.h"

#include "SireBase/property.h"

using namespace SireSystem;
using namespace SireMol;
using namespace SireFF;
using namespace SireBase;

/** Protected constructor used by derived classes to provide
    the SimSystem with the components of the System to be
    simulated */
SimSystem::SimSystem(SystemData &sysdata,
                     ForceFieldsBase &ffields,
                     SystemMonitors &monitors)
          : QuerySystem(sysdata, ffields, monitors)
{}

/** Destructor */
SimSystem::~SimSystem()
{}

/** Commit the current configuration of this system - this
    tells the monitors to update themselves using the current
    system configuration. */
void SimSystem::commit()
{
    monitors.update(*this);
    sysdata.incrementMinorVersion();
}

/** Set the property 'name' in all forcefields to the value 'value'

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
void SimSystem::setProperty(const QString &name, const Property &value)
{
    ffields.setProperty(name, value);
    sysdata.incrementMinorVersion();
}

/** Set the property 'name' to the value 'value' in the
    forcefield with ID == ffid

    \throw SireFF::missing_forcefield
    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
void SimSystem::setProperty(ForceFieldID ffid, const QString &name,
                            const Property &value)
{
    ffields.setProperty(ffid, name, value);
    sysdata.incrementMinorVersion();
}

/** Set the property 'name' to the value 'value' in the
    forcefield whose IDs are in 'ffids'

    \throw SireFF::missing_forcefield
    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
void SimSystem::setProperty(const QSet<ForceFieldID> &ffids,
                            const QString &name, const Property &value)
{
    ffields.setProperty(ffids, name, value);
    sysdata.incrementMinorVersion();
}

/** Change the molecule 'molecule' in this system */
void SimSystem::change(const PartialMolecule &molecule)
{
    //change the molecule in a copy of the SystemData
    //(this is to maintain the invariant if something goes
    // wrong when we change the molecule in the forcefields)
    SystemData new_data(sysdata);

    //change the molecule in the data - this also maps it
    //into the system space
    PartialMolecule mapped_mol = new_data.change(molecule);

    //change the molecule in the forcefields
    if (ffields.change(mapped_mol))
        //increment the version of the system
        new_data.incrementMinorVersion();

    //everything was ok - copy the new data to the original
    sysdata = new_data;
}

/** Change lots of molecules in this system */
void SimSystem::change(const QHash<MoleculeID,PartialMolecule> &molecules)
{
    if (molecules.isEmpty())
        return;
    else if (molecules.count() == 1)
    {
        this->change( *(molecules.constBegin()) );
        return;
    }

    //change the molecules in a copy of the SystemData
    //(this is to maintain the invariant if something goes
    // wrong when we change the molecules in the forcefields)
    SystemData new_data(sysdata);

    QHash<MoleculeID,PartialMolecule> mapped_mols = new_data.change(molecules);

    //change the molecules in the forcefields
    if (ffields.change(mapped_mols))
        //increment the version of the system
        new_data.incrementMinorVersion();

    //everything was ok - copy the new data to the original
    sysdata = new_data;
}

/** Change lots of molecules in this system */
void SimSystem::change(const QList<PartialMolecule> &molecules)
{
    if (molecules.isEmpty())
        return;
    else if (molecules.count() == 1)
    {
        this->change( molecules.first() );
        return;
    }

    //change the molecules in a copy of the SystemData
    //(this is to maintain the invariant if something goes
    // wrong when we change the molecules in the forcefields)
    SystemData new_data(sysdata);

    QHash<MoleculeID,PartialMolecule> mapped_mols = new_data.change(molecules);

    //change the molecules in the forcefields
    if (ffields.change(mapped_mols))
        //increment the version of the system
        new_data.incrementMinorVersion();

    //everything was ok - copy the new data to the original
    sysdata = new_data;
}

/** Add the molecule 'molecule' to the molecule groups whose
    IDs are in 'molgroupids'

    \throw SireMol::missing_group
*/
void SimSystem::add(const PartialMolecule &molecule,
                    const QSet<MoleculeGroupID> &molgroupids)
{
#warning - I am here, ready to update to the new new add functions!

    if (not molgroupids.isEmpty())
    {
        //map this molecule into the System space
        PartialMolecule mapped_mol = sysdata.mapIntoSystemSpace(molecule);

        //select only the atoms to be added - the API of MoleculeGroup and
        //forcefield will ensure that the whole molecule is changed if
        //it has been mapped
        mapped_mol.setSelection(molecules.selectedAtoms()); 

        sysdata.add(mapped_mol, molgroupids);
        sysdata.incrementMinorVersion();
    }
}

/** Add the molecule 'molecule' to the molecule groups whose
    IDs are in 'molgroupids' and to the forcefield groups
    whose IDs are in ffgroupids

    \throw SireMol::missing_group
    \throw SireFF::missing_forcefield
    \throw SireFF::invalid_group
*/
void SimSystem::add(const PartialMolecule &molecule,
                    const QSet<FFGroupID> &ffgroupids,
                    const QSet<MoleculeGroupID> &molgroupids)
{
    if (ffgroupids.isEmpty())
        this->add(molecule, molgroupids);
    else
    {
        PartialMolecule mapped_mol = sysdata.mapIntoSystemSpace(molecule);

        //select only the atoms to be added - the API of MoleculeGroup and
        //forcefield will ensure that the whole molecule is changed if
        //it has been mapped
        mapped_mol.setSelection(molecules.selectedAtoms()); 

        SystemData new_data(sysdata);

        new_data.add(mapped_mol, molgroupids);
        ffields.add(mapped_mol, ffgroupids);

        new_data.incrementMinorVersion();
        sysdata = new_data;
    }
}

/** Add lots of molecules to the molecule groups whose
    IDs are in 'molgroupids'

    \throw SireMol::missing_group
*/
void SimSystem::add(const QHash<MoleculeID,PartialMolecule> &molecules,
                    const QSet<MoleculeGroupID> &molgroupids)
{
    if ( not (molgroupids.isEmpty() or molecules.isEmpty()) )
    {
        if (molecules.count() == 1)
        {
            this->add( *(molecules.constBegin()), molgroupids );
            return;
        }
        
        //map this molecule into the System space
        QHash<MoleculeID,PartialMolecule> mapped_mols
                          = sysdata.mapIntoSystemSpace(molecules);

        //select only the atoms to be added - the API of MoleculeGroup and
        //forcefield will ensure that the whole molecule is changed if
        //it has been mapped
        for (QHash<MoleculeID,PartialMolecule>::iterator it = mapped_mols.begin();
             it != mapped_mols.end();
             ++it)
        {
            BOOST_ASSERT( molecules.contains(it.key()) );
        
            it->setSelection( molecules.constFind(it.key())->selectedAtoms() );
        }
        
        sysdata.add(mapped_mols, molgroupids);
        sysdata.incrementMinorVersion();
    }
}

/** Add lots of molecules to the molecule groups whose
    IDs are in 'molgroupids' and to the forcefield groups
    whose IDs are in ffgroupids

    \throw SireMol::missing_group
    \throw SireFF::missing_forcefield
    \throw SireFF::invalid_group
*/
void SimSystem::add(const QHash<MoleculeID,PartialMolecule> &molecules,
                    const QSet<FFGroupID> &ffgroupids,
                    const QSet<MoleculeGroupID> &molgroupids)
{
    if (molecules.isEmpty())
        return;
    else if (ffgroupids.isEmpty())
        this->add(molecule, molgroupids);
    else if (molecules.count() == 1)
        this->add( *(molecules.constBegin()), ffgroupids, molgroupids );
    else
    {
        QHash<MoleculeID,PartialMolecule> mapped_mols
                          = sysdata.mapIntoSystemSpace(molecules);

        SystemData new_data(sysdata);

        //select only the atoms to be added - the API of MoleculeGroup and
        //forcefield will ensure that the whole molecule is changed if
        //it has been mapped
        for (QHash<MoleculeID,PartialMolecule>::iterator it = mapped_mols.begin();
             it != mapped_mols.end();
             ++it)
        {
            BOOST_ASSERT( molecules.contains(it.key()) );
            
            it->setSelection( molecules.constFind(it.key())->selectedAtoms() );
        }
        
        new_data.add( mapped_mols, molgroupids );
        ffields.add( mapped_mols, ffgroupids );

        new_data.incrementMinorVersion();
        sysdata = new_data;
    }
}

/** Add lots of molecules to the molecule groups whose
    IDs are in 'molgroupids'

    \throw SireMol::missing_group
*/
void SimSystem::add(const QList<PartialMolecule> &molecules,
                    const QSet<MoleculeGroupID> &molgroupids)
{
    if ( not (molgroupids.isEmpty() or molecules.isEmpty()) )
    {
        if (molecules.count() == 1)
        {
            this->add(molecules.first(), molgroupids);
        }
        else
        {
            //convert the list into a hash and add that
            QHash<MoleculeID,PartialMolecule> mols;
            mols.reserve(molecules.count());
        
            for (QList<PartialMolecule>::const_iterator it = molecules.begin();
                 it != molecules.end();
                 ++it)
            {
                mols.insert( it->key(), *it );
            }
            
            this->add(mols, molgroupids);
        }
    }
}

/** Add lots of molecules to the molecule groups whose
    IDs are in 'molgroupids' and to the forcefield groups
    whose IDs are in ffgroupids

    \throw SireMol::missing_group
    \throw SireFF::missing_forcefield
    \throw SireFF::invalid_group
*/
void SimSystem::add(const QList<PartialMolecule> &molecules,
                    const QSet<FFGroupID> &ffgroupids,
                    const QSet<MoleculeGroupID> &molgroupids)
{
    if (molecules.isEmpty())
        return;
    else if (ffgroupids.isEmpty())
        this->add(molecule, molgroupids);
    else
    {
        if (molecules.count() == 1)
        {
            this->add(molecules.first(), ffgroupids, molgroupids);
        }
        else
        {
            //convert the list into a hash and add that
            QHash<MoleculeID,PartialMolecule> mols;
            mols.reserve(molecules.count());
        
            for (QList<PartialMolecule>::const_iterator it = molecules.begin();
                 it != molecules.end();
                 ++it)
            {
                mols.insert( it->key(), *it );
            }
            
            this->add(mols, ffgroupids, molgroupids);
        }
    }
}

/** Completely remove all selected atoms of any version
    of the molecule 'molecule' from the system */
void SimSystem::remove(const PartialMolecule &molecule)
{
    SystemData new_data( sysdata );

    new_data.remove(molecule);
    ffields.remove(molecule);

    new_data.incrementMinorVersion();
    sysdata = newdata;
}

/** Completely remove all selected atoms of any version
    of the molecules in 'molecules' */
void SimSystem::remove(const QList<PartialMolecule> &molecules)
{
    if (not molecules.isEmpty())
    {
        SystemData new_data( sysdata );

        new_data.remove(molecules);
        ffields.remove(molecules);

        new_data.incrementMinorVersion();
        sysdata = newdata;
    }
}

/** Completely remove all selected atoms from any version
    of the molecule 'molecule' from all of the molecule groups
    whose IDs are in 'molgroupids'

    \throw SireMol::missing_group
*/
void SimSystem::remove(const PartialMolecule &molecule,
                       const QSet<MoleculeGroupID> &molgroupids)
{
    if (not molgroupids.isEmpty())
    {
        new_data.remove(molecule, molgroupids);
        new_data.incrementMinorVersion();
    }
}

/** Completely remove all selected atoms from any version
    of the molecule 'molecule' from all of the molecule groups
    and forcefields whose IDs are in 'ffgroupids' and
    'molgroupids'

    \throw SireFF::missing_forcefield
    \throw SireFF::invalid_group
    \throw SireMol::missing_group
*/
void SimSystem::remove(const PartialMolecule &molecule,
                       const QSet<FFGroupID> &ffgroupids,
                       const QSet<MoleculeGroupID> &molgroupids)
{
    if (ffgroupids.isEmpty())
        this->remove(molecule, molgroupids);
    else
    {
        SystemData new_data(sysdata);
        new_data.remove(molecule, molgroupids);
        ffields.remove(molecule, ffgroupids);

        new_data.incrementMinorVersion();
        sysdata = new_data;
    }
}

/** Completely remove all selected atoms from any version
    of the molecules in 'molecules' from all of the molecule groups
    whose IDs are in 'molgroupids'

    \throw SireMol::missing_group
*/
void SimSystem::remove(const QList<PartialMolecule> &molecules,
                       const QSet<MoleculeGroupID> &molgroupids)
{
    if ( not (molecules.isEmpty() or molgroupids.isEmpty()) )
    {
        sysdata.remove(molecules, molgroupids);
        sysdata.incrementMinorVersion();
    }
}

/** Completely remove all selected atoms from any version
    of the molecules in 'molecules' from all of the molecule groups
    and forcefields whose IDs are in 'ffgroupids' and
    'molgroupids'

    \throw SireFF::missing_forcefield
    \throw SireFF::invalid_group
    \throw SireMol::missing_group
*/
void SimSystem::remove(const QList<PartialMolecule> &molecules,
                       const QSet<FFGroupID> &ffgroupids,
                       const QSet<MoleculeGroupID> &molgroupids)
{
    if (molecules.isEmpty())
        return;
    else if (ffgroupids.isEmpty())
        this->remove(molecules,molgroupids);
    else
    {
        SystemData new_data(sysdata);

        new_data.remove(molecules, molgroupids);
        ffields.remove(molecules, ffgroupids);

        new_data.incrementMinorVersion();
        sysdata = new_data;
    }
}
