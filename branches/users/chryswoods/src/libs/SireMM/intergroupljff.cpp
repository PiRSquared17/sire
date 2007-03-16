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

#include "SireCAS/qhash_sirecas.h"
#include "SireMol/qhash_siremol.h"

#include "intergroupljff.h"

#include "SireMol/errors.h"
#include "SireFF/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireFF;
using namespace SireStream;

///////////
/////////// Implementation of InterGroupLJFF::Components
///////////

/** Constructor */
InterGroupLJFF::Components::Components() : LJFF::Components()
{}

/** Constructor that just passes its arguments up to the parent's constructor */
InterGroupLJFF::Components::Components(const FFBase &ffbase, const Symbols &symbols)
               : LJFF::Components(ffbase,symbols)
{}

/** Copy constructor */
InterGroupLJFF::Components::Components(const Components &other)
               : LJFF::Components(other)
{}

/** Destructor */
InterGroupLJFF::Components::~Components()
{}

///////////
/////////// Implementation of InterGroupLJFF::Parameters
///////////

/** Constructor */
InterGroupLJFF::Parameters::Parameters() : LJFF::Parameters()
{}

/** Copy constructor */
InterGroupLJFF::Parameters::Parameters(const Parameters &other)
               : LJFF::Parameters(other)
{}

/** Destructor */
InterGroupLJFF::Parameters::~Parameters()
{}

///////////
/////////// Implementation of InterGroupLJFF::Groups
///////////

/** Constructor */
InterGroupLJFF::Groups::Groups() : LJFF::Groups()
{
    a = this->getUniqueID();
    b = this->getUniqueID();
}

/** Copy constructor */
InterGroupLJFF::Groups::Groups(const Groups &other)
               : LJFF::Groups(other), a(other.a), b(other.b)
{}

/** Destructor */
InterGroupLJFF::Groups::~Groups()
{}

/** Static instance of this class returned by all InterGroupLJFF objects */
InterGroupLJFF::Groups InterGroupLJFF::Groups::default_group;

///////////
/////////// Implementation of InterGroupLJFF
///////////

static const RegisterMetaType<InterGroupLJFF> r_interljff;

/** Serialise to a binary data stream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const InterGroupLJFF &interljff)
{
    writeHeader(ds, r_interljff, 1);

    SharedDataStream sds(ds);

    sds << interljff.mols[0]
        << interljff.mols[1]
        << interljff.changed_mols[0]
        << interljff.changed_mols[1]
        << static_cast<const LJFF&>(interljff);

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      InterGroupLJFF &interljff)
{
    VersionID v = readHeader(ds, r_interljff);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> interljff.mols[0]
            >> interljff.mols[1]
            >> interljff.changed_mols[0]
            >> interljff.changed_mols[1]
            >> static_cast<LJFF&>(interljff);

        //rebuild mols and changed_mols

        for (int group=0; group<2; ++group)
        {
            uint nmols = interljff.mols[group].count();

            if (nmols > 0)
            {
                QHash<MoleculeID,uint> molid_to_index;
                molid_to_index.reserve(nmols);

                const LJFF::LJMolecule *mol_array = interljff.mols[group].constData();

                for (uint i=0; i<nmols; ++i)
                {
                    molid_to_index.insert( mol_array[i].molecule().ID(), i );
                }

                interljff.molid_to_index[group] = molid_to_index;
            }
            else
                interljff.molid_to_index[group].clear();

            nmols = interljff.changed_mols[group].count();

            if (nmols > 0)
            {
                QHash<MoleculeID,uint> molid_to_changedindex;
                molid_to_changedindex.reserve(nmols);

                const LJFF::ChangedLJMolecule *mol_array
                                = interljff.changed_mols[group].constData();

                for (uint i=0; i<nmols; ++i)
                {
                    MoleculeID molid = mol_array[i].oldMolecule().molecule().ID();

                    if (molid == 0)
                        molid = mol_array[i].newMolecule().molecule().ID();

                    BOOST_ASSERT(molid != 0);

                    molid_to_changedindex.insert( molid, i );
                }

                interljff.molid_to_changedindex[group] = molid_to_changedindex;
            }
            else
                interljff.molid_to_changedindex[group].clear();
        }
    }
    else
        throw version_error(v, "1", r_interljff, CODELOC);

    return ds;
}

/** Constructor */
InterGroupLJFF::InterGroupLJFF() : LJFF()
{}

/** Construct a LJ forcefield using the passed Space, combining rules and
    switching function */
InterGroupLJFF::InterGroupLJFF(const Space &space,
                               const SwitchingFunction &switchfunc)
               : LJFF(space, switchfunc)
{}

/** Copy constructor */
InterGroupLJFF::InterGroupLJFF(const InterGroupLJFF &other)
               : LJFF(other)
{
    for (int i=0; i<2; ++i)
    {
        mols[i] = other.mols[i];
        molid_to_index[i] = other.molid_to_index[i];

        changed_mols[i] = other.changed_mols[i];
        molid_to_changedindex[i] = other.molid_to_changedindex[i];
    }
}

/** Destructor */
InterGroupLJFF::~InterGroupLJFF()
{}

/** Assignment operator */
InterGroupLJFF& InterGroupLJFF::operator=(const InterGroupLJFF &other)
{
    for (int i=0; i<2; ++i)
    {
        mols[i] = other.mols[i];
        molid_to_index[i] = other.molid_to_index[i];

        changed_mols[i] = other.changed_mols[i];
        molid_to_changedindex[i] = other.molid_to_changedindex[i];
    }

    return *this;
}

/** Return the ID of the other group

    \throw SireFF::invalid_group
*/
int InterGroupLJFF::otherGroup(int group_id) const
{
    BOOST_ASSERT( group_id == 0 or group_id == 1 );
    return (group_id == 0);
}

/** Assert that there is no molecule with the ID == molid in
    any other group than the group with ID == group_id,
    and that there is indeed a group available with this ID number.

    \throw SireFF::invalid_group
*/
void InterGroupLJFF::assertValidGroup(int group_id, MoleculeID molid) const
{
    if ( molid_to_index[ otherGroup(group_id) ].contains(molid) )
    {
        QString newgroup = "A";
        QString oldgroup = "B";

        if ( FFBase::Group(group_id) == groups().B() )
        {
            newgroup = "B";
            oldgroup = "A";
        }

        throw SireFF::invalid_group( QObject::tr(
              "The molecule with ID == %1 cannot be put in the group %2 "
              "as group %3 already contains a molecule with this ID number.")
                  .arg(molid).arg(newgroup,oldgroup), CODELOC );
    }
}

/** Recalculate the total energy of this forcefield from scratch */
void InterGroupLJFF::recalculateTotalEnergy()
{
    //calculate the total LJ energy between pairs of
    //molecules between the two groups

    double ljnrg = 0;

    int nmols0 = mols[0].count();
    int nmols1 = mols[1].count();

    if (nmols0 > 0 and nmols1 > 0)
    {
        const LJMolecule *mols0_array = mols[0].constData();
        const LJMolecule *mols1_array = mols[1].constData();

        for (int i=0; i<nmols0; ++i)
        {
            const LJMolecule &mol0 = mols0_array[i];

            for (int j=0; j<nmols1; ++j)
            {
                const LJMolecule &mol1 = mols1_array[j];

                ljnrg += LJFF::calculateEnergy( mol0, mol1,
                                                space(),
                                                switchingFunction(),
                                                distanceMatrix(),
                                                ljMatrix() );
            }
        }
    }

    this->setComponent( components().lj(), ljnrg );
    this->setComponent( components().total(), ljnrg );

    //clear the list of changed molecules
    for (int i=0; i<2; ++i)
    {
        changed_mols[i].clear();
        molid_to_changedindex[i].clear();
    }
}

/** Recalculate the energy knowing that only the group at index 'changed_idx'
    has changed */
double InterGroupLJFF::recalculateWithOneChangedGroup(int changed_idx)
{
    BOOST_ASSERT(changed_idx == 0 or changed_idx == 1);

    //just calculate the change in energy of the changed parts
    //of the changed group with all of the unchanged group

    //get the index of the unchanged group
    //(0 if 1 has changed, 1 if 0 has changed)
    int unchanged_idx = (changed_idx == 0);

    BOOST_ASSERT(unchanged_idx == 0 or unchanged_idx == 1);

    int nchanged = changed_mols[changed_idx].count();
    const ChangedLJMolecule *changed_array = changed_mols[changed_idx].constData();

    int nmols = mols[unchanged_idx].count();
    const LJMolecule *mols_array = mols[unchanged_idx].constData();

    double ljnrg = 0;

    for (int i=0; i<nchanged; ++i)
    {
        const ChangedLJMolecule &changedmol = changed_array[i];

        for (int j=0; j<nmols; ++j)
        {
            const LJMolecule &mol = mols_array[j];

            //calculate the change in energy between this molecule and the changed
            //parts of the changed molecule
            ljnrg += calculateEnergy(mol, changedmol.newParts(),
                                     space(), switchingFunction(),
                                     distanceMatrix(), ljMatrix());

            ljnrg -= calculateEnergy(mol, changedmol.oldParts(),
                                     space(), switchingFunction(),
                                     distanceMatrix(), ljMatrix());
        }
    }

    return ljnrg;
}

/** Recalculate the energy of the unchanged parts of group 'unchanged_idx'
    with the changed parts of group 'changed_idx' */
double InterGroupLJFF::recalculateChangedWithUnchanged(int unchanged_idx,
                                                       int changed_idx)
{
    BOOST_ASSERT( ( unchanged_idx == 0 or unchanged_idx == 1 ) and
                  ( changed_idx == 0 or changed_idx == 1) and
                  ( unchanged_idx != changed_idx ) );

    int nchanged = changed_mols[changed_idx].count();
    const ChangedLJMolecule *changed_array = changed_mols[changed_idx].constData();

    int nmols = mols[unchanged_idx].count();
    const LJMolecule *mols_array = mols[unchanged_idx].constData();

    double ljnrg = 0;

    for (int i=0; i<nmols; ++i)
    {
        const LJMolecule &mol = mols_array[i];

        //has this molecule changed? if so, then skip it
        if ( molid_to_changedindex[unchanged_idx].contains(mol.molecule().ID()) )
        {
            continue;
        }

        //this molecule hasn't changed, so calculate its energy with all of
        //the changed molecules in the other group
        for (int j=0; j<nchanged; ++j)
        {
            const ChangedLJMolecule &changedmol = changed_array[j];

            ljnrg += calculateEnergy( mol, changedmol.newParts(),
                                      space(), switchingFunction(),
                                      distanceMatrix(), ljMatrix() );

            ljnrg -= calculateEnergy( mol, changedmol.oldParts(),
                                      space(), switchingFunction(),
                                      distanceMatrix(), ljMatrix() );
        }
    }

    return ljnrg;
}

/** Recalculate the energy of the changed parts of one group with the
    changed parts of the other group */
double InterGroupLJFF::recalculateChangedWithChanged()
{
    int nchanged0 = changed_mols[0].count();
    int nchanged1 = changed_mols[1].count();

    const ChangedLJMolecule *changed_array0 = changed_mols[0].constData();
    const ChangedLJMolecule *changed_array1 = changed_mols[1].constData();

    double ljnrg = 0;

    for (int i=0; i<nchanged0; ++i)
    {
        const ChangedLJMolecule &changedmol0 = changed_array0[i];

        if (changedmol0.changedAll())
        {
            //the whole of mol0 has changed - calculate the change in energy
            //of the whole of this molecule with the whole of the changed molecules
            //in group 1
            for (int j=0; j<nchanged1; ++j)
            {
                const ChangedLJMolecule &changedmol1 = changed_array1[j];

                ljnrg += calculateEnergy( changedmol0.newMolecule(),
                                          changedmol1.newMolecule(),
                                          space(), switchingFunction(),
                                          distanceMatrix(), ljMatrix() );

                ljnrg -= calculateEnergy( changedmol0.oldMolecule(),
                                          changedmol1.oldMolecule(),
                                          space(), switchingFunction(),
                                          distanceMatrix(), ljMatrix() );
            }
        }
        else
        {
            //only part of mol0 has changed - what we do now depends on how much
            //the molecules in the other group have changed...
            for (int j=0; j<nchanged1; ++j)
            {
                const ChangedLJMolecule &changedmol1 = changed_array1[j];

                if (changedmol1.changedAll())
                {
                    //all of mol1 has changed - calculate the change in energy
                    //of the whole of this molecule with the whole of mol0
                    ljnrg += calculateEnergy( changedmol0.newMolecule(),
                                              changedmol1.newMolecule(),
                                              space(), switchingFunction(),
                                              distanceMatrix(), ljMatrix() );

                    ljnrg -= calculateEnergy( changedmol0.oldMolecule(),
                                              changedmol1.oldMolecule(),
                                              space(), switchingFunction(),
                                              distanceMatrix(), ljMatrix() );
                }
                else
                {
                    //this is the complicated, but hopefully much rarer edge case
                    // - part of this molecule *and* part of the other molecule have
                    // changed. For now, I will do the less efficient method of calculating
                    // the change in energy, namely calculating the energy of the part
                    // of mol0 that changed with all of mol1, then the energy of the part
                    // of mol1 that changed with all of mol0, then finally, subtracting the
                    // energy of the part of mol0 with the part of mol1 (which would
                    // otherwise be double-counted).

                    //first, the changed parts of mol0 with all of mol1
                    ljnrg += calculateEnergy( changedmol0.newParts(),
                                              changedmol1.newMolecule(),
                                              space(), switchingFunction(),
                                              distanceMatrix(), ljMatrix() );

                    ljnrg -= calculateEnergy( changedmol0.oldParts(),
                                              changedmol1.oldMolecule(),
                                              space(), switchingFunction(),
                                              distanceMatrix(), ljMatrix() );

                    //now the changed parts of mol1 with the whole of mol0
                    ljnrg += calculateEnergy( changedmol0.newMolecule(),
                                              changedmol1.newParts(),
                                              space(), switchingFunction(),
                                              distanceMatrix(), ljMatrix() );

                    ljnrg -= calculateEnergy( changedmol0.oldMolecule(),
                                              changedmol1.oldParts(),
                                              space(), switchingFunction(),
                                              distanceMatrix(), ljMatrix() );

                    //finally, remove the double-counted interaction of the parts
                    //of mol0 that changed with the parts of mol1 that changed
                    ljnrg -= calculateEnergy( changedmol0.newParts(),
                                              changedmol1.newParts(),
                                              space(), switchingFunction(),
                                              distanceMatrix(), ljMatrix() );

                    ljnrg += calculateEnergy( changedmol0.oldParts(),
                                              changedmol1.oldParts(),
                                              space(), switchingFunction(),
                                              distanceMatrix(), ljMatrix() );
                }
            }
        }
    }

    return ljnrg;
}

/** Recalculate the energy knowing that both of the groups have changes */
double InterGroupLJFF::recalculateWithTwoChangedGroups()
{
    //first calculate the energy of the changed parts
    //of each group with the unchanged parts of the other group
    double ljnrg = recalculateChangedWithUnchanged(0,1);
    ljnrg += recalculateChangedWithUnchanged(1,0);

    //now calculate the change in energy between the two changed
    //parts of each group
    ljnrg += recalculateChangedWithChanged();

    return ljnrg;
}

/** Recalculate the energy by using a delta from the old configuration */
void InterGroupLJFF::recalculateViaDelta()
{
    int nmols0 = mols[0].count();
    int nmols1 = mols[1].count();

    int nchanged0 = changed_mols[0].count();
    int nchanged1 = changed_mols[1].count();

    if ( (nchanged0 == 0 and nchanged1 == 0) or nmols0 == 0 or nmols1 == 0 )
    {
        //no changed molecules, or no molecules,
        //so just recalculate the total energy
        this->recalculateTotalEnergy();
        return;
    }

    double ljnrg = 0;

    //calculate the change in energy of group 0 with group 1
    if ( nchanged0 == 0 )
    {
        ljnrg = this->recalculateWithOneChangedGroup(1);
    }
    else if ( nchanged1 == 0 )
    {
        ljnrg = this->recalculateWithOneChangedGroup(0);
    }
    else
    {
        ljnrg = this->recalculateWithTwoChangedGroups();
    }

    this->setComponent( components().lj(), ljnrg );
    this->setComponent( components().total(), ljnrg );

    //clear the list of changed molecules
    for (int i=0; i<2; ++i)
    {
        changed_mols[i].clear();
        molid_to_changedindex[i].clear();
    }
}

/** Calculate the total energy of this forcefield. This will either
    calculate the energy from scratch, or it will calculate it as a
    delta based on the parts of the forcefield that have moved since
    the last evaluation. */
void InterGroupLJFF::recalculateEnergy()
{
    if (changed_mols[0].isEmpty() and changed_mols[1].isEmpty())
    {
        this->recalculateTotalEnergy();
    }
    else
    {
        this->recalculateViaDelta();
    }
}

/** Update the current state of 'molecule' in the forcefield */
void InterGroupLJFF::updateCurrentState(int group_idx,
                                        const LJFF::LJMolecule &molecule)
{
    MoleculeID molid = molecule.molecule().ID();

    //get the group that contains this molecule (one of the groups *must*
    //contain this molecule!)
    BOOST_ASSERT( group_idx == 0 or group_idx == 1 );
    BOOST_ASSERT( molid_to_index[group_idx].contains(molid) );
    BOOST_ASSERT( not molid_to_index[ group_idx == 0 ].contains(molid) );

    uint idx = molid_to_index[group_idx].value(molid);

    mols[group_idx][idx] = molecule;
}

/** Add the 'molecule' to the current state of group 'group_idx' */
void InterGroupLJFF::addToCurrentState(int group_idx,
                                       const LJFF::LJMolecule &molecule)
{
    BOOST_ASSERT( group_idx == 0 or group_idx == 1 );

    MoleculeID molid = molecule.molecule().ID();

    BOOST_ASSERT( not (molid_to_index[0].contains(molid) or
                       molid_to_index[1].contains(molid)) );

    uint idx = mols[group_idx].count();
    molid_to_index[group_idx].insert( molid, idx );
    mols[group_idx].append(molecule);
}

/** Remove the molecule with ID == molid from the current state
    of the forcefield */
void InterGroupLJFF::removeFromCurrentState(int group_idx, MoleculeID molid)
{
    BOOST_ASSERT( group_idx == 0 or group_idx == 1 );

    if (molid_to_index[group_idx].contains(molid))
    {
        // need to remove from list, reindex hash, remove from hash
        int idx = molid_to_index[group_idx][molid];

        // remove this item from the array of molecules...
        mols[group_idx].remove(idx);

        // reindex the vector
        molid_to_index[group_idx].clear();

        uint nmols = mols[group_idx].count();

        if (nmols > 0)
        {
            molid_to_index[group_idx].reserve(nmols);

            const LJMolecule *mols_array = mols[group_idx].constData();

            for (uint i=0; i<nmols; ++i)
            {
                molid_to_index[group_idx].insert(
                                    mols_array[i].molecule().ID(), i );
            }
        }

        return;
    }
}

/** Return the change record for the molecule with ID == molid. This
    returns an empty record if this molecule is not in this forcefield */
LJFF::ChangedLJMolecule InterGroupLJFF::changeRecord(MoleculeID molid) const
{
    for (int i=0; i<2; ++i)
    {
        //only return the changerecord if the molecule is currently in the
        //group - beware of old changerecords of molecules that were deleted
        //from a group - these will hang around until the next energy
        //evaluation
        if (molid_to_index[i].contains(molid))
        {
            if (molid_to_changedindex[i].contains(molid))
                return changed_mols[i].at( molid_to_changedindex[i][molid] );
            else
                return ChangedLJMolecule( mols[i].at(molid_to_index[i][molid]) );
        }
    }

    //the molecule is not in this forcefield
    return ChangedLJMolecule();
}

/** Record the change described by 'changed_mol' */
void InterGroupLJFF::recordChange(int group_idx, MoleculeID molid,
                                  const ChangedLJMolecule &changed_mol)
{
    //save the change
    QHash<MoleculeID,uint>::const_iterator it
                        = molid_to_changedindex[group_idx].constFind(molid);

    if (it != molid_to_changedindex[group_idx].constEnd())
    {
        //the molecule has been changed before
        changed_mols[group_idx][it.value()] = changed_mol;
    }
    else
    {
        //this is a new change
        uint idx = changed_mols[group_idx].count();

        changed_mols[group_idx].append(changed_mol);
        molid_to_changedindex[group_idx].insert(molid, idx);
    }
}

/** Apply the change described by 'changed_mol' */
bool InterGroupLJFF::applyChange(MoleculeID molid,
                                 const ChangedLJMolecule &changed_mol)
{
    if (changed_mol.isEmpty())
        return false;

    //get the current group in which the molecule resides
    int group_idx = molid_to_index[1].contains(molid);

    if (not molid_to_index[0].contains(molid))
    {
        //this function cannot be used to add molecules!
        throw SireMol::missing_molecule( QObject::tr(
                "You cannot add parts to a molecule that doesn't yet exist "
                "in this forcefield without first saying to which group it "
                "should belong!"), CODELOC );
    }

    BOOST_ASSERT( not changed_mol.oldMolecule().isEmpty() );
    BOOST_ASSERT( group_idx == 0 or group_idx == 1 );

    if (changed_mol.newMolecule().isEmpty())
    {
        //the molecule has been completely removed
        BOOST_ASSERT( molid == changed_mol.oldMolecule().molecule().ID() );
        this->removeFromCurrentState(group_idx, molid);
    }
    else
    {
        BOOST_ASSERT( molid == changed_mol.newMolecule().molecule().ID() );
        this->updateCurrentState(group_idx, changed_mol.newMolecule());
    }

    this->recordChange(group_idx, molid, changed_mol);

    return true;
}

/** Private class used by the "change" functions to actually change the
    molecule or part of molecule */
template<class T>
bool InterGroupLJFF::_pvt_change(const T &mol)
{
    MoleculeID molid = mol.ID();

    ChangedLJMolecule new_molecule = this->changeRecord(molid).change(mol);

    if (this->applyChange(molid, new_molecule))
    {
        this->incrementMinorVersion();
        return true;
    }
    else
        return isDirty();
}

/** Change the molecule 'molecule' */
bool InterGroupLJFF::change(const Molecule &molecule)
{
    return this->_pvt_change<Molecule>(molecule);
}

/** Change the residue 'residue' */
bool InterGroupLJFF::change(const Residue &residue)
{
    return this->_pvt_change<Residue>(residue);
}

/** Change the atom 'atom' */
bool InterGroupLJFF::change(const NewAtom &atom)
{
    return this->_pvt_change<NewAtom>(atom);
}

/** Private class used by the "remove" functions to actually remove the
    molecule or part of molecule */
template<class T>
bool InterGroupLJFF::_pvt_remove(const T &mol)
{
    MoleculeID molid = mol.ID();

    ChangedLJMolecule new_molecule = changeRecord(molid).remove(mol);

    if (this->applyChange(molid, new_molecule))
    {
        this->incrementMajorVersion();
        return true;
    }
    else
        return isDirty();
}

/** Remove the molecule 'molecule' */
bool InterGroupLJFF::remove(const Molecule &molecule)
{
    return this->_pvt_remove<Molecule>(molecule);
}

/** Remove the residue 'residue' */
bool InterGroupLJFF::remove(const Residue &residue)
{
    return this->_pvt_remove<Residue>(residue);
}

/** Remove the atom 'atom' */
bool InterGroupLJFF::remove(const NewAtom &atom)
{
    return this->_pvt_remove<NewAtom>(atom);
}

/** Remove the selected atoms from 'molecule' from this forcefield */
bool InterGroupLJFF::remove(const Molecule &molecule, const AtomSelection &selected_atoms)
{
    MoleculeID molid = molecule.ID();

    ChangedLJMolecule new_molecule = changeRecord(molid).remove(selected_atoms);

    if (this->applyChange(molid, new_molecule))
    {
        this->incrementMajorVersion();
        return true;
    }
    else
        return isDirty();
}

/** Private class used by the "add" functions to actually add the
    molecule or part of molecule */
template<class T>
bool InterGroupLJFF::_pvt_add(const T &mol, const ParameterMap &map)
{
    //get the molecule's ID
    MoleculeID molid = mol.ID();

    ChangedLJMolecule new_molecule =
                changeRecord(molid).add( mol, map.source(this->parameters().lj()) );

    if (this->applyChange(molid, new_molecule))
    {
        this->incrementMajorVersion();
        return true;
    }
    else
        return isDirty();
}

/** Add the molecule 'molecule' to this forcefield using the optionally
    supplied map to find the forcefield parameters amongst the
    molecule's properties.

    \throw SireMol::missing_property
    \throw SireError::invalid_cast
*/
bool InterGroupLJFF::add(const Molecule &molecule, const ParameterMap &map)
{
    return this->_pvt_add<Molecule>(molecule, map);
}

/** Add the residue 'residue' to this forcefield using the optionally
    supplied map to find the forcefield parameters amongst the
    residue's properties. Note that the property used must agree
    with the rest of the molecule if it is already in this forcefield.

    \throw SireMol::missing_property
    \throw SireMol::invalid_cast
*/
bool InterGroupLJFF::add(const Residue &residue, const ParameterMap &map)
{
    return this->_pvt_add<Residue>(residue, map);
}

/** Add the atom 'atom' to this forcefield using the optionally
    supplied map to find the forcefield parameters amongst the
    atom's properties. Note that the property used must agree
    with the rest of the molecule if it is already in this forcefield.

    \throw SireMol::missing_property
    \throw SireMol::invalid_cast
*/
bool InterGroupLJFF::add(const NewAtom &atom, const ParameterMap &map)
{
    return this->_pvt_add<NewAtom>(atom, map);
}

/** Add the selected atoms from 'molecule' to this forcefield
    using the supplied map to find the forcefield parameters
    amongst the molecule's properties. Note that the property
    used must agree with the rest of the molecule if it is already
    in this forcefield.

    \throw SireMol::missing_property
    \throw SireMol::invalid_cast
*/
bool InterGroupLJFF::add(const Molecule &molecule, const AtomSelection &selected_atoms,
                         const ParameterMap &map)
{
    //get the molecule's ID
    MoleculeID molid = molecule.ID();

    ChangedLJMolecule new_molecule =
                changeRecord(molid).add( selected_atoms,
                                         map.source(this->parameters().lj()) );

    if (this->applyChange(molid, new_molecule))
    {
        this->incrementMajorVersion();
        return true;
    }
    else
        return isDirty();
}

/** Return the index of the group in this forcefield

    \throw SireFF::invalid_group
*/
int InterGroupLJFF::groupIndex(FFBase::Group group) const
{
    if (group == groups().A())
        return 0;
    else if (group == groups().B())
        return 1;
    else
    {
        throw SireFF::invalid_group( QObject::tr(
                    "InterGroupLJFF can only place molecules into group A or B! %1")
                          .arg(group), CODELOC );

        return -1;
    }
}

/** Private class used by the "addTo" functions to actually add the
    molecule or part of molecule */
template<class T>
bool InterGroupLJFF::_pvt_addTo(FFBase::Group group, const T &mol, const ParameterMap &map)
{
    //get the index for this group
    int group_idx = this->groupIndex(group);

    //get the molecule's ID
    MoleculeID molid = mol.ID();

    this->assertValidGroup(group_idx, molid);

    ChangedLJMolecule new_molecule = changeRecord(molid);

    if (new_molecule.oldMolecule().isEmpty())
    {
        new_molecule = ChangedLJMolecule( LJMolecule(),
                        LJMolecule(mol, map.source(parameters().lj())) );

        //this is a freshly added molecule - add it to the current state
        this->addToCurrentState(group_idx, new_molecule.newMolecule());
        this->recordChange(group_idx, molid, new_molecule);

        this->incrementMajorVersion();

        return true;
    }
    else
    {
        new_molecule = new_molecule.add( mol, map.source(parameters().lj()) );

        if (this->applyChange(molid, new_molecule))
        {
            this->incrementMajorVersion();
            return true;
        }
        else
            return this->isDirty();
    }
}

/** Add the molecule 'molecule' to the group with index 'group' using the
    supplied map to find the forcefield parameters amongst the molecule's
    properties. Note that it is an error to try to add this to more than
    one group in the forcefield

    \throw SireMol::missing_property
    \throw SireMol::invalid_cast
    \throw SireFF::invalid_group
*/
bool InterGroupLJFF::addTo(FFBase::Group group, const Molecule &molecule,
                           const ParameterMap &map)
{
    return this->_pvt_addTo<Molecule>(group, molecule, map);
}

/** Add the residue 'residue' to the group with index 'group' using the
    supplied map to find the forcefield parameters amongst the residue's
    properties. Note that it is an error to try to add this to more than
    one group in the forcefield. Note that the property
    used must agree with the rest of the molecule if it is already
    in this forcefield.

    \throw SireMol::missing_property
    \throw SireMol::invalid_cast
    \throw SireFF::invalid_group
*/
bool InterGroupLJFF::addTo(FFBase::Group group, const Residue &residue,
                           const ParameterMap &map)
{
    return this->_pvt_addTo<Residue>(group, residue, map);
}

/** Add the atom 'atom' to the group with index 'group' using the
    supplied map to find the forcefield parameters amongst the atom's
    properties. Note that it is an error to try to add this to more than
    one group in the forcefield. Note that the property
    used must agree with the rest of the molecule if it is already
    in this forcefield.

    \throw SireMol::missing_property
    \throw SireMol::invalid_cast
    \throw SireFF::invalid_group
*/
bool InterGroupLJFF::addTo(FFBase::Group group, const NewAtom &atom,
                           const ParameterMap &map)
{
    return this->_pvt_addTo<NewAtom>(group, atom, map);
}

/** Add the selected atoms of 'molecule' to the group with index 'group' using the
    supplied map to find the forcefield parameters amongst the molecule's
    properties. Note that it is an error to try to add this to more than
    one group in the forcefield. Note that the property
    used must agree with the rest of the molecule if it is already
    in this forcefield.

    \throw SireMol::missing_property
    \throw SireMol::invalid_cast
    \throw SireFF::invalid_group
*/
bool InterGroupLJFF::addTo(FFBase::Group group, const Molecule &molecule,
                           const AtomSelection &selected_atoms,
                           const ParameterMap &map)
{
    //get the index for this group
    int group_idx = this->groupIndex(group);

    //get the molecule's ID
    MoleculeID molid = molecule.ID();

    this->assertValidGroup(group_idx, molid);

    ChangedLJMolecule new_molecule = changeRecord(molid);

    if (new_molecule.oldMolecule().isEmpty())
    {
        new_molecule = ChangedLJMolecule( LJMolecule(),
                                  LJMolecule( molecule, selected_atoms,
                                              map.source(this->parameters().lj()) ) );

        //this is a freshly added molecule - add it to the current state
        this->addToCurrentState(group_idx, new_molecule.newMolecule());

        //record the change
        this->recordChange(group_idx, molid, new_molecule);

        this->incrementMajorVersion();

        return true;
    }
    else
    {
        new_molecule = new_molecule.add( selected_atoms,
                                         map.source(this->parameters().lj()) );

        if (this->applyChange(molid, new_molecule))
        {
            this->incrementMajorVersion();
            return true;
        }
        else
            return isDirty();
    }
}
