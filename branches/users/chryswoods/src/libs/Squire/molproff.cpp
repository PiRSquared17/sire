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

#include "molproff.h"
#include "molproprocessor.h"
#include "molprosession.h"

#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"

#include "SireMM/atomiccharges.h"

#include "SireUnits/convert.h"
#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace Squire;
using namespace SireMM;
using namespace SireFF;
using namespace SireMol;
using namespace SireMaths;
using namespace SireUnits;
using namespace SireStream;

///////////
/////////// Implementation of MolproFF::Parameters
///////////

/** Constructor - by default the coulomb properties come from the 'charges'
    property */
MolproFF::Parameters::Parameters()
         : FFBase::Parameters(), coulomb_params("coulomb", "charges")
{}

/** Copy constructor */
MolproFF::Parameters::Parameters(const MolproFF::Parameters &other)
         : FFBase::Parameters(other), coulomb_params(other.coulomb_params)
{}

/** Destructor */
MolproFF::Parameters::~Parameters()
{}

/** Static object returned by MolproFF::parameters() */
MolproFF::Parameters MolproFF::Parameters::default_sources;

///////////
/////////// Implementation of MolproFF::Components
///////////

/** Constructor */
MolproFF::Components::Components() : FFBase::Components()
{}

/** Construct using the supplied forcefield */
MolproFF::Components::Components(const FFBase &ffbase, const Symbols &symbols)
         : FFBase::Components( ffbase, addToSymbols(symbols, x(), y(), z()) ),
           e_qm(ffbase, "qm", x(), y(), z())
{
    this->registerComponent(e_qm);
}

/** Copy constructor */
MolproFF::Components::Components(const MolproFF::Components &other)
         : FFBase::Components(other), e_qm(other.e_qm)
{}

/** Destructor */
MolproFF::Components::~Components()
{}

/** Assignment operator */
MolproFF::Components& MolproFF::Components::operator=(const MolproFF::Components &other)
{
    e_qm = other.e_qm;

    FFBase::Components::operator=(other);

    return *this;
}

/** Set the forcefield */
void MolproFF::Components::setForceField(const FFBase &ffbase)
{
    *this = MolproFF::Components(ffbase);
}

/** Describe the coulomb component */
QString MolproFF::Components::describe_qm()
{
    return QObject::tr("The total QM energy (includes the QM part of the QM/MM energy).");
}

///////////
/////////// Implementation of MolproFF
///////////

static const RegisterMetaType<MolproFF> r_molproff;

/** Serialise to a binary data stream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const MolproFF &molproff)
{
    writeHeader(ds, r_molproff, 1);
    
    SharedDataStream sds(ds);
    
    sds << molproff.spce << molproff.switchfunc
        << static_cast<const FFBase&>(molproff);
    
    //don't stream qm_coords and mm_coords_and_charges as
    //they can be rebuilt
    
    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, MolproFF &molproff)
{
    VersionID v = readHeader(ds, r_molproff);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> molproff.spce >> molproff.switchfunc
            >> static_cast<FFBase&>(molproff);
    
        //need to rebuild qm_coords and mm_coords_and_charges
    
    }
    else
        throw version_error(v, "1", r_molproff, CODELOC);
    
    return ds;
}

/** Construct an empty MolproFF */
MolproFF::MolproFF() : FFBase()
{
    this->registerComponents();
}

/** Construct a MolproFF with a specfied Space and switching function */
MolproFF::MolproFF(const Space &space, const SwitchingFunction &switchingfunction)
         : FFBase(), spce(space), switchfunc(switchingfunction)
{
    this->registerComponents();
}

/** Copy constructor */
MolproFF::MolproFF(const MolproFF &other)
         : FFBase(other), spce(other.spce), switchfunc(other.switchfunc),
           qm_coords(other.qm_coords), 
           mm_coords_and_charges(other.mm_coords_and_charges)
{
    //get the pointer from the base class...
    components_ptr = dynamic_cast<const MolproFF::Components*>( &(FFBase::components()) );
    BOOST_ASSERT( components_ptr != 0 );
}

/** Destructor */
MolproFF::~MolproFF()
{}

/** Register the components of this forcefield */
void MolproFF::registerComponents()
{
    std::auto_ptr<MolproFF::Components> ptr( new MolproFF::Components(*this) );

    FFBase::registerComponents(ptr.get());

    components_ptr = ptr.release();
}

/** Update the QM and MM arrays so that they contain the latest
    data, and return whether or not the arrays have changed since the
    last update.
*/
bool MolproFF::updateArrays()
{
    //if moved QM molecules
    //  return rebuild_all
    //else if moved MM molecules
    //  return rebuild_mm

    return false;
}

// /** Reconstruct the QM coordinate array */
// void MolproFF::reconstructQMArray()
// {
//     int nmols = qm_molecules.count();
//     const Molecule *molecule_array = qm_molecules.constData();
//
//     int nats = 0;
//
//     //get the number of atoms
//     for (int i=0; i<nmols; ++i)
//     {
//         nats += molecule_array[i].nAtoms();
//     }
//
//     if (nats == 0)
//     {
//         //clear the coordinates array
//         qm_coords.clear();
//         return;
//     }
//
//     //reserve sufficient space (3 times number of atoms)
//     QVector<double> coords( nats*3 );
//     double *coords_array = coords.data();
//
//     int atomid = 0;
//
//     //copy the atoms' coordinates, and convert to bohr radii
//     for (int i=0; i<nmols; ++i)
//     {
//         QVector<Vector> molcoords = molecule_array[i].coordinates();
//
//         int ncoords = molcoords.count();
//         const Vector *molcoords_array = molcoords.constData();
//
//         for (int j=0; j<ncoords; ++j)
//         {
//             const Vector &atmcoords = molcoords_array[j];
//
//             coords_array[atomid] = convertTo(atmcoords.x(), bohr_radii);
//             coords_array[atomid+1] = convertTo(atmcoords.y(), bohr_radii);
//             coords_array[atomid+2] = convertTo(atmcoords.z(), bohr_radii);
//
//             atomid += 3;
//         }
//     }
//
//     qm_coords = coords;
// }
//
// /** Reconstruct the MM coordinate and charge array */
// void MolproFF::reconstructMMArray()
// {
//     int nmols = mm_molecules.count();
//     const Molecule *molecule_array = mm_molecules.constData();
//     const QVector<ChargeParameter> *molcharge_array = mm_charges.constData();
//
//     //count the number of atoms that have a non-zero charge
//     //(same as the number of non-zero charge parameters...)
//     int nats = 0;
//
//     for (int i=0; i<nmols; ++i)
//     {
//         const QVector<ChargeParameter> &molcharges = molcharge_array[i];
//
//         int ncharges = molcharges.count();
//         const ChargeParameter *charge_array = molcharges.constData();
//
//         for (int j=0; j<nmols; ++j)
//         {
//             if ( not isZero(charge_array[j].charge()) )
//                ++nats;
//         }
//     }
//
//     if (nats == 0)
//     {
//         //there are no MM charges
//         mm_coords_and_charges.clear();
//         return;
//     }
//
//     //reserve sufficient space for the coordinates and charges (4*nats)
//     QVector<double> coords( nats*4 );
//     double *coords_array = coords.data();
//
//     int atomid = 0;
//
//     for (int i=0; i<nmols; ++i)
//     {
//         QVector<Vector> molcoords = molecule_array[i].coordinates();
//
//         int ncoords = molcoords.count();
//         const Vector *molcoords_array = molcoords.constData();
//
//         const ChargeParameter *charge_array = molcharge_array[i].constData();
//
//         for (int j=0; j<ncoords; ++j)
//         {
//             if ( not isZero(charge_array[j].charge()) )
//             {
//                 const Vector &atmcoords = molcoords_array[j];
//
//                 //copy (and convert) the coordinates
//                 coords_array[atomid] = convertTo(atmcoords.x(), bohr_radii);
//                 coords_array[atomid+1] = convertTo(atmcoords.y(), bohr_radii);
//                 coords_array[atomid+2] = convertTo(atmcoords.z(), bohr_radii);
//
//                 //copy (and convert) the charge
//                 coords_array[atomid+3] = convertTo(charge_array[j].charge(),
//                                                    mod_electrons);
//
//                 atomid += 4;
//             }
//         }
//     }
//
//     mm_coords_and_charges = coords;
// }

/** Add a molecule to the QM region

    \throw SireMol::duplicate_molecule
*/
void MolproFF::_pvt_addToQM(const Molecule &molecule)
{
    MoleculeID molid = molecule.ID();

    //check that the molecule is not in the MM region

    //add the molecule to the qm list, and also to the changed list
}

/** Add a molecule to the MM region.

    \throw SireMol::duplicate_molecule
    \throw SireError::incompatible_error
*/
void MolproFF::_pvt_addToMM(const Molecule &molecule, const ParameterMap &map)
{
    //get the charges for this molecule from the coulomb property
    AtomicCharges charges = molecule.getProperty( map.source(parameters().coulomb()) );

    MoleculeID molid = molecule.ID();

    //check that the molecule is not in the QM region


    //add the molecule to the list, and also add it to the changed list

}

/** Add a molecule to the QM region

    \throw SireMol::duplicate_molecule
*/
void MolproFF::addToQM(const Molecule &molecule)
{
    this->_pvt_addToQM(molecule);

    this->incrementMajorVersion();
}

/** Add a load of molecules to the QM region

    \throw SireMol::duplicate_molecule
*/
void MolproFF::addToQM(const QList<Molecule> &molecules)
{
    // Add the molecules to a copy of this forcefield - this
    // is to maintain the invariant

    MolproFF copy(*this);

    for (QList<Molecule>::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        copy._pvt_addToQM(*it);
    }

    copy.incrementMajorVersion();

    //everything's ok - copy back to the original
    *this = copy;
}

/** Add a molecule to the MM region.

    \throw SireMol::duplicate_molecule
*/
void MolproFF::addToMM(const Molecule &molecule, const ParameterMap &map)
{
    this->_pvt_addToMM(molecule, map);
    this->incrementMajorVersion();
}

/** Add a load of molecules to the MM region

    \throw SireMol::duplicate_molecule
*/
void MolproFF::addToMM(const QList<Molecule> &molecules,
                       const ParameterMap &map)
{
    //add the Molecules to a copy of this forcefield -
    // this maintains the invariant

    MolproFF copy(*this);

    for (QList<Molecule>::const_iterator it = molecules.begin();
         it != molecules.end();
         ++it)
    {
        copy._pvt_addToMM(*it, map);
    }

    copy.incrementMajorVersion();

    //everything's ok - copy back to the original
    *this = copy;
}

/** Does nothing... */
bool MolproFF::change(const Molecule&)
{
    return false;
}

/** Does nothing... */
bool MolproFF::change(const Residue&)
{
    return false;
}

/** Does nothing... */
bool MolproFF::change(const NewAtom&)
{
    return false;
}

/** Use the passed MolproSession to recalculate the energy of
    this forcefield - this throws an exception if the session is
    incompatible with this forcefield

    \throw SireError::incompatible_error
*/
Values MolproFF::recalculateEnergy(MolproSession &session)
{
    session.assertCompatibleWith(*this);

    //update the QM and MM arrays to account for any moves
    if (this->updateArrays())
    {
        //the arrays have changed! - pass them to the session
        session.setArrays(qm_coords, mm_coords_and_charges);

        //calculate the HF energy of the system
        double hf_nrg = session.calculateEnergy("hf");

        setComponent( components().total(), hf_nrg );
        setComponent( components().qm(), hf_nrg );
    }

    return currentEnergies();
}

/** The MolproFF forcefield is unable to recalculate its own
    energy - it has to be placed on a MolproProcessor to be
    able to do this. */
void MolproFF::recalculateEnergy()
{
    throw SireError::program_bug( QObject::tr(
            "It is not possible to calculate the energy of a MolproFF "
            "directly. You must place the MolproFF forcefield onto a "
            "MolproProcessor and ask the MolproProcessor to calculate "
            "the energy."), CODELOC );
}
