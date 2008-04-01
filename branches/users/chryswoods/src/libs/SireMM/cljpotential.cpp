/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "SireBase/sparsematrix.hpp"

#include "cljpotential.h"
#include "ljparameter.h"
#include "switchingfunction.h"

#include "SireMol/mover.hpp"
#include "SireMol/atomcoords.h"

#include "SireVol/cartesian.h"

#include "SireMaths/maths.h"

#include "SireUnits/units.h"

#include "SireBase/errors.h"
#include "SireError/errors.h"
#include "SireFF/errors.h"

#include "SireStream/datastream.h"

using namespace SireMM;
using namespace SireMM::detail;

using namespace SireFF;
using namespace SireFF::detail;

using namespace SireMol;
using namespace SireVol;

using namespace SireMaths;

using namespace SireBase;

using namespace SireStream;

///////
/////// Completely instantiate the CLJPotential ancillary classes
///////

template
class AtomicParameters3D<CLJParameter>;

template
class IntraScaledParameters<CLJNBPairs>;

template
class IntraScaledAtomicParameters< AtomicParameters3D<CLJParameter>,
                                   IntraScaledParameters<CLJNBPairs> >;

template
class FFMolecule3D<InterCLJPotential>;

template
class FFMolecules3D<InterCLJPotential>;

template
class ChangedMolecule<InterCLJPotential::Molecule>;

template
class FFMolecule3D<IntraCLJPotential>;

template
class FFMolecules3D<IntraCLJPotential>;

template
class ChangedMolecule<IntraCLJPotential::Molecule>;

/** Internal function used to get the charge and LJ parameters from a molecule
    and convert them into a PackedArray of reduced charges and LJ parameter IDs */
static PackedArray2D<CLJParameter> getCLJParameters(const PartialMolecule &molecule,
                                                    const PropertyName &charge_property,
                                                    const PropertyName &lj_property)
{
    const AtomCharges &chgs = molecule.property(charge_property)->asA<AtomCharges>();
    const AtomLJs &ljs = molecule.property(lj_property)->asA<AtomLJs>();
    
    const AtomSelection &selected_atoms = molecule.selection();
    
    if (selected_atoms.selectedNone())
        return PackedArray2D<CLJParameter>();
    
    //create space for the parameters - only need space for CutGroups
    //that contain at least one selected atom
    QVector< QVector<CLJParameter> > cljparams( selected_atoms.nSelectedCutGroups() );
    QVector<CLJParameter>* cljparams_array = cljparams.data();

    const double sqrt_4pieps0 = std::sqrt(SireUnits::one_over_four_pi_eps0);

    try
    {

    LJParameterDB::lock();

    if (selected_atoms.selectedAllCutGroups())
    {
        const int ncg = molecule.data().info().nCutGroups();
    
        for (CGIdx i(0); i<ncg; ++i)
        {
            const int nats = molecule.data().info().nAtoms(i);
            
            QVector<CLJParameter> group_cljs(nats);
            CLJParameter *group_cljs_array = group_cljs.data();
            
            //get the arrays containing the charge and LJ parameters
            //for this CutGroup
            const SireUnits::Dimension::Charge *group_chgs = chgs.constData(i);
            const LJParameter *group_ljs = ljs.constData(i);
            
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<nats; ++j)
                {
                    group_cljs_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                    group_cljs_array[j].ljid = 
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    group_cljs_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                    group_cljs_array[j].ljid =
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            
            cljparams_array[i] = group_cljs;
        }
    }
    else
    {
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            const int nats = molecule.data().info().nAtoms(i);
            
            QVector<CLJParameter> group_cljs(nats);
            CLJParameter *group_cljs_array = group_cljs.data();
            
            //get the arrays containing the charge and LJ parameters
            //for this CutGroup
            const SireUnits::Dimension::Charge *group_chgs = chgs.constData(i);
            const LJParameter *group_ljs = ljs.constData(i);
            
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<nats; ++j)
                {
                    group_cljs_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                    group_cljs_array[j].ljid = 
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    group_cljs_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                    group_cljs_array[j].ljid =
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            
            cljparams_array[i] = group_cljs;
        }
    }
    
    LJParameterDB::unlock();
    
    }
    catch(...)
    {
        LJParameterDB::unlock();
        throw;
    }
    
    return PackedArray2D<CLJParameter>( cljparams );
}

/////////////
///////////// Implementation of CLJPotential
/////////////

static const Space default_space = Cartesian();
static const SwitchingFunction default_switchfunc = NoCutoff();

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const CLJPotential &cljpot)
{
    ds << cljpot.props;
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      CLJPotential &cljpot)
{
    ds >> cljpot.props;
    
    //extract all of the properties
    cljpot.spce = cljpot.props.property("space")->asA<SpaceBase>();
    cljpot.switchfunc = cljpot.props.property("switchingFunction")
                                          ->asA<SwitchFunc>();
    
    cljpot.combining_rules = LJParameterDB::interpret(
                                cljpot.props.property("combiningRules")
                                    ->asA<VariantProperty>().convertTo<QString>() );

    cljpot.use_electrostatic_shifting = cljpot.props.property("shiftElectrostatics")
                                    ->asA<VariantProperty>().convertTo<bool>();
                                    
    cljpot.need_update_ljpairs = true;
    
    return ds;
}

/** Constructor */
CLJPotential::CLJPotential()
             : spce(default_space), switchfunc(default_switchfunc),
               combining_rules( LJParameterDB::interpret("arithmetic") ),
               need_update_ljpairs(true), use_electrostatic_shifting(false)
{
    //record the defaults
    props.setProperty( "space", spce );
    props.setProperty( "switchingFunction", switchfunc );
    props.setProperty( "combiningRules", 
                       VariantProperty( LJParameterDB::toString(combining_rules) ) );
    props.setProperty( "shiftElectrostatics",
                       VariantProperty(use_electrostatic_shifting) );
}

/** Copy constructor */
CLJPotential::CLJPotential(const CLJPotential &other)
             : ljpairs(other.ljpairs), props(other.props),
               spce(other.spce), switchfunc(other.switchfunc),
               combining_rules(other.combining_rules),
               need_update_ljpairs(other.need_update_ljpairs),
               use_electrostatic_shifting(other.use_electrostatic_shifting)
{}

/** Destructor */
CLJPotential::~CLJPotential()
{}

/** Copy assignment operator */
CLJPotential& CLJPotential::operator=(const CLJPotential &other)
{
    if (this != &other)
    {
        ljpairs = other.ljpairs;
        props = other.props;
        spce = other.spce;
        switchfunc = other.switchfunc;
        combining_rules = other.combining_rules;
        need_update_ljpairs = other.need_update_ljpairs;
        use_electrostatic_shifting = other.use_electrostatic_shifting;
    }
    
    return *this;
}

/** You need to call this function before you start a block of 
    energy of force evaluation using this forcefield. You should
    also call 'finishedEvaluation()' once you have finished. */
void CLJPotential::startEvaluation()
{
    if (need_update_ljpairs)
    {
        //get the LJPairs array from the LJParameterDB
        ljpairs = LJParameterDB::getLJPairs(combining_rules);
        need_update_ljpairs = false;
    }
}

/** You should call this function once you have finished a block of
    force or energy evaluation using this potential */
void CLJPotential::finishedEvaluation()
{}

/** Return all of the properties set in this forcefield */
const Properties& CLJPotential::properties() const
{
    return props;
}

/** Return whether or not this potential has a property called 'name' */
bool CLJPotential::containsProperty(const QString &name) const
{
    return props.hasProperty(name);
}

/** Return the property with name 'name'

    \throw SireBase::missing_property
*/
const Property& CLJPotential::property(const QString &name) const
{
    return props.property(name);
}

/** Set the 3D space in which the molecules in this potential are evaluated */
bool CLJPotential::setSpace(const Space &new_space)
{
    if (spce != new_space)
    {
        spce = new_space;
        props.setProperty( "space", new_space );
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set the switching function used to scale the interactions between
    CutGroups to zero at the cutoff */
bool CLJPotential::setSwitchingFunction(const SwitchingFunction &new_switchfunc)
{
    if (switchfunc != new_switchfunc)
    {
        switchfunc = new_switchfunc;
        props.setProperty( "switchingFunction", new_switchfunc );
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set whether or not to shift the electrostatics between CutGroups so that
    the group-group net electrostatic interaction energy between CutGroups
    is zero at the cutoff */
bool CLJPotential::setShiftElectrostatics(bool switchelectro)
{
    if (switchelectro != use_electrostatic_shifting)
    {
        use_electrostatic_shifting = switchelectro;
        props.setProperty( "shiftElectrostatics", VariantProperty(switchelectro) );
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set the combining rules to use to obtain mixed LJ parameters */
bool CLJPotential::setCombiningRules(const QString &combiningrules)
{
    LJParameterDB::CombiningRules new_rules = LJParameterDB::interpret(combiningrules);
    
    if (new_rules != combining_rules)
    {
        combining_rules = new_rules;
        need_update_ljpairs = true;
        props.setProperty( "combiningRules", VariantProperty(combiningrules) );
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set the property 'name' to the value 'value'. Returns whether or not
    this changes this forcefield.

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
bool CLJPotential::setProperty(const QString &name, const PropertyBase &value)
{
    if (name == QLatin1String("space"))
    {
        return this->setSpace( value.asA<SpaceBase>() );
    }
    else if (name == QLatin1String("switchingFunction"))
    {
        return this->setSwitchingFunction( value.asA<SwitchFunc>() );
    }
    else if (name == QLatin1String("shiftElectrostatics"))
    {
        return this->setShiftElectrostatics( value.asA<VariantProperty>()
                                                     .convertTo<bool>() );
    }
    else if (name == QLatin1String("combiningRules"))
    {
        return this->setCombiningRules( value.asA<VariantProperty>()
                                                     .convertTo<QString>() );
    }
    else
        throw SireBase::missing_property( QObject::tr(
            "The CLJ potentials do not have a property called \"%1\" that "
            "can be changed. Available properties are [ %2 ].")
                .arg(name, QStringList(props.propertyKeys()).join(", ")), CODELOC );
                
    return false;
}

/** Return the 3D space in which this potential is evaluated */
const SpaceBase& CLJPotential::space() const
{
    return *spce;
}

/** Return the switching function used to scale the group-group
    interactions to zero at the cutoff */
const SwitchFunc& CLJPotential::switchingFunction() const
{
    return *switchfunc;
}

/** Return whether or not the net group-group electrostatic interaction
    energy is shifted so that it is zero at the cutoff */
bool CLJPotential::shiftElectrostatics() const
{
    return use_electrostatic_shifting;
}

/** Return the string identifying the combining rules used to 
    obtain the mixed LJ parameters */
const QString& CLJPotential::combiningRules() const
{
    return LJParameterDB::toString(combining_rules);
}

/////////////
///////////// Implementation of InterCLJPotential
/////////////

static const RegisterMetaType<InterCLJPotential> r_interclj( MAGIC_ONLY,
                                            InterCLJPotential::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const InterCLJPotential &interclj)
{
    writeHeader(ds, r_interclj, 1);
    
    ds << static_cast<const CLJPotential&>(interclj);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      InterCLJPotential &interclj)
{
    VersionID v = readHeader(ds, r_interclj);
    
    if (v == 1)
    {
        ds >> static_cast<CLJPotential&>(interclj);
    }
    else
        throw version_error(v, "1", r_interclj, CODELOC);
        
    return ds;
}

/** Constructor */
InterCLJPotential::InterCLJPotential() : CLJPotential()
{}

/** Copy constructor */
InterCLJPotential::InterCLJPotential(const InterCLJPotential &other)
                  : CLJPotential(other)
{}

/** Destructor */
InterCLJPotential::~InterCLJPotential()
{}

/** Copy assignment operator */
InterCLJPotential& InterCLJPotential::operator=(const InterCLJPotential &other)
{
    CLJPotential::operator=(other);
    return *this;
}

void InterCLJPotential::throwMissingForceComponent(const Symbol &symbol,
                              const IntraCLJPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no force component in potential %1 - available "
        "components are %1, %2 and %3.")
            .arg(this->what())
            .arg(components.total().toString(), components.coulomb().toString(),
                 components.lj().toString()), CODELOC );
}

/** Return all of the parameters needed by this potential for 
    the molecule 'molecule', using the supplied property map to
    find the properties that contain those parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCLJPotential::Parameters 
InterCLJPotential::getParameters(const PartialMolecule &molecule,
                                 const PropertyMap &map)
{
    return Parameters( molecule, map[parameters().coordinates()],
                       getCLJParameters(molecule, map[parameters().charge()],
                                        map[parameters().lj()]) );
}

/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCLJPotential::Parameters
InterCLJPotential::updateParameters(const InterCLJPotential::Parameters &old_params,
                                    const PartialMolecule &old_molecule,
                                    const PartialMolecule &new_molecule,
                                    const PropertyMap &map)
{
    if (old_molecule.selection() != new_molecule.selection())
        //the selection has changed - just get completely new parameters
        return this->getParameters(new_molecule, map);

    Parameters new_params = old_params;

    //get the property names
    const PropertyName &coords_property = map[parameters().coordinates()];
    const PropertyName &chg_property = map[parameters().charge()];
    const PropertyName &lj_property = map[parameters().lj()];
    
    //get what has changed
    bool new_coords = old_molecule.version(coords_property) !=
                         new_molecule.version(coords_property);
                             
    bool new_clj = ( old_molecule.version(chg_property) !=
                         new_molecule.version(chg_property) ) or
                   ( old_molecule.version(lj_property) !=
                         new_molecule.version(lj_property) );

    if (new_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        coords_property) );

    if (new_clj)
        new_params.setAtomicParameters( getCLJParameters(new_molecule,
                                                         chg_property, lj_property) );

    return new_params;
}
                 
/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', also while the parameters of 'old_molecule'
    where found in 'old_map', now get the parameters using 'new_map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCLJPotential::Parameters
InterCLJPotential::updateParameters(const InterCLJPotential::Parameters &old_params,
                                    const PartialMolecule &old_molecule,
                                    const PartialMolecule &new_molecule,
                                    const PropertyMap &old_map, 
                                    const PropertyMap &new_map)
{
    if (old_molecule.selection() != new_molecule.selection())
        //the selection has changed - just get completely new parameters
        return this->getParameters(new_molecule, new_map);

    Parameters new_params = old_params;

    //get the property names
    const PropertyName &old_coords = old_map[parameters().coordinates()];
    const PropertyName &old_chg = old_map[parameters().charge()];
    const PropertyName &old_lj = old_map[parameters().lj()];
    
    const PropertyName &new_coords = new_map[parameters().coordinates()];
    const PropertyName &new_chg = new_map[parameters().charge()];
    const PropertyName &new_lj = new_map[parameters().lj()];
    
    //get what has changed
    bool changed_coords = (new_coords != old_coords) or
                           old_molecule.version(old_coords) !=
                           new_molecule.version(old_coords);
                             
    bool changed_clj = (new_chg != old_chg or new_lj != old_lj) or
                       ( old_molecule.version(old_chg) !=
                         new_molecule.version(old_chg) ) or
                       ( old_molecule.version(old_lj) !=
                         new_molecule.version(old_lj) );

    if (changed_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        new_coords) );

    if (changed_clj)
        new_params.setAtomicParameters( getCLJParameters(new_molecule,
                                                         new_chg, new_lj) );

    return new_params;
}

/** Return the InterCLJPotential::Molecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCLJPotential::Molecule
InterCLJPotential::parameterise(const PartialMolecule &molecule,
                                const PropertyMap &map)
{
    return InterCLJPotential::Molecule(molecule, *this, map);
}

/** Concert the passed group of molecules into InterCLJPotential::Molecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCLJPotential::Molecules 
InterCLJPotential::parameterise(const MolGroup &molecules,
                                const PropertyMap &map)
{
    return InterCLJPotential::Molecules(molecules, *this, map);
}

/** Return the total charge of the parameters for the group in 'params' */
double InterCLJPotential::totalCharge(
                        const InterCLJPotential::Parameters::Array &params) const
{
    int n = params.count();
    const Parameter *params_array = params.constData();
    
    double chg = 0;
    
    for (int i=0; i<n; ++i)
    {
        chg += params_array[i].reduced_charge;
    }
    
    return chg;
}

/** Calculate the coulomb and LJ energy between the passed pair
    of molecules and add these energies onto 'energy'. This uses
    the passed workspace to perform the calculation */
void InterCLJPotential::_pvt_calculateEnergy(const InterCLJPotential::Molecule &mol0,
                                             const InterCLJPotential::Molecule &mol1,
                                             InterCLJPotential::Energy &energy,
                                             InterCLJPotential::EnergyWorkspace &distmat,
                                             double scale_energy) const
{
    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    BOOST_ASSERT( mol0.parameters().atomicParameters().count() == int(ngroups0) );
    BOOST_ASSERT( mol1.parameters().atomicParameters().count() == int(ngroups1) );
    
    const Parameters::Array *molparams0_array 
                                = mol0.parameters().atomicParameters().constData();
    const Parameters::Array *molparams1_array 
                                = mol1.parameters().atomicParameters().constData();
    
    double cnrg = 0;
    double ljnrg = 0;
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is either the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDist(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;
                
            double icnrg = 0;
            double iljnrg = 0;
            
            //loop over all interatomic pairs and calculate the energies
            const quint32 nats1 = group1.count();
            const Parameter *params1_array = params1.constData();
            
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        icnrg += param0.reduced_charge * 
                                 params1_array[j].reduced_charge / distmat[j];
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                        
                        const double invdist = double(1) / distmat[j];
                        
                        icnrg += param0.reduced_charge * param1.reduced_charge 
                                      * invdist;
                              
                        if (param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                      ljpairs.map(param0.ljid,
                                                                  param1.ljid)];
                        
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            iljnrg += ljpair.epsilon() * (sig_over_dist12 - 
                                                          sig_over_dist6);
                        }
                    }
                }
            }
            
            //are we shifting the electrostatic potential?
            if (use_electrostatic_shifting)
                icnrg -= this->totalCharge(params0) * this->totalCharge(params1)
                              / switchfunc->electrostaticCutoffDistance();
            
            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            if (mindist > switchfunc->featherDistance())
            {
                cnrg += switchfunc->electrostaticScaleFactor(mindist) * icnrg;
                ljnrg += switchfunc->vdwScaleFactor(mindist) * iljnrg;
            }
            else
            {
                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
    }
    
    //add this molecule pair's energy onto the total
    energy += Energy(scale_energy * cnrg, scale_energy * ljnrg);
}

/** Add to the forces in 'forces0' the forces acting on 'mol0' caused
    by 'mol1' */
void InterCLJPotential::_pvt_calculateForce(const InterCLJPotential::Molecule &mol0, 
                                            const InterCLJPotential::Molecule &mol1,
                                            MolForceTable &forces0, 
                                            InterCLJPotential::ForceWorkspace &distmat,
                                            double scale_force) const
{
    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == forces0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == forces0.molNum() );

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    BOOST_ASSERT(mol0.parameters().atomicParameters().count() == int(ngroups0));
    BOOST_ASSERT(mol1.parameters().atomicParameters().count() == int(ngroups1));
    
    const Parameters::Array *molparams0_array 
                                    = mol0.parameters().atomicParameters().constData();
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    const MolForceTable::Array *forces0_array = forces0.constData();
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        //get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the forces array
        int force0_idx = forces0.map(cgidx_igroup);
        
        if (force0_idx == -1)
            //there is no space for the forces on this CutGroup in 
            //the forcetable - were are therefore not interested in
            //this CutGroup
            continue;

        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        BOOST_ASSERT(forces0_array[force0_idx].count() == int(nats0));
    
        Vector *group_forces0_array = forces0.data(force0_idx);

        //ok, we are interested in the forces acting on this CutGroup
        // - calculate all of the forces on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the forces taking into account
                //the derivative of the switching function!
                
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_coul = switchfunc->electrostaticScaleFactor(mindist);
                const double scl_lj = switchfunc->vdwScaleFactor(mindist);
                
                Vector group_sep = (group1.aaBox().center() - aabox0.center())
                                        .normalise();
                
                Vector dscl_coul = switchfunc->dElectrostaticScaleFactor(mindist) 
                                     * group_sep;
                                     
                Vector dscl_lj = switchfunc->dVDWScaleFactor(mindist)
                                     * group_sep;
                
                double shift_coul = 0;

                if (use_electrostatic_shifting)
                    shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                                    / switchfunc->electrostaticCutoffDistance();
                
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];
                
                    Vector total_force;
                
                    if (param0.ljid == 0)
                    {
                        //null LJ parameter - only add on the coulomb energy
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const double q2 = param0.reduced_charge *
                                              params1_array[j].reduced_charge;
                            
                            if (q2 != 0)
                            {
                                //calculate the coulomb energy
                                const double cnrg = q2 / distmat[j].length();
                                               
                                //calculate the coulomb force
                                Vector cforce = (scl_coul * -cnrg / distmat[j].length() *
                                                 distmat[j].direction()) +
                                             
                                                ((cnrg-shift_coul) * dscl_coul);
                        
                                total_force -= cforce;
                            }
                        }
                    }
                    else
                    {
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            //do both coulomb and LJ
                            const Parameter &param1 = params1_array[j];
                        
                            const double invdist = double(1) / distmat[j].length();
                        
                            Vector force;
                            
                            const double q2 = param0.reduced_charge *
                                              param1.reduced_charge;
                        
                            if (q2 != 0)
                            {
                                //calculate the energy
                                const double cnrg = q2 * invdist;
                        
                                //calculate the force
                                force = (scl_coul * -cnrg / distmat[j].length() *
                                         distmat[j].direction()) +
                                             
                                         ((cnrg-shift_coul) * dscl_coul);
                            }
                              
                            if (param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                //calculate the energy
                                const double ljnrg = ljpair.epsilon() *
                                                      (sig_over_dist12 - sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                force += ((scl_lj * ljpair.epsilon() * 
                                            (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                            * distmat[j].direction())
                                            
                                          + (ljnrg * dscl_lj);
                            }

                            total_force -= force;
                        }
                    }
                    
                    //update the forces array
                    group_forces0_array[i] += scale_force * total_force;
                }
            }
            else
            {
                //not in the feather region, so can calculate the forces
                //directly (also, no need to calculate shift, as 
                //the shifting function is constant, so does not
                //affect the gradient)
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];

                    Vector total_force;
                
                    if (param0.ljid == 0)
                    {
                        //null LJ parameter - only add on the coulomb energy
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const double q2 = param0.reduced_charge * 
                                              params1_array[j].reduced_charge;
                        
                            //calculate the coulomb force
                            if (q2 != 0)
                            {
                                Vector cforce = -(q2 / distmat[j].length2()) *
                                                    distmat[j].direction();
                        
                                total_force -= cforce;
                            }
                        }
                    }
                    else
                    {
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            //do both coulomb and LJ
                            const Parameter &param1 = params1_array[j];
                        
                            const double invdist = double(1) / distmat[j].length();
                            const double invdist2 = pow_2(invdist);
                        
                            //calculate the force
                            Vector force = -(param0.reduced_charge * 
                                             param1.reduced_charge * invdist2) 
                                            
                                            * distmat[j].direction();
                              
                            if (param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                force += (ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                              12.0*sig_over_dist12))
                                        * distmat[j].direction();
                            }
                        
                            total_force -= force;
                        }
                    }
                    
                    group_forces0_array[i] += scale_force * total_force;

                } // end of loop over i atoms

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

/** Add to the forces in 'forces0' the forces acting on 'mol0' caused
    by 'mol1' */
void InterCLJPotential::_pvt_calculateCoulombForce(
                                            const InterCLJPotential::Molecule &mol0, 
                                            const InterCLJPotential::Molecule &mol1,
                                            MolForceTable &forces0, 
                                            InterCLJPotential::ForceWorkspace &distmat,
                                            double scale_force) const
{
    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == forces0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == forces0.molNum() );

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    BOOST_ASSERT(mol0.parameters().atomicParameters().count() == int(ngroups0));
    BOOST_ASSERT(mol1.parameters().atomicParameters().count() == int(ngroups1));
    
    const Parameters::Array *molparams0_array 
                                    = mol0.parameters().atomicParameters().constData();
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    const MolForceTable::Array *forces0_array = forces0.constData();
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        //get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the forces array
        int force0_idx = forces0.map(cgidx_igroup);
        
        if (force0_idx == -1)
            //there is no space for the forces on this CutGroup in 
            //the forcetable - were are therefore not interested in
            //this CutGroup
            continue;

        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        BOOST_ASSERT(forces0_array[force0_idx].count() == int(nats0));
    
        Vector *group_forces0_array = forces0.data(force0_idx);

        //ok, we are interested in the forces acting on this CutGroup
        // - calculate all of the forces on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the forces taking into account
                //the derivative of the switching function!
                
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_coul = switchfunc->electrostaticScaleFactor(mindist);
                
                Vector group_sep = (group1.aaBox().center() - aabox0.center())
                                        .normalise();
                
                Vector dscl_coul = switchfunc->dElectrostaticScaleFactor(mindist) 
                                     * group_sep;
                
                double shift_coul = 0;

                if (use_electrostatic_shifting)
                    shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                                    / switchfunc->electrostaticCutoffDistance();
                
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];
                
                    Vector total_force;
                
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double q2 = param0.reduced_charge *
                                          params1_array[j].reduced_charge;
                          
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = q2 / distmat[j].length();
                                               
                            //calculate the coulomb force
                            Vector cforce = (scl_coul * -cnrg / distmat[j].length() *
                                             distmat[j].direction()) +
                                             
                                            ((cnrg-shift_coul) * dscl_coul);
                        
                            total_force -= cforce;
                        }
                    }
                    
                    //update the forces array
                    group_forces0_array[i] += scale_force * total_force;
                }
            }
            else
            {
                //not in the feather region, so can calculate the forces
                //directly (also, no need to calculate shift, as 
                //the shifting function is constant, so does not
                //affect the gradient)
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];

                    Vector total_force;
                
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double q2 = param0.reduced_charge * 
                                          params1_array[j].reduced_charge;
                        
                        //calculate the coulomb force
                        if (q2 != 0)
                        {
                            Vector cforce = -(q2 / distmat[j].length2()) *
                                                distmat[j].direction();
                        
                            total_force -= cforce;
                        }
                    }
                    
                    group_forces0_array[i] += scale_force * total_force;

                } // end of loop over i atoms

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

/** Add to the forces in 'forces0' the forces acting on 'mol0' caused
    by 'mol1' */
void InterCLJPotential::_pvt_calculateLJForce(
                                            const InterCLJPotential::Molecule &mol0, 
                                            const InterCLJPotential::Molecule &mol1,
                                            MolForceTable &forces0, 
                                            InterCLJPotential::ForceWorkspace &distmat,
                                            double scale_force) const
{
    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == forces0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == forces0.molNum() );

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    BOOST_ASSERT(mol0.parameters().atomicParameters().count() == int(ngroups0));
    BOOST_ASSERT(mol1.parameters().atomicParameters().count() == int(ngroups1));
    
    const Parameters::Array *molparams0_array 
                                    = mol0.parameters().atomicParameters().constData();
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    const MolForceTable::Array *forces0_array = forces0.constData();
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        //get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the forces array
        int force0_idx = forces0.map(cgidx_igroup);
        
        if (force0_idx == -1)
            //there is no space for the forces on this CutGroup in 
            //the forcetable - were are therefore not interested in
            //this CutGroup
            continue;

        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        BOOST_ASSERT(forces0_array[force0_idx].count() == int(nats0));
    
        Vector *group_forces0_array = forces0.data(force0_idx);

        //ok, we are interested in the forces acting on this CutGroup
        // - calculate all of the forces on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the forces taking into account
                //the derivative of the switching function!
                
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_lj = switchfunc->vdwScaleFactor(mindist);
                
                Vector group_sep = (group1.aaBox().center() - aabox0.center())
                                        .normalise();
                
                Vector dscl_lj = switchfunc->dVDWScaleFactor(mindist)
                                     * group_sep;
                
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];
                
                    Vector total_force;
                
                    if (param0.ljid != 0)
                    {
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const Parameter &param1 = params1_array[j];

                            if (param1.ljid != 0)
                            {
                                const double invdist = double(1) / distmat[j].length();

                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                //calculate the energy
                                const double ljnrg = ljpair.epsilon() *
                                                      (sig_over_dist12 - sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                Vector force = ((scl_lj * ljpair.epsilon() * 
                                            (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                               * distmat[j].direction())
                                            
                                             + (ljnrg * dscl_lj);

                                total_force -= force;
                            }
                        }
                    }
                    
                    //update the forces array
                    group_forces0_array[i] += scale_force * total_force;
                }
            }
            else
            {
                //not in the feather region, so can calculate the forces
                //directly
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];

                    Vector total_force;
                
                    if (param0.ljid != 0)
                    {
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            //do both coulomb and LJ
                            const Parameter &param1 = params1_array[j];
                              
                            if (param1.ljid != 0)
                            {
                                const double invdist = double(1) / distmat[j].length();

                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                Vector force = (ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                                12.0*sig_over_dist12))
                                                * distmat[j].direction();

                                total_force -= force;
                            }
                        }
                    }
                    
                    group_forces0_array[i] += scale_force * total_force;

                } // end of loop over i atoms

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

/////////////
///////////// Implementation of IntraCLJPotential
/////////////

static const RegisterMetaType<IntraCLJPotential> r_intraclj( MAGIC_ONLY,
                                            IntraCLJPotential::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const IntraCLJPotential &intraclj)
{
    writeHeader(ds, r_intraclj, 1);
    
    ds << static_cast<const CLJPotential&>(intraclj);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      IntraCLJPotential &intraclj)
{
    VersionID v = readHeader(ds, r_intraclj);
    
    if (v == 1)
    {
        ds >> static_cast<CLJPotential&>(intraclj);
    }
    else
        throw version_error(v, "1", r_intraclj, CODELOC);
        
    return ds;
}

/** Constructor */
IntraCLJPotential::IntraCLJPotential() : CLJPotential()
{}

/** Copy constructor */
IntraCLJPotential::IntraCLJPotential(const IntraCLJPotential &other)
                  : CLJPotential(other)
{}

/** Destructor */
IntraCLJPotential::~IntraCLJPotential()
{}

/** Copy assignment operator */
IntraCLJPotential& IntraCLJPotential::operator=(const IntraCLJPotential &other)
{
    CLJPotential::operator=(other);
    return *this;
}

void IntraCLJPotential::throwMissingForceComponent(const Symbol &symbol,
                              const IntraCLJPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no force component in potential %1 - available "
        "components are %1, %2 and %3.")
            .arg(this->what())
            .arg(components.total().toString(), components.coulomb().toString(),
                 components.lj().toString()), CODELOC );
}

/** Return all of the parameters needed by this potential for 
    the molecule 'molecule', using the supplied property map to
    find the properties that contain those parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCLJPotential::Parameters 
IntraCLJPotential::getParameters(const PartialMolecule &molecule,
                                 const PropertyMap &map)
{
    return Parameters( AtomicParameters3D<CLJParameter>(
                               molecule, map[parameters().coordinates()],
                               getCLJParameters(molecule, 
                                                map[parameters().charge()],
                                                map[parameters().lj()]) ),
              IntraScaledParameters<CLJNBPairs>(
                               molecule, map[parameters().intraScaleFactors()] )
                     );
}

/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCLJPotential::Parameters
IntraCLJPotential::updateParameters(const IntraCLJPotential::Parameters &old_params,
                                    const PartialMolecule &old_molecule,
                                    const PartialMolecule &new_molecule,
                                    const PropertyMap &map)
{
    if (old_molecule.selection() != new_molecule.selection())
        //the selection has changed - just get completely new parameters
        return this->getParameters(new_molecule, map);

    Parameters new_params = old_params;

    //get the property names
    const PropertyName &coords_property = map[parameters().coordinates()];
    const PropertyName &chg_property = map[parameters().charge()];
    const PropertyName &lj_property = map[parameters().lj()];
    const PropertyName &scl_property = map[parameters().intraScaleFactors()];
    
    //get what has changed
    bool new_coords = old_molecule.version(coords_property) !=
                         new_molecule.version(coords_property);
                             
    bool new_clj = ( old_molecule.version(chg_property) !=
                         new_molecule.version(chg_property) ) or
                   ( old_molecule.version(lj_property) !=
                         new_molecule.version(lj_property) );

    bool new_scl = ( old_molecule.version(scl_property) !=
                         new_molecule.version(scl_property) );

    if (new_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        coords_property) );

    if (new_clj)
        new_params.setAtomicParameters( getCLJParameters(new_molecule,
                                            chg_property, lj_property) );

    if (new_scl)
        new_params.setIntraScaleFactors( 
                IntraScaledParameters<CLJNBPairs>(new_molecule, scl_property) );

    return new_params;
}
                 
/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', also while the parameters of 'old_molecule'
    where found in 'old_map', now get the parameters using 'new_map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCLJPotential::Parameters
IntraCLJPotential::updateParameters(const IntraCLJPotential::Parameters &old_params,
                                    const PartialMolecule &old_molecule,
                                    const PartialMolecule &new_molecule,
                                    const PropertyMap &old_map, 
                                    const PropertyMap &new_map)
{
    if (old_molecule.selection() != new_molecule.selection())
        //the selection has changed - just get completely new parameters
        return this->getParameters(new_molecule, new_map);

    Parameters new_params = old_params;

    //get the property names
    const PropertyName &old_coords = old_map[parameters().coordinates()];
    const PropertyName &old_chg = old_map[parameters().charge()];
    const PropertyName &old_lj = old_map[parameters().lj()];
    const PropertyName &old_scl = old_map[parameters().intraScaleFactors()];
    
    const PropertyName &new_coords = new_map[parameters().coordinates()];
    const PropertyName &new_chg = new_map[parameters().charge()];
    const PropertyName &new_lj = new_map[parameters().lj()];
    const PropertyName &new_scl = new_map[parameters().intraScaleFactors()];
    
    //get what has changed
    bool changed_coords = (new_coords != old_coords) or
                           old_molecule.version(old_coords) !=
                           new_molecule.version(old_coords);
                             
    bool changed_clj = (new_chg != old_chg or new_lj != old_lj) or
                       ( old_molecule.version(old_chg) !=
                         new_molecule.version(old_chg) ) or
                       ( old_molecule.version(old_lj) !=
                         new_molecule.version(old_lj) );

    bool changed_scl = (new_scl != old_scl) or
                        old_molecule.version(old_scl) !=
                        new_molecule.version(old_scl);

    if (changed_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, new_coords) );

    if (changed_clj)
        new_params.setAtomicParameters( getCLJParameters(new_molecule,
                                                         new_chg, new_lj) );

    if (changed_scl)
        new_params.setIntraScaleFactors( 
                        IntraScaledParameters<CLJNBPairs>(new_molecule, new_scl) );

    return new_params;
}

/** Return the IntraCLJPotential::Molecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCLJPotential::Molecule
IntraCLJPotential::parameterise(const PartialMolecule &molecule,
                                const PropertyMap &map)
{
    return IntraCLJPotential::Molecule(molecule, *this, map);
}

/** Concert the passed group of molecules into IntraCLJPotential::Molecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCLJPotential::Molecules 
IntraCLJPotential::parameterise(const MolGroup &molecules,
                                const PropertyMap &map)
{
    return IntraCLJPotential::Molecules(molecules, *this, map);
}

/** Return the total charge of the parameters for the group in 'params' */
double IntraCLJPotential::totalCharge(
                        const IntraCLJPotential::Parameters::Array &params) const
{
    int n = params.count();
    const Parameter *params_array = params.constData();
    
    double chg = 0;
    
    for (int i=0; i<n; ++i)
    {
        chg += params_array[i].reduced_charge;
    }
    
    return chg;
}

void IntraCLJPotential::calculateEnergy(const CLJNBPairs::CGPairs &group_pairs, 
                            IntraCLJPotential::EnergyWorkspace &distmat,
                            const IntraCLJPotential::Parameter *params0_array, 
                            const IntraCLJPotential::Parameter *params1_array,
                            const quint32 nats0, const quint32 nats1, 
                            double &icnrg, double &iljnrg) const
{
    if (group_pairs.isEmpty())
    {
        //there are no scale factors between groups...
        for (quint32 i=0; i<nats0; ++i)
        {
            distmat.setOuterIndex(i);
            const Parameter &param0 = params0_array[i];
                
            if (param0.ljid == 0)
            {
                //null LJ parameter - only add on the coulomb energy
                for (quint32 j=0; j<nats1; ++j)
                {
                    icnrg += param0.reduced_charge * 
                             params1_array[j].reduced_charge / distmat[j];
                }
            }
            else
            {
                for (quint32 j=0; j<nats1; ++j)
                {
                    //do both coulomb and LJ
                    const Parameter &param1 = params1_array[j];
                        
                    const double invdist = double(1) / distmat[j];
                        
                    icnrg += param0.reduced_charge * param1.reduced_charge
                                                   * invdist;
                              
                    if (param1.ljid != 0)
                    {
                        const LJPair &ljpair = ljpairs.constData()[
                                                  ljpairs.map(param0.ljid,
                                                              param1.ljid)];
                        
                        double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                        double sig_over_dist12 = pow_2(sig_over_dist6);

                        iljnrg += ljpair.epsilon() * (sig_over_dist12 - 
                                                      sig_over_dist6);
                    }
                }
            }
        } 
    }
    else
    {
        //there are different nb scale factors between
        //the atoms. We need to calculate the energies using
        //them...
        for (quint32 i=0; i<nats0; ++i)
        {
            distmat.setOuterIndex(i);
            const Parameter &param0 = params0_array[i];
                
            if (param0.ljid == 0)
            {
                //null LJ parameter - only add on the coulomb energy
                for (quint32 j=0; j<nats1; ++j)
                {
                    const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                    if (cljscl.coulomb() != 0)
                           icnrg += cljscl.coulomb() * 
                                    param0.reduced_charge * 
                                    params1_array[j].reduced_charge / distmat[j];
                }
            }
            else
            {
                for (quint32 j=0; j<nats1; ++j)
                {
                    //do both coulomb and LJ
                    const CLJScaleFactor &cljscl = group_pairs(i,j);

                    if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                    {
                        const Parameter &param1 = params1_array[j];
                        
                        const double invdist = double(1) / distmat[j];
                        
                        icnrg += cljscl.coulomb() *  
                                 param0.reduced_charge * 
                                 param1.reduced_charge * invdist;
                              
                        if (cljscl.lj() != 0 and param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                     ljpairs.map(param0.ljid,
                                                                 param1.ljid)];
                        
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            iljnrg += cljscl.lj() * ljpair.epsilon() * 
                                       (sig_over_dist12 - sig_over_dist6);
                        }
                    }
                }
            }
        }
    }
}

/** Calculate the intramolecular CLJ energy of the passed molecule, and
    add this onto 'energy'. This uses the passed workspace when
    performing the calculation */
void IntraCLJPotential::calculateEnergy(const IntraCLJPotential::Molecule &mol,
                                        IntraCLJPotential::Energy &energy,
                                        IntraCLJPotential::EnergyWorkspace &distmat,
                                        double scale_energy) const
{
    if (scale_energy == 0 or mol.isEmpty())
        return;

    const quint32 ngroups = mol.nCutGroups();
    
    const CoordGroup *groups_array = mol.coordinates().constData();
    
    BOOST_ASSERT( mol.parameters().atomicParameters().count() == int(ngroups) );
    const Parameters::Array *molparams_array 
                            = mol.parameters().atomicParameters().constData();

    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    double cnrg = 0;
    double ljnrg = 0;
      
    //loop over all pairs of CutGroups in the molecule
    for (quint32 igroup=0; igroup<ngroups; ++igroup)
    {
        const Parameters::Array &params0 = molparams_array[igroup];

        const CoordGroup &group0 = groups_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
    
        for (quint32 jgroup=igroup; jgroup<ngroups; ++jgroup)
        {
            const CoordGroup &group1 = groups_array[jgroup];
            const Parameters::Array &params1 = molparams_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(don't test igroup==jgroup as this is the same CutGroup
            // and definitely within cutoff!)
            const bool within_cutoff = (igroup == jgroup) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDist(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;
                
            CGIdx cgidx_jgroup = mol.cgIdx(jgroup);
                
            //get the non-bonded scale factors for all pairs of atoms
            //between these groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            double icnrg = 0;
            double iljnrg = 0;
            
            //loop over all intraatomic pairs and calculate the energies
            const quint32 nats1 = group1.count();
            const Parameter *params1_array = params1.constData();
            
            calculateEnergy(group_pairs, distmat, params0_array, params1_array,
                            nats0, nats1, icnrg, iljnrg);
            
            //if this is the same group then half the energies to 
            //correct for double-counting
            if (igroup == jgroup)
            {
                icnrg *= 0.5;
                iljnrg *= 0.5;
            }

            //are we shifting the electrostatic potential?
            if (use_electrostatic_shifting and igroup != jgroup)
                icnrg -= this->totalCharge(params0) * this->totalCharge(params1)
                              / switchfunc->electrostaticCutoffDistance();

            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            if (mindist > switchfunc->featherDistance())
            {
                cnrg += switchfunc->electrostaticScaleFactor(mindist) * icnrg;
                ljnrg += switchfunc->vdwScaleFactor(mindist) * iljnrg;
            }
            else
            {
                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
    }
    
    //add this molecule pair's energy onto the total
    energy += Energy(scale_energy * cnrg, scale_energy * ljnrg);
}

static bool intersects(const QSet<Index> &atoms0, const QSet<Index> &atoms1)
{
    if (atoms0.count() <= atoms1.count())
    {
        for (QSet<Index>::const_iterator it = atoms0.constBegin();
             it != atoms0.constEnd();
             ++it)
        {
            if (atoms1.contains(*it))
                return true;
        }
    }
    else
    {
        for (QSet<Index>::const_iterator it = atoms1.constBegin();
             it != atoms1.constEnd();
             ++it)
        {
            if (atoms0.contains(*it))
                return true;
        }
    }
    
    return false;
}

/** Calculate the intramolecular CLJ energy of the passed molecule
    interacting with the rest of the same molecule in 'rest_of_mol', and
    add this onto 'energy'. This uses the passed workspace when
    performing the calculation. Note that mol and rest_of_mol should
    not contain any overlapping atoms, and that they should both be
    part of the same molecule (albeit potentially at different versions,
    but with the same layout UID)
    
    \throw SireError::incompatible_error
*/
void IntraCLJPotential::calculateEnergy(const IntraCLJPotential::Molecule &mol,
                                        const IntraCLJPotential::Molecule &rest_of_mol,
                                        IntraCLJPotential::Energy &energy,
                                        IntraCLJPotential::EnergyWorkspace &distmat,
                                        double scale_energy) const
{
    if (scale_energy == 0 or mol.isEmpty() or rest_of_mol.isEmpty())
        return;

    //ensure that this is the same molecule, with the same layout UID
    this->assertCompatible(mol, rest_of_mol);

    const quint32 ngroups0 = mol.nCutGroups();
    const CoordGroup *groups0_array = mol.coordinates().constData();
    
    BOOST_ASSERT( mol.parameters().atomicParameters().count() == int(ngroups0) );
    const Parameters::Array *molparams0_array 
                            = mol.parameters().atomicParameters().constData();

    const quint32 ngroups1 = rest_of_mol.nCutGroups();
    const CoordGroup *groups1_array = rest_of_mol.coordinates().constData();
    
    BOOST_ASSERT( rest_of_mol.parameters().atomicParameters().count() == int(ngroups1) );
    const Parameters::Array *molparams1_array
                            = rest_of_mol.parameters().atomicParameters().constData();

    //the CLJNBPairs must be the same in both molecules - this is checked
    //as part of assertCompatible(..)
    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    double cnrg = 0;
    double ljnrg = 0;
    
    //calculate the energy of all of the atoms in 'mol' interacting with
    //all of the atoms in 'rest_of_mol'
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
    
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            CGIdx cgidx_jgroup = rest_of_mol.cgIdx(jgroup);

            //check first that these two CoordGroups could be within cutoff
            //(don't test igroup==jgroup as this is the same CutGroup
            // and definitely within cutoff!)
            const bool within_cutoff = (cgidx_igroup == cgidx_jgroup) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDist(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;
                
            //get the non-bonded scale factors for all pairs of atoms
            //between these groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            double icnrg = 0;
            double iljnrg = 0;
            
            //loop over all intraatomic pairs and calculate the energies
            const quint32 nats1 = group1.count();
            const Parameter *params1_array = params1.constData();
            
            if (cgidx_igroup == cgidx_jgroup)
            {
                //this is the same CutGroup - we must be careful not to 
                //calculate the energy between overlapping atoms
                
                //get the atoms from this CutGroup that are contained in 
                //each part of the molecule
                QSet<Index> atoms0 = mol.molecule()
                                        .selection().selectedAtoms(cgidx_igroup);
                                        
                QSet<Index> atoms1 = rest_of_mol.molecule()
                                        .selection().selectedAtoms(cgidx_jgroup);
                                        
                BOOST_ASSERT( not ::intersects(atoms0, atoms1) );
                
                foreach (Index i, atoms0)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];

                    if (param0.ljid == 0)
                    {
                        //null LJ parameter - only add on the coulomb energy
                        foreach (Index j, atoms1)
                        {
                            const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                            if (cljscl.coulomb() != 0)
                                icnrg += cljscl.coulomb() * 
                                         param0.reduced_charge * 
                                         params1_array[j].reduced_charge / distmat[j];
                        }
                    }
                    else
                    {
                        foreach (Index j, atoms1)
                        {
                            //do both coulomb and LJ
                            const CLJScaleFactor &cljscl = group_pairs(i,j);

                            if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                            {
                                const Parameter &param1 = params1_array[j];
                        
                                const double invdist = double(1) / distmat[j];
                        
                                icnrg += cljscl.coulomb() *  
                                        param0.reduced_charge * 
                                        param1.reduced_charge * invdist;
                              
                                if (cljscl.lj() != 0 and param1.ljid != 0)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param1.ljid)];
                        
                                    double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);

                                    iljnrg += cljscl.lj() * ljpair.epsilon() * 
                                                (sig_over_dist12 - sig_over_dist6);
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                calculateEnergy(group_pairs, distmat,
                                params0_array, params1_array,
                                nats0, nats1, icnrg, iljnrg);
            }

            //are we shifting the electrostatic potential?
            if (use_electrostatic_shifting and cgidx_igroup != cgidx_jgroup)
                icnrg -= this->totalCharge(params0) * this->totalCharge(params1)
                              / switchfunc->electrostaticCutoffDistance();

            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            if (mindist > switchfunc->featherDistance())
            {
                cnrg += switchfunc->electrostaticScaleFactor(mindist) * icnrg;
                ljnrg += switchfunc->vdwScaleFactor(mindist) * iljnrg;
            }
            else
            {
                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
    }
}

void IntraCLJPotential::calculateForce(const CLJNBPairs::CGPairs &group_pairs,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraCLJPotential::ForceWorkspace &distmat,
                             const IntraCLJPotential::Parameter *params0_array,
                             const IntraCLJPotential::Parameter *params1_array,
                             const quint32 nats0, const quint32 nats1,
                             const double shift_coul,
                             Vector *group_forces0_array,
                             const double scale_force) const
{
    if (mindist > switchfunc->featherDistance())
    {
        //we need to calculate the forces taking into account
        //the derivative of the switching function!
                
        //calculate the switching scale factors and their 
        //derivatives
        const double scl_coul = switchfunc->electrostaticScaleFactor(mindist);
        const double scl_lj = switchfunc->vdwScaleFactor(mindist);
                
        Vector group_sep = (group1.aaBox().center() -
                            group0.aaBox().center()).normalise(); 
                
        Vector dscl_coul = switchfunc->dElectrostaticScaleFactor(mindist) 
                                     * group_sep;
                                     
        Vector dscl_lj = switchfunc->dVDWScaleFactor(mindist)
                                     * group_sep;

        if (group_pairs.isEmpty())
        {
            //there are no scale factors between atoms in these groups
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];

                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double q2 = param0.reduced_charge *
                                          params1_array[j].reduced_charge;
                                
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = scl_coul * q2 /
                                                distmat[j].length();

                            //calculate the coulomb force
                            Vector cforce = (-cnrg / distmat[j].length() *
                                                     distmat[j].direction()) +
                                             
                                                    ((cnrg-shift_coul) * dscl_coul);

                            total_force -= cforce;
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                        
                        const double invdist = double(1) / distmat[j].length();
                        
                        Vector force;

                        const double q2 = param0.reduced_charge *
                                          param1.reduced_charge;
                                
                        if (q2 != 0)
                        {
                            //calculate the energy
                            const double cnrg = scl_coul * q2 * invdist;

                            //calculate the force
                            force = (-cnrg / distmat[j].length() *
                                             distmat[j].direction()) +
                                             
                                             ((cnrg-shift_coul) * dscl_coul);
                        }
                              
                        if (param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                        
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = ljpair.epsilon() *
                                                   (sig_over_dist12 - sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            force += ((scl_lj * ljpair.epsilon() * 
                                       (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                           * distmat[j].direction())
                                            
                                     + (ljnrg * dscl_lj);
                        }
                        
                        total_force -= force;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //the atoms. We need to calculate the forces
            //using them...
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                        if (cljscl.coulomb() != 0)
                        {
                            const double q2 = param0.reduced_charge * 
                                              params1_array[j].reduced_charge;
                                                      
                            if (q2 != 0)
                            {
                                //calculate the coulomb energy
                                const double cnrg = cljscl.coulomb() *
                                                      scl_coul * q2
                                                      / distmat[j].length();
                                               
                                //calculate the coulomb force
                                Vector cforce = (-cnrg 
                                                  / distmat[j].length() *
                                                    distmat[j].direction()) +
                                             
                                                  ((cnrg-shift_coul) * dscl_coul);
                            
                                total_force -= cforce;
                            }
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];
                        
                            const double invdist = double(1) 
                                                     / distmat[j].length();
                        
                            Vector force;
                                    
                            const double q2 = param0.reduced_charge *
                                              param1.reduced_charge;

                            if (q2 != 0)
                            {
                                //calculate the energy
                                const double cnrg = cljscl.coulomb() *
                                                            scl_coul * q2 * invdist;
                        
                                //calculate the force
                                force = (scl_coul * -cnrg 
                                           / distmat[j].length() *
                                             distmat[j].direction()) +
                                             
                                        ((cnrg-shift_coul) * dscl_coul);
                            }
                            
                            if (cljscl.lj() != 0 and param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                            
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                //calculate the energy
                                const double ljnrg = cljscl.lj() *
                                                      ljpair.epsilon() *
                                                   (sig_over_dist12 - sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                force += ((scl_lj * ljpair.epsilon() * 
                                         (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                          * distmat[j].direction())
                                            
                                          + (ljnrg * dscl_lj);
                            }
                        
                            total_force -= force;
                        }
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
    }
    else
    {
        //not in the feather region, so can calculate the forces
        //directly (also, no need to calculate shift, as 
        //the shifting function is constant, so does not
        //affect the gradient)
                
        if (group_pairs.isEmpty())
        {
            //no nb scale factors to worry about
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //calculate the coulomb force
                        Vector cforce = -(param0.reduced_charge *
                                          params1_array[j].reduced_charge / 
                                          distmat[j].length2()) *
                                             
                                          distmat[j].direction();

                        total_force -= cforce;
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                            
                        const double invdist = double(1) / distmat[j].length();
                        const double invdist2 = pow_2(invdist);
                        
                        //calculate the force
                        Vector force = -(param0.reduced_charge * 
                                         param1.reduced_charge * invdist2) 
                                            
                                         * distmat[j].direction();
                              
                        if (param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                        
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            force += (ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                         12.0*sig_over_dist12))
                                      * distmat[j].direction();
                        }

                        total_force -= force;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //different atoms...
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                        if (cljscl.coulomb() != 0)
                        {
                            //calculate the coulomb force
                            Vector cforce = -(cljscl.coulomb() *
                                              param0.reduced_charge *
                                              params1_array[j].reduced_charge / 
                                              distmat[j].length2()) *
                                             
                                              distmat[j].direction();

                            total_force -= cforce;
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];
                            
                            const double invdist = double(1)    
                                                    / distmat[j].length();
                            const double invdist2 = pow_2(invdist);
                        
                            //calculate the force
                            Vector force = -(cljscl.coulomb() *
                                             param0.reduced_charge * 
                                             param1.reduced_charge * invdist2) 
                                            
                                             * distmat[j].direction();
                              
                            if (cljscl.lj() != 0 and param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                         ljpairs.map(param0.ljid,
                                                                     param1.ljid)];
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                force += (cljscl.lj() *
                                          ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                             12.0*sig_over_dist12))
                                           * distmat[j].direction();
                            }

                            total_force -= force;
                        }
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;

            } // end of loop over i atoms

        } // end of whether there are intra scale factors

    } // end of whether within feather region
}

/** Calculate the coulomb and LJ forces between the atoms in the molecule 'mol'
    and add these forces onto 'forces'. This uses
    the passed workspace to perform the calculation */
void IntraCLJPotential::calculateForce(const IntraCLJPotential::Molecule &mol,
                                       MolForceTable &forces, 
                                       IntraCLJPotential::ForceWorkspace &distmat,
                                       double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty())
        return;
    
    const quint32 ngroups = mol.nCutGroups();
    
    const CoordGroup *groups_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams_array 
                            = mol.parameters().atomicParameters().constData();
    
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();

    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    //loop over all pairs of CutGroups in the molecule
    for (quint32 igroup=0; igroup<ngroups; ++igroup)
    {
        const Parameters::Array &params0 = molparams_array[igroup];

        const CoordGroup &group0 = groups_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
    
        //get the forces acting on this group from all other groups
        //(yes, we are doing a full n2 loop, and not taking advantages
        // of equal and opposite forces)
        for (quint32 jgroup=0; jgroup<ngroups; ++jgroup)
        {
            const CoordGroup &group1 = groups_array[jgroup];
            const Parameters::Array &params1 = molparams_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = igroup != jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            CGIdx cgidx_jgroup = mol.cgIdx(jgroup);

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            double shift_coul = 0;

            if (use_electrostatic_shifting and igroup != jgroup)
                shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                              / switchfunc->electrostaticCutoffDistance();

            //now calculate the forces acting on group0 caused by group1
            calculateForce(group_pairs, group0, group1,
                           mindist, distmat, params0_array, params1_array,
                           nats0, nats1, shift_coul, group_forces0_array,
                           scale_force);
             
        } // end of loop over CutGroups (jgroup)

    } // end of loop over CutGroups (igroup)
}

/** Calculate the total forces acting on the atoms in 'mol' caused by the 
    other atoms in the same molecule contained in 'rest_of_mol'. This calculates
    the forces and adds them onto 'forces' (which are for 'mol'). Note that
    'mol' and 'rest_of_mol' must not overlap, and also that they must
    use the same layout UID and same intra-nonbonded scaling factors
    
    \throw SireError::incompatible_error
*/
void IntraCLJPotential::calculateForce(const IntraCLJPotential::Molecule &mol,
                                       const IntraCLJPotential::Molecule &rest_of_mol,
                                       MolForceTable &forces,
                                       IntraCLJPotential::ForceWorkspace &distmat,
                                       double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty() or rest_of_mol.isEmpty())
        return;
        
    const quint32 ngroups0 = mol.nCutGroups();
    const CoordGroup *groups0_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams0_array
                            = mol.parameters().atomicParameters().constData();
                            
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();
    
    this->assertCompatible(mol, rest_of_mol);
    
    const quint32 ngroups1 = rest_of_mol.nCutGroups();
    const CoordGroup *groups1_array = rest_of_mol.coordinates().constData();
    
    const Parameters::Array *molparams1_array
                            = rest_of_mol.parameters().atomicParameters().constData();

    //this is now guaranteed to be the same in both passed parts of 
    //the molecule
    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();

    //calculate the forces acting on the atoms in 'mol' caused by 
    //the atoms in 'rest_of_mol'
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
        
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            CGIdx cgidx_jgroup = rest_of_mol.cgIdx(jgroup);

            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = cgidx_igroup != cgidx_jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);
            
            if (cgidx_igroup == cgidx_jgroup)
            {
                //we have to be careful to not calculate the forces of
                //overlapping atoms...
                
                QSet<Index> atoms0 = mol.molecule().selection()
                                               .selectedAtoms(cgidx_igroup);
                                               
                QSet<Index> atoms1 = rest_of_mol.molecule().selection()
                                               .selectedAtoms(cgidx_jgroup);
                                         
                BOOST_ASSERT( not ::intersects(atoms0, atoms1) );
                                                           
                foreach (Index i, atoms0)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];
                
                    Vector total_force;
                
                    if (param0.ljid == 0)
                    {
                        //null LJ parameter - only add on the coulomb energy
                        foreach (Index j, atoms1)
                        {
                            const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                            if (cljscl.coulomb() != 0)
                            {
                                //calculate the coulomb force
                                Vector cforce = -(cljscl.coulomb() *
                                                  param0.reduced_charge *
                                                  params1_array[j].reduced_charge / 
                                                  distmat[j].length2()) *
                                             
                                                  distmat[j].direction();

                                total_force -= cforce;
                            }
                        }
                    }
                    else
                    {
                        foreach (Index j, atoms1)
                        {
                            //do both coulomb and LJ
                            const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                            if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                            {
                                const Parameter &param1 = params1_array[j];
                            
                                const double invdist = double(1)    
                                                        / distmat[j].length();
                                const double invdist2 = pow_2(invdist);
                        
                                //calculate the force
                                Vector force = -(cljscl.coulomb() *
                                                 param0.reduced_charge * 
                                                 param1.reduced_charge * invdist2) 
                                            
                                                 * distmat[j].direction();
                              
                                if (cljscl.lj() != 0 and param1.ljid != 0)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param1.ljid)];
                        
                                    double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);

                                    // dU/dr requires an extra power of r
                                    sig_over_dist6 *= invdist;
                                    sig_over_dist12 *= invdist;

                                    force += (cljscl.lj() *
                                            ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                                12.0*sig_over_dist12))
                                            * distmat[j].direction();
                                }

                                total_force -= force;
                            }
                        }
                    }
                        
                    group_forces0_array[i] += scale_force * total_force;
                }
            }
            else
            {
                double shift_coul = 0;
            
                if (use_electrostatic_shifting)
                    shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                                / switchfunc->electrostaticCutoffDistance();

                calculateForce(group_pairs, group0, group1,
                               mindist, distmat,
                               params0_array, params1_array,
                               nats0, nats1, shift_coul,
                               group_forces0_array, scale_force);
            }
        }
    }
}

void IntraCLJPotential::calculateCoulombForce(const CLJNBPairs::CGPairs &group_pairs,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraCLJPotential::ForceWorkspace &distmat,
                             const IntraCLJPotential::Parameter *params0_array,
                             const IntraCLJPotential::Parameter *params1_array,
                             const quint32 nats0, const quint32 nats1,
                             const double shift_coul,
                             Vector *group_forces0_array,
                             const double scale_force) const
{
    if (mindist > switchfunc->featherDistance())
    {
        //we need to calculate the forces taking into account
        //the derivative of the switching function!
                
        //calculate the switching scale factors and their 
        //derivatives
        const double scl_coul = switchfunc->electrostaticScaleFactor(mindist);
                
        Vector group_sep = (group1.aaBox().center() -
                            group0.aaBox().center()).normalise(); 
                
        Vector dscl_coul = switchfunc->dElectrostaticScaleFactor(mindist) 
                              * group_sep;

        if (group_pairs.isEmpty())
        {
            //there are no scale factors between atoms in these groups
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];

                if (param0.reduced_charge == 0)
                    return;

                distmat.setOuterIndex(i);

                Vector total_force;
                
                for (quint32 j=0; j<nats1; ++j)
                {
                    const double q2 = param0.reduced_charge *
                                      params1_array[j].reduced_charge;
                                
                    if (q2 != 0)
                    {
                        //calculate the coulomb energy
                        const double cnrg = scl_coul * q2 /
                                            distmat[j].length();

                        //calculate the coulomb force
                        Vector cforce = (-cnrg / distmat[j].length() *
                                                 distmat[j].direction()) +
                                             
                                                ((cnrg-shift_coul) * dscl_coul);

                        total_force -= cforce;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //the atoms. We need to calculate the forces
            //using them...
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];

                if (param0.reduced_charge == 0)
                    continue;
                
                distmat.setOuterIndex(i);

                Vector total_force;
                
                for (quint32 j=0; j<nats1; ++j)
                {
                    const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                    if (cljscl.coulomb() != 0)
                    {
                        const double q2 = param0.reduced_charge * 
                                          params1_array[j].reduced_charge;
                                                      
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = cljscl.coulomb() *
                                                  scl_coul * q2
                                                  / distmat[j].length();
                                               
                            //calculate the coulomb force
                            Vector cforce = (-cnrg 
                                             / distmat[j].length() *
                                               distmat[j].direction()) +
                                             
                                              ((cnrg-shift_coul) * dscl_coul);
                            
                            total_force -= cforce;
                        }
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
    }
    else
    {
        //not in the feather region, so can calculate the forces
        //directly (also, no need to calculate shift, as 
        //the shifting function is constant, so does not
        //affect the gradient)
                
        if (group_pairs.isEmpty())
        {
            //no nb scale factors to worry about
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];

                if (param0.reduced_charge == 0)
                    continue;
                    
                distmat.setOuterIndex(i);
                
                Vector total_force;
                
                for (quint32 j=0; j<nats1; ++j)
                {
                    if (params1_array[j].reduced_charge != 0)
                    {
                        //calculate the coulomb force
                        Vector cforce = -(param0.reduced_charge *
                                          params1_array[j].reduced_charge / 
                                          distmat[j].length2()) *
                                             
                                          distmat[j].direction();

                        total_force -= cforce;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //different atoms...
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];
                
                if (param0.reduced_charge == 0)    
                    continue;
                
                distmat.setOuterIndex(i);

                Vector total_force;
                
                for (quint32 j=0; j<nats1; ++j)
                {
                    const CLJScaleFactor &cljscl = group_pairs(i,j);
                          
                    double cscale = cljscl.coulomb() 
                                      * params1_array[j].reduced_charge;
                           
                    if (cscale != 0)
                    {
                        //calculate the coulomb force
                        Vector cforce = -(cscale *
                                          param0.reduced_charge /
                                          distmat[j].length2()) *
                                             
                                          distmat[j].direction();

                        total_force -= cforce;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;

            } // end of loop over i atoms
        } // end of whether there are intra scale factors
    } // end of whether within feather region
}

/** Calculate the coulomb forces between the atoms in the molecule 'mol'
    and add these forces onto 'forces'. This uses
    the passed workspace to perform the calculation */
void IntraCLJPotential::calculateCoulombForce(
                                       const IntraCLJPotential::Molecule &mol,
                                       MolForceTable &forces, 
                                       IntraCLJPotential::ForceWorkspace &distmat,
                                       double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty())
        return;
    
    const quint32 ngroups = mol.nCutGroups();
    
    const CoordGroup *groups_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams_array 
                            = mol.parameters().atomicParameters().constData();
    
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();

    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    //loop over all pairs of CutGroups in the molecule
    for (quint32 igroup=0; igroup<ngroups; ++igroup)
    {
        const Parameters::Array &params0 = molparams_array[igroup];

        const CoordGroup &group0 = groups_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
    
        //get the forces acting on this group from all other groups
        //(yes, we are doing a full n2 loop, and not taking advantages
        // of equal and opposite forces)
        for (quint32 jgroup=0; jgroup<ngroups; ++jgroup)
        {
            const CoordGroup &group1 = groups_array[jgroup];
            const Parameters::Array &params1 = molparams_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = igroup != jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            CGIdx cgidx_jgroup = mol.cgIdx(jgroup);

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            double shift_coul = 0;

            if (use_electrostatic_shifting and igroup != jgroup)
               shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                           / switchfunc->electrostaticCutoffDistance();

            //calculate the forces acting on group0 caused by group1
            calculateCoulombForce(group_pairs, group0, group1,
                                  mindist, distmat, 
                                  params0_array, params1_array,
                                  nats0, nats1, shift_coul, 
                                  group_forces0_array, scale_force);

        } // end of loop over CutGroups (jgroup)

    } // end of loop over CutGroups (igroup)
}

/** Calculate the coulomb force acting on the part of the molecule
    in 'mol' caused by the rest of the molecule in 'rest_of_mol'. Note
    that these must both be of the same molecule, with the same
    layout UID and same nonbonded scale factors
    
    \throw SireError::incompatible_error
*/
void IntraCLJPotential::calculateCoulombForce(
                                        const IntraCLJPotential::Molecule &mol,
                                        const IntraCLJPotential::Molecule &rest_of_mol,
                                        MolForceTable &forces,
                                        IntraCLJPotential::ForceWorkspace &distmat,
                                        double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty() or rest_of_mol.isEmpty())
        return;
        
    const quint32 ngroups0 = mol.nCutGroups();
    const CoordGroup *groups0_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams0_array
                            = mol.parameters().atomicParameters().constData();
                            
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();
    
    this->assertCompatible(mol, rest_of_mol);
    
    const quint32 ngroups1 = rest_of_mol.nCutGroups();
    const CoordGroup *groups1_array = rest_of_mol.coordinates().constData();
    
    const Parameters::Array *molparams1_array
                            = rest_of_mol.parameters().atomicParameters().constData();

    //this is now guaranteed to be the same in both passed parts of 
    //the molecule
    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();

    //calculate the forces acting on the atoms in 'mol' caused by 
    //the atoms in 'rest_of_mol'
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
        
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            CGIdx cgidx_jgroup = rest_of_mol.cgIdx(jgroup);

            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = cgidx_igroup != cgidx_jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);
            
            if (cgidx_igroup == cgidx_jgroup)
            {
                //we have to be careful to not calculate the forces of
                //overlapping atoms...
                
                QSet<Index> atoms0 = mol.molecule().selection()
                                               .selectedAtoms(cgidx_igroup);
                                               
                QSet<Index> atoms1 = rest_of_mol.molecule().selection()
                                               .selectedAtoms(cgidx_jgroup);
                                         
                BOOST_ASSERT( not ::intersects(atoms0, atoms1) );
                                                           
                foreach (Index i, atoms0)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];
                
                    Vector total_force;
                
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                        if (cljscl.coulomb() != 0)
                        {
                            //calculate the coulomb force
                            Vector cforce = -(cljscl.coulomb() *
                                              param0.reduced_charge *
                                              params1_array[j].reduced_charge / 
                                              distmat[j].length2()) *
                                             
                                              distmat[j].direction();

                            total_force -= cforce;
                        }
                    }
                        
                    group_forces0_array[i] += scale_force * total_force;
                }
            }
            else
            {
                double shift_coul = 0;
            
                if (use_electrostatic_shifting)
                    shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                                / switchfunc->electrostaticCutoffDistance();

                calculateCoulombForce(group_pairs, group0, group1,
                                      mindist, distmat,
                                      params0_array, params1_array,
                                      nats0, nats1, shift_coul,
                                      group_forces0_array, scale_force);
            }
        }
    }
}

void IntraCLJPotential::calculateLJForce(const CLJNBPairs::CGPairs &group_pairs,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraCLJPotential::ForceWorkspace &distmat,
                             const IntraCLJPotential::Parameter *params0_array,
                             const IntraCLJPotential::Parameter *params1_array,
                             const quint32 nats0, const quint32 nats1,
                             Vector *group_forces0_array,
                             const double scale_force) const
{
    if (mindist > switchfunc->featherDistance())
    {
        //we need to calculate the forces taking into account
        //the derivative of the switching function!
        
        //calculate the switching scale factors and their 
        //derivatives
        const double scl_lj = switchfunc->vdwScaleFactor(mindist);
        
        Vector group_sep = (group1.aaBox().center() -
                            group0.aaBox().center()).normalise(); 
                             
        Vector dscl_lj = switchfunc->dVDWScaleFactor(mindist)
                             * group_sep;

        if (group_pairs.isEmpty())
        {
            //there are no scale factors between atoms in these groups
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];

                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    distmat.setOuterIndex(i);

                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                      
                        if (param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();
                
                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = ljpair.epsilon() *
                                              (sig_over_dist12 - sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force = 
                                    ((scl_lj * ljpair.epsilon() * 
                                    (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                    * distmat[j].direction())
                                    
                                    + (ljnrg * dscl_lj);

                            total_force -= force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //the atoms. We need to calculate the forces
            //using them...
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];
        
                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    distmat.setOuterIndex(i);

                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                        
                        if (cljscl.lj() != 0 and param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();
                            
                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = cljscl.lj() *
                                                 ljpair.epsilon() *
                                          (sig_over_dist12 - sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force = 
                                   ((scl_lj * ljpair.epsilon() * 
                                    (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                   * distmat[j].direction())
                                   
                                   + (ljnrg * dscl_lj);

                            total_force -= force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;
            }
        }
    }
    else
    {
        //not in the feather region, so can calculate the forces
        //directly
        if (group_pairs.isEmpty())
        {
            //no nb scale factors to worry about
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
        
                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                
                        if (param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();
                            
                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force =
                                (ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                     12.0*sig_over_dist12))
                                 * distmat[j].direction();

                            total_force -= force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //different atoms...
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];
        
                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    distmat.setOuterIndex(i);

                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                        
                        if (cljscl.lj() != 0 and param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();

                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force = 
                                  (cljscl.lj() *
                                   ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                       12.0*sig_over_dist12))
                                   * distmat[j].direction();

                            total_force -= force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;

            } // end of loop over i atoms

        } // end of whether there are intra scale factors

    } // end of whether within feather region
}

/** Calculate the LJ forces between the atoms in the molecule 'mol'
    and add these forces onto 'forces'. This uses
    the passed workspace to perform the calculation */
void IntraCLJPotential::calculateLJForce(const IntraCLJPotential::Molecule &mol,
                                         MolForceTable &forces, 
                                         IntraCLJPotential::ForceWorkspace &distmat,
                                         double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty())
        return;
    
    const quint32 ngroups = mol.nCutGroups();
    
    const CoordGroup *groups_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams_array 
                            = mol.parameters().atomicParameters().constData();
    
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();

    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    //loop over all pairs of CutGroups in the molecule
    for (quint32 igroup=0; igroup<ngroups; ++igroup)
    {
        const Parameters::Array &params0 = molparams_array[igroup];

        const CoordGroup &group0 = groups_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
    
        //get the forces acting on this group from all other groups
        //(yes, we are doing a full n2 loop, and not taking advantages
        // of equal and opposite forces)
        for (quint32 jgroup=0; jgroup<ngroups; ++jgroup)
        {
            const CoordGroup &group1 = groups_array[jgroup];
            const Parameters::Array &params1 = molparams_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = igroup != jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            CGIdx cgidx_jgroup = mol.cgIdx(jgroup);

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            //calculate the forces acting on group0 caused by group1
            calculateLJForce(group_pairs, group0, group1,
                             mindist, distmat, 
                             params0_array, params1_array,
                             nats0, nats1, group_forces0_array, scale_force);
            
        } // end of loop over CutGroups (jgroup)

    } // end of loop over CutGroups (igroup)
}

/** Calculate the LJ force acting on the part of the molecule
    in 'mol' caused by the rest of the molecule in 'rest_of_mol'. Note
    that these must both be of the same molecule, with the same
    layout UID and same nonbonded scale factors
    
    \throw SireError::incompatible_error
*/
void IntraCLJPotential::calculateLJForce(
                                        const IntraCLJPotential::Molecule &mol,
                                        const IntraCLJPotential::Molecule &rest_of_mol,
                                        MolForceTable &forces,
                                        IntraCLJPotential::ForceWorkspace &distmat,
                                        double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty() or rest_of_mol.isEmpty())
        return;
        
    const quint32 ngroups0 = mol.nCutGroups();
    const CoordGroup *groups0_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams0_array
                            = mol.parameters().atomicParameters().constData();
                            
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();
    
    this->assertCompatible(mol, rest_of_mol);
    
    const quint32 ngroups1 = rest_of_mol.nCutGroups();
    const CoordGroup *groups1_array = rest_of_mol.coordinates().constData();
    
    const Parameters::Array *molparams1_array
                            = rest_of_mol.parameters().atomicParameters().constData();

    //this is now guaranteed to be the same in both passed parts of 
    //the molecule
    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();

    //calculate the forces acting on the atoms in 'mol' caused by 
    //the atoms in 'rest_of_mol'
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
        
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            CGIdx cgidx_jgroup = rest_of_mol.cgIdx(jgroup);

            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = cgidx_igroup != cgidx_jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);
            
            if (cgidx_igroup == cgidx_jgroup)
            {
                //we have to be careful to not calculate the forces of
                //overlapping atoms...
                
                QSet<Index> atoms0 = mol.molecule().selection()
                                               .selectedAtoms(cgidx_igroup);
                                               
                QSet<Index> atoms1 = rest_of_mol.molecule().selection()
                                               .selectedAtoms(cgidx_jgroup);
                                         
                BOOST_ASSERT( not ::intersects(atoms0, atoms1) );
                                                           
                foreach (Index i, atoms0)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];
                
                    Vector total_force;
                
                    if (param0.ljid != 0)
                    {
                        foreach (Index j, atoms1)
                        {
                            const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                            const Parameter &param1 = params1_array[j];
                            
                            if (cljscl.lj() != 0 and param1.ljid != 0)
                            {
                                const double invdist = double(1) /
                                                  distmat(j).length();
                            
                                const LJPair &ljpair = ljpairs.constData()[
                                                         ljpairs.map(param0.ljid,
                                                                     param1.ljid)];
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                const Vector force = (cljscl.lj() *
                                          ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                             12.0*sig_over_dist12))
                                         * distmat[j].direction();

                                total_force -= force;
                            }
                        }
                        
                        group_forces0_array[i] += scale_force * total_force;
                    }
                }
            }
            else
            {
                calculateLJForce(group_pairs, group0, group1,
                                 mindist, distmat,
                                 params0_array, params1_array,
                                 nats0, nats1, group_forces0_array, scale_force);
            }
        }
    }
}
