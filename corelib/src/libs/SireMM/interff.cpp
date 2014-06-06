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

#include "interff.h"
#include "cljshiftfunction.h"
#include "cljcalculator.h"

#include "SireBase/booleanproperty.h"
#include "SireBase/lengthproperty.h"

#include "SireError/errors.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/molresid.h"
#include "SireMol/residue.h"
#include "SireMol/atomselection.h"
#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QElapsedTimer>
#include <QDebug>

using namespace SireMM;
using namespace SireMol;
using namespace SireFF;
using namespace SireBase;
using namespace SireStream;

namespace SireMM
{
    namespace detail
    {
        class InterFFData : public QSharedData
        {
        public:
            InterFFData() : QSharedData(), fixed_only(false)
            {}
            
            InterFFData(const InterFFData &other)
                 : maps_for_mol(other.maps_for_mol),
                   mol_to_idx(other.mol_to_idx),
                   molres_to_idx(other.molres_to_idx),
                   atom_selections(other.atom_selections),
                   cljfunc(other.cljfunc),
                   fixed_atoms(other.fixed_atoms),
                   cljcomps(other.cljcomps),
                   props(other.props),
                   fixed_only(other.fixed_only)
            {}
            
            ~InterFFData()
            {}
            
            /** The property map used to add each molecule */
            QHash<SireMol::MolNum,PropertyMap> maps_for_mol;

            /** The index of each molecule into the array. This is fixed
                throughout the lifetime of the forcefield */
            QHash<SireMol::MolNum,quint32> mol_to_idx;

            /** The index of each molecule/residue combination into the 
                array. This used for molecules that contain more than
                one residue, and is fixed throughout the lifetime of the forcefield */
            QHash<SireMol::MolResNum,quint32> molres_to_idx;

            /** The atom selections of each molecule/residue in the system. If the entry
                is empty, then this means that the whole molecule/residue is in the 
                forcefield */
            ChunkedVector<AtomSelection> atom_selections;

            /** The function used to calculate energies */
            CLJFunctionPtr cljfunc;
            
            /** All of the fixed atoms */
            CLJGrid fixed_atoms;

            /** The energy components available for this forcefield */
            CLJComponent cljcomps;
            
            /** All of the properties in this forcefield */
            Properties props;
            
            /** Whether or not to only calculate the energy with
                the fixed atoms */
            bool fixed_only;
        };
    }
}

static RegisterMetaType<InterFF> r_interff;

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const InterFF &interff)
{
    writeHeader(ds, r_interff, 1);
    
    SharedDataStream sds(ds);
    sds << static_cast<const G1FF&>(interff);
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, InterFF &interff)
{
    VersionID v = readHeader(ds, r_interff);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> static_cast<G1FF&>(interff);
        
        interff._pvt_updateName();
    }
    else
        throw version_error(v, "1", r_interff, CODELOC);
    
    return ds;
}

/** Constructor */
InterFF::InterFF() : ConcreteProperty<InterFF,G1FF>(), needs_reboxing(false)
{
    d = new detail::InterFFData();
    this->_pvt_updateName();
    this->setCLJFunction( CLJShiftFunction::defaultShiftFunction() );
}

/** Construct, specifying the name of the forcefield */
InterFF::InterFF(const QString &name) : ConcreteProperty<InterFF, G1FF>(), needs_reboxing(false)
{
    d = new detail::InterFFData();
    G1FF::setName(name);
    this->setCLJFunction( CLJShiftFunction::defaultShiftFunction() );
}

/** Copy constructor */
InterFF::InterFF(const InterFF &other)
        : ConcreteProperty<InterFF,G1FF>(other),
          atom_locs(other.atom_locs), wspace(other.wspace),
          cljboxes(other.cljboxes), d(other.d),
          needs_reboxing(other.needs_reboxing)
{}

/** Destructor */
InterFF::~InterFF()
{}

/** Function used to set the CLJFunction used to calculate the energy */
void InterFF::setCLJFunction(const CLJFunction &func)
{
    if (not d.constData()->cljfunc.read().equals(func))
    {
        d->fixed_atoms.setCLJFunction(func);
        d->cljfunc = func;
        rebuildProps();
        this->mustNowRecalculateFromScratch();
    }
}

/** Return the function used to calculate the energy */
const CLJFunction& InterFF::cljFunction() const
{
    return d.constData()->cljfunc.read();
}

/** Internal function called when the name of the forcefield changes */
void InterFF::_pvt_updateName()
{
    d->cljcomps = CLJComponent( this->name() );
    G1FF::_pvt_updateName();
}

const char* InterFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<InterFF>() );
}

const char* InterFF::what() const
{
    return InterFF::typeName();
}

/** Copy assignment operator */
InterFF& InterFF::operator=(const InterFF &other)
{
    //qDebug() << "COPY ASSIGNMENT";

    if (this != &other)
    {
        atom_locs = other.atom_locs;
        wspace = other.wspace;
        cljboxes = other.cljboxes;
        d = other.d;
        needs_reboxing = other.needs_reboxing;
        G1FF::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool InterFF::operator==(const InterFF &other) const
{
    return (this == &other) or
           (G1FF::operator==(other) and d->fixed_atoms == other.d->fixed_atoms);
}

/** Comparison operator */
bool InterFF::operator!=(const InterFF &other) const
{
    return not operator==(other);
}

InterFF* InterFF::clone() const
{
    return new InterFF(*this);
}

/** Return the energy components of this forcefield */
const CLJComponent& InterFF::components() const
{
    return d->cljcomps;
}

/** Internal function used to rebuild the properties object that 
    stores all of the properties of this forcefield */
void InterFF::rebuildProps()
{
    d->props = d->cljfunc.read().properties();
    d->props.setProperty("cljFunction", d->cljfunc);
    d->props.setProperty("useGrid", BooleanProperty(d->fixed_atoms.usesGrid()));
    d->props.setProperty("gridBuffer", LengthProperty(d->fixed_atoms.gridBuffer()));
    d->props.setProperty("gridSpacing", LengthProperty(d->fixed_atoms.gridSpacing()));
    d->props.setProperty("fixedOnly", BooleanProperty(d->fixed_only));
}

/** Set the forcefield property called 'name' to the value 'property' */
bool InterFF::setProperty(const QString &name, const Property &property)
{
    if (name == "cljFunction")
    {
        if (not d.constData()->cljfunc.read().equals(property))
        {
            this->setCLJFunction(property.asA<CLJFunction>());
            return true;
        }
        else
            return false;
    }
    else if (name == "useGrid")
    {
        bool use_grid = property.asA<BooleanProperty>().value();
        
        if (use_grid != d.constData()->fixed_atoms.usesGrid())
        {
            d->fixed_atoms.setUseGrid(use_grid);
            d->props.setProperty("useGrid", BooleanProperty(use_grid));
            return true;
        }
        else
            return false;
    }
    else if (name == "gridBuffer")
    {
        Length buffer = property.asA<LengthProperty>().value();
        
        if (buffer != d.constData()->fixed_atoms.gridBuffer())
        {
            d->fixed_atoms.setGridBuffer(buffer);
            d->props.setProperty("gridBuffer", property);
            return true;
        }
        else
            return false;
    }
    else if (name == "gridSpacing")
    {
        Length spacing = property.asA<LengthProperty>().value();
        
        if (spacing != d.constData()->fixed_atoms.gridSpacing())
        {
            d->fixed_atoms.setGridSpacing(spacing);
            d->props.setProperty("gridSpacing", property);
            return true;
        }
        else
            return false;
    }
    else if (name == "fixedOnly")
    {
        bool fixed_only = property.asA<BooleanProperty>().value();
        
        if (fixed_only != d.constData()->fixed_only)
        {
            d->fixed_only = fixed_only;
            d->props.setProperty("fixedOnly", property);
            this->mustNowRecalculateFromScratch();
            return true;
        }
        else
            return false;
    }
    else
    {
        PropertyPtr old_prop = d.constData()->cljfunc.read().property(name);
    
        if (not property.equals(old_prop.read()))
        {
            //need to set the property
            this->setCLJFunction( d.constData()->cljfunc.read().setProperty(name, property) );
            return true;
        }
        else
            return false;
    }
}

/** Return the value of the forcefield property with name 'name' */
const Property& InterFF::property(const QString &name) const
{
    return d->props.property(name);
}

/** Return whether or not this forcefield contains the property 'property' */
bool InterFF::containsProperty(const QString &name) const
{
    return d->props.hasProperty(name);
}

/** Return all of the properties of this function */
const Properties& InterFF::properties() const
{
    return d->props;
}

/** Add the passed atoms as fixed atoms to the forcefield */
void InterFF::addFixedAtoms(const CLJAtoms &atoms)
{
    d->fixed_atoms.addFixedAtoms(atoms);
    this->mustNowRecalculateFromScratch();
}

/** Add the passed atoms as fixed atoms to the forcefield */
void InterFF::addFixedAtoms(const MoleculeView &molecule, const PropertyMap &map)
{
    this->addFixedAtoms( CLJAtoms(molecule,map) );
}

/** Add the passed molecules as fixed atoms to the forcefield */
void InterFF::addFixedAtoms(const Molecules &molecules, const PropertyMap &map)
{
    this->addFixedAtoms( CLJAtoms(molecules,map) );
}

/** Set the fixed atoms equal to 'atoms' */
void InterFF::setFixedAtoms(const CLJAtoms &atoms)
{
    d->fixed_atoms.setFixedAtoms(atoms);
    this->mustNowRecalculateFromScratch();
}

/** Set the fixed atoms equal to 'molecule' */
void InterFF::setFixedAtoms(const MoleculeView &molecule, const PropertyMap &map)
{
    this->setFixedAtoms( CLJAtoms(molecule,map) );
}

/** Set the fixed atoms equal to 'molecules' */
void InterFF::setFixedAtoms(const Molecules &molecules, const PropertyMap &map)
{
    this->setFixedAtoms( CLJAtoms(molecules,map) );
}

/** Set whether or not the energy calculation is only between the mobile and 
    fixed atoms (i.e. the mobile-mobile interaction is ignored) */
void InterFF::setFixedOnly(bool on)
{
    if (d.constData()->fixed_only != on)
    {
        d->fixed_only = on;
        d->props.setProperty("fixedOnly", BooleanProperty(on));
        this->mustNowRecalculateFromScratch();
    }
}

/** Set whether or not a grid is used to optimise energy calculations with the fixed atoms */
void InterFF::setUseGrid(bool on)
{
    if (this->usesGrid() != on)
    {
        d->fixed_atoms.setUseGrid(on);
        
        if (this->usesGrid() == on)
        {
            this->mustNowRecalculateFromScratch();
            d->props.setProperty("useGrid", BooleanProperty(on));
        }
        else
        {
            //enabling the grid may not actually enable it (cljfunction may not support the grid)
            qDebug() << "Switching on the grid failed as the function"
                     << d->fixed_atoms.cljFunction().toString() << "doesn't appear "
                     << "to support using the grid.";
        }
    }
}

/** Turn on the use of the grid */
void InterFF::enableGrid()
{
    this->setUseGrid(true);
}

/** Turn off use of the grid */
void InterFF::disableGrid()
{
    this->setUseGrid(false);
}

/** Return whether or not the grid is used */
bool InterFF::usesGrid() const
{
    return d->fixed_atoms.usesGrid();
}

/** Return whether or not only the energy between the mobile and fixed
    atoms is being calculated */
bool InterFF::fixedOnly() const
{
    return d->fixed_only;
}

/** Set the buffer used when using a grid. This is the distance
    added around the maximum extent of the atoms when working out the
    dimension of the grid */
void InterFF::setGridBuffer(Length buffer)
{
    if (d.constData()->fixed_atoms.gridBuffer() != buffer)
    {
        d->fixed_atoms.setGridBuffer(buffer);
        d->props.setProperty("gridBuffer", LengthProperty(buffer));
        
        if (usesGrid())
            this->mustNowRecalculateFromScratch();
    }
}

/** Return the buffer used when working out the dimension of the grid */
Length InterFF::gridBuffer() const
{
    return d->fixed_atoms.gridBuffer();
}

/** Set the spacing between grid points */
void InterFF::setGridSpacing(Length spacing)
{
    if (d.constData()->fixed_atoms.gridSpacing() != spacing)
    {
        d->fixed_atoms.setGridSpacing(spacing);
        d->props.setProperty("gridSpacing", LengthProperty(spacing));
        
        if (usesGrid())
            this->mustNowRecalculateFromScratch();
    }
}

/** Return spacing between grid points */
Length InterFF::gridSpacing() const
{
    return d->fixed_atoms.gridSpacing();
}

/** Return the grid used to calculate the energy with fixed atoms. This will
    only be set after the first energy calculation that uses the grid */
GridInfo InterFF::grid() const
{
    return d->fixed_atoms.grid();
}

/** Internal function used to rebox changed atoms */
void InterFF::reboxChangedAtoms()
{
    const CLJDelta *deltas = wspace.constData();

    for (int i=0; i<wspace.nDeltas(); ++i)
    {
        atom_locs[ deltas[i].ID() ] = cljboxes.apply( deltas[i] );
    }

    wspace.clear();
    needs_reboxing = false;
}

/** Internal function used to regrid the atoms */
void InterFF::regridAtoms()
{
    if (usesGrid() and not d.constData()->fixed_atoms.isEmpty())
    {
        d->fixed_atoms.setGridDimensions( cljboxes.atoms() );
    }
}

/** Signal that this forcefield must now be recalculated from scratch */
void InterFF::mustNowRecalculateFromScratch()
{
    if (not wspace.isEmpty())
    {
        reboxChangedAtoms();
    }
    
    this->setDirty();
}

/** Signal to completely do everything from scratch */
void InterFF::mustNowReallyRecalculateFromScratch()
{
    atom_locs.clear();
    d->mol_to_idx.clear();
    d->molres_to_idx.clear();
    d->atom_selections.clear();
    wspace = CLJWorkspace();
    needs_reboxing = false;
    mustNowRecalculateFromScratch();
}

void InterFF::reextractAtoms(const PartialMolecule &pmol)
{
    if (pmol.selectedAll())
    {
        Molecule m = pmol.molecule();
        
        if (m.nResidues() == 1)
        {
            //this is a single-residue molecule, e.g. a solvent or ligand
            CLJAtoms mol( m, d.constData()->maps_for_mol.value(m.number(),PropertyMap()) );
            
            QVector<CLJBoxIndex> idxs = cljboxes.add(mol);
            
            const quint32 id = atom_locs.count();
            
            atom_locs.append(idxs);
            d->atom_selections.append( AtomSelection() );
            
            d->mol_to_idx.insert(m.number(), id);
        }
        else
        {
            //add each residue individually. This is to optimise memory
            //access and the energy calculation when only one residue is moved
            //at a time
            PropertyMap map = d.constData()->maps_for_mol.value(m.number(), PropertyMap());
            
            for (int i=0; i<m.nResidues(); ++i)
            {
                Residue r = m.residue(ResIdx(i));
                
                CLJAtoms res( r, map );
                
                QVector<CLJBoxIndex> idxs = cljboxes.add(res);
                
                const quint32 id = atom_locs.count();
                
                atom_locs.append(idxs);
                d->atom_selections.append( AtomSelection() );
                
                if (i == 0)
                    //save the index to the first residue's atoms
                    d->mol_to_idx.insert(m.number(), id);
                
                d->molres_to_idx.insert( MolResNum(m.number(),r.number()), id );
            }
        }
    }
    else
    {
        if (pmol.molecule().nResidues() == 1)
        {
            //we have selected part of a single-residue molecule
            //this is a single-residue molecule, e.g. a solvent or ligand
            CLJAtoms mol( pmol, d.constData()->maps_for_mol.value(pmol.number(),PropertyMap()) );
            
            QVector<CLJBoxIndex> idxs = cljboxes.add(mol);
            
            const quint32 id = atom_locs.count();
            
            atom_locs.append(idxs);
            d->atom_selections.append(pmol.selection());
            
            d->mol_to_idx.insert(pmol.number(), id);
        }
        else
        {
            //loop through and add all of the selected residues
            AtomSelection selected_atoms = pmol.selection();
            PropertyMap map = d.constData()->maps_for_mol.value(pmol.number(), PropertyMap());

            bool have_first = false;

            foreach (ResIdx residx, selected_atoms.selectedResidues())
            {
                if (selected_atoms.selectedAll(residx))
                {
                    Residue r = pmol.residue(residx);
                
                    CLJAtoms res( r, map );
                    
                    QVector<CLJBoxIndex> idxs = cljboxes.add(res);
                    
                    const quint32 id = atom_locs.count();
                    
                    atom_locs.append(idxs);
                    d->atom_selections.append( AtomSelection() );
                    
                    d->molres_to_idx.insert( MolResNum(pmol.number(),r.number()), id );
                    
                    if (not have_first)
                    {
                        d->mol_to_idx.insert(pmol.number(), id);
                        have_first = false;
                    }
                }
                else
                {
                    AtomSelection res_selection(selected_atoms);
                    res_selection = res_selection.mask(residx);
                    
                    CLJAtoms res( PartialMolecule(pmol.data(),res_selection), map );
                    
                    QVector<CLJBoxIndex> idxs = cljboxes.add(res);
                    
                    const quint32 id = atom_locs.count();
                    
                    atom_locs.append(idxs);
                    d->atom_selections.append( res_selection );
                    
                    d->molres_to_idx.insert( MolResNum(pmol.number(),
                                                       pmol.data().info().number(residx)), id );

                    if (not have_first)
                    {
                        d->mol_to_idx.insert(pmol.number(), id);
                        have_first = false;
                    }
                }
            }
        }
    }
}

/** Re-extract all of the atoms - this also reboxes them all */
void InterFF::reextractAtoms()
{
    QElapsedTimer t;
    t.start();

    mustNowReallyRecalculateFromScratch();
    
    //we know that we are only a single molecule group
    const Molecules &mols = this->group( MGIdx(0) ).molecules();
    
    if (mols.isEmpty())
        return;
    
    atom_locs.clear();
    atom_locs.reserve(mols.nMolecules());
    d->mol_to_idx.clear();
    d->mol_to_idx.reserve(mols.nMolecules());
    d->molres_to_idx.clear();
    d->atom_selections.clear();
    d->atom_selections.reserve(mols.nMolecules());
    cljboxes = CLJBoxes();
    
    for (Molecules::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        reextractAtoms(*it);
    }
    
    qint64 ns = t.nsecsElapsed();
    
    qDebug() << "Re-extracting took" << (0.000001*ns) << "ms";
}

/** Recalculate the energy of this forcefield */
void InterFF::recalculateEnergy()
{
    if (needs_reboxing)
        this->reboxChangedAtoms();

    if (atom_locs.isEmpty())
    {
        //everything needs to be recalculated from scratch
        //extract all of the atoms
        reextractAtoms();
            
        //if there are no atoms, then nothing to be done
        if (atom_locs.isEmpty())
        {
            d.constData()->cljcomps.setEnergy(*this, CLJEnergy(0.0));
            this->setClean();
            return;
        }
        
        //calculate the energy from scratch
        tuple<double,double> nrgs(0,0);
        
        if (not d.constData()->fixed_only)
        {
            QElapsedTimer t;
            t.start();
            CLJCalculator calc;
            nrgs = calc.calculate(*(d.constData()->cljfunc), cljboxes);
            //tuple<double,double> nrgs = d.constData()->cljfunc->calculate(cljboxes);
            qint64 ns = t.nsecsElapsed();
            
            qDebug() << "Calculating total energy took" << (0.000001*ns) << "ms";
        }

        if (not d.constData()->fixed_atoms.isEmpty())
        {
            this->regridAtoms();
            tuple<double,double> grid_nrgs = d.constData()->fixed_atoms.calculate(cljboxes);
            nrgs.get<0>() += grid_nrgs.get<0>();
            nrgs.get<1>() += grid_nrgs.get<1>();
        }
        
        d.constData()->cljcomps.setEnergy(*this, CLJEnergy(nrgs.get<0>(), nrgs.get<1>()));
        this->setClean();
    }
    else if (not wspace.isEmpty())
    {
        tuple<double,double> delta_nrgs(0,0);
        
        CLJDelta delta = wspace.merge();
        
        if (not d.constData()->fixed_only)
        {
            //calculate the change in energy using the molecules in changed_atoms
            CLJCalculator calc;
            delta_nrgs = calc.calculate(*(d.constData()->cljfunc),
                                        delta, cljboxes);
        }
        
        if (not d.constData()->fixed_atoms.isEmpty())
        {
            tuple<double,double> grid_deltas = d.constData()
                                                    ->fixed_atoms.calculate(delta.changedAtoms());
            
            delta_nrgs.get<0>() += grid_deltas.get<0>();
            delta_nrgs.get<1>() += grid_deltas.get<1>();
        }
        
        //changed molecules will need reboxing. We will defer this, in case
        //the move that caused this change is rejected
        needs_reboxing = true;
        
        d.constData()->cljcomps.changeEnergy(*this,
                                    CLJEnergy(delta_nrgs.get<0>(), delta_nrgs.get<1>()));
        this->setClean();
    }
    else
    {
        //recalculate everything from scratch as this has been requested
        //calculate the energy from scratch
        tuple<double,double> nrgs(0,0);
        
        if (not d.constData()->fixed_only)
        {
            CLJCalculator calc;
            nrgs = calc.calculate(*(d.constData()->cljfunc), cljboxes);
            //tuple<double,double> nrgs = d.constData()->cljfunc->calculate(cljboxes);
        }

        if (not d.constData()->fixed_atoms.isEmpty())
        {
            this->regridAtoms();
            tuple<double,double> grid_nrgs = d.constData()->fixed_atoms.calculate(cljboxes);
            nrgs.get<0>() += grid_nrgs.get<0>();
            nrgs.get<1>() += grid_nrgs.get<1>();
        }
        
        d.constData()->cljcomps.setEnergy(*this, CLJEnergy(nrgs.get<0>(), nrgs.get<1>()));
        this->setClean();
    }
}

/** Function called to add a molecule to this forcefield */
void InterFF::_pvt_added(const SireMol::PartialMolecule &mol, const SireBase::PropertyMap &map)
{
    if (needs_reboxing)
        reboxChangedAtoms();

    if (atom_locs.isEmpty())
    {
        //we can wait until everything is reboxed
        mustNowReallyRecalculateFromScratch();
        setDirty();
        return;
    }

    //do we already contain part of this molecule?
    if (d.constData()->mol_to_idx.contains(mol.number()))
    {
        //we cannot easily cope with adding new parts of a molecule that already
        //exists in this forcefield. Let give up and do everything from scratch
        if (d.constData()->maps_for_mol.value(mol.number(), PropertyMap()) != map)
            throw SireError::incompatible_error( QObject::tr(
                    "Trying to add a molecule with a new PropertyMap (%1) when it already "
                    "exists in the forcefield %2 with a different PropertyMap (%3)")
                        .arg(map.toString())
                        .arg(this->toString())
                        .arg(d.constData()->maps_for_mol.value(mol.number(),PropertyMap())
                                .toString()), CODELOC );
        
        mustNowReallyRecalculateFromScratch();
        setDirty();
        return;
    }
    else if (mol.molecule().nResidues() > 1)
    {
        //adding in multiple-residue molecules is too difficult to debug at the moment
        //do it the slow way
        mustNowReallyRecalculateFromScratch();
        setDirty();
        return;
    }
    else
    {
        //we are adding a new single-residue molecule. Just add this and box it with the rest
        if (not map.isDefault())
            d->maps_for_mol.insert(mol.number(), map);
        
        if (mol.selectedAll())
            d->atom_selections.append( AtomSelection() );
        else
            d->atom_selections.append( mol.selection() );
        
        atom_locs.append( QVector<CLJBoxIndex>() );
        const quint32 id = atom_locs.count() - 1;
        d->mol_to_idx.insert(mol.number(), id);
        
        const bool was_single_id = wspace.isSingleID();

        wspace.push( id, cljboxes, QVector<CLJBoxIndex>(),
                     this->group(MGIdx(0)).molecules()[mol.number()], map );

        if (was_single_id)
        {
            if (not wspace.isSingleID())
            {
                //we need to remove all old atoms in the delta from cljboxes, or the
                //delta energy calculation won't work
                for (int i=0; i<wspace.nDeltas(); ++i)
                {
                    cljboxes.remove( wspace.constData()[i].oldIndicies() );
                    atom_locs[ wspace.constData()[i].ID() ] = QVector<CLJBoxIndex>();
                }
            }
        }
        
        this->setDirty();
    }
}

/** Function called to remove a molecule from this forcefield */
void InterFF::_pvt_removed(const SireMol::PartialMolecule &mol)
{
    if (needs_reboxing)
        reboxChangedAtoms();

    //be lazy for the moment - recalculate everything!
    mustNowReallyRecalculateFromScratch();
    
    if (d.constData()->maps_for_mol.contains(mol.number()))
    {
        if (not this->contains(mol.number()))
        {
            d->maps_for_mol.remove(mol.number());
        }
    }
    
    this->setDirty();
}

/** Function called to indicate that a molecule in this forcefield has changed */
void InterFF::_pvt_changed(const SireMol::Molecule &molecule, bool auto_update)
{
    if (auto_update or needs_reboxing)
        reboxChangedAtoms();

    //get the ID number of this molecule in the forcefield
    MolNum molnum = molecule.number();

    if (not d.constData()->mol_to_idx.contains(molnum))
        throw SireError::program_bug( QObject::tr(
                "It should not be possible that the molecule %1 is not in the forcefield already?")
                    .arg(molnum.value()), CODELOC );
    
    //is this a multi-residue molecule?
    if (molecule.nResidues() > 1)
    {
        //for now, do this the painful way!
        mustNowReallyRecalculateFromScratch();
        setDirty();
        return;
    }
    else
    {
        const bool was_single_id = wspace.isSingleID();
    
        //does this forcefield contain the whole molecule
        const quint32 id = d.constData()->mol_to_idx.value(molnum);
        
        //have we changed this molecule before? If so, then it is safest that
        //we redo everything from scratch
        for (int i=0; i<wspace.nDeltas(); ++i)
        {
            if (wspace.constData()[i].ID() == id)
            {
                //yes, we have done this before
                mustNowReallyRecalculateFromScratch();
                setDirty();
                return;
            }
        }
        
        if (d.constData()->atom_selections[id].isEmpty())
        {
            //we have the whole molecule
            wspace.push( id, cljboxes, atom_locs.at(id), molecule,
                         d.constData()->maps_for_mol.value(molnum,PropertyMap()) );
        }
        else
        {
            //we have only a part of the molecule
            wspace.push( id, cljboxes, atom_locs.at(id),
                         PartialMolecule(molecule.data(), d.constData()->atom_selections[id]),
                         d.constData()->maps_for_mol.value(molnum,PropertyMap()) );
        }

        if (was_single_id)
        {
            if (not wspace.isSingleID())
            {
                //we need to remove all old atoms in the delta from cljboxes, or the
                //delta energy calculation won't work
                for (int i=0; i<wspace.nDeltas(); ++i)
                {
                    cljboxes.remove( wspace.constData()[i].oldIndicies() );
                    atom_locs[ wspace.constData()[i].ID() ] = QVector<CLJBoxIndex>();
                }
            }
        }
        else
        {
            //we need to remove the old atoms from the cljboxes, or the delta energy
            //calculation won't work
            cljboxes.remove( atom_locs.at(id) );
            atom_locs[id] = QVector<CLJBoxIndex>();
        }
    
        this->setDirty();
    }
    
    if (auto_update and this->needsAccepting())
        this->accept();
}

/** Function called to indicate that a list of molecules in this forcefield have changed */
void InterFF::_pvt_changed(const QList<SireMol::Molecule> &molecules, bool auto_update)
{
    if (not atom_locs.isEmpty())
    {
        foreach (const Molecule &molecule, molecules)
        {
            this->_pvt_changed(molecule, auto_update);
            
            if (atom_locs.isEmpty())
                break;
        }
    }
    
    setDirty();
}

/** Function called to indicate that all molecules in this forcefield have been removed */
void InterFF::_pvt_removedAll()
{
    if (needs_reboxing)
        reboxChangedAtoms();

    atom_locs.clear();
    cljboxes = CLJBoxes();
    wspace = CLJWorkspace();
    d->maps_for_mol.clear();
    d->mol_to_idx.clear();

    mustNowReallyRecalculateFromScratch();
    
    this->setDirty();
}

/** Function called to query whether or not a change in source properties would
    change the properties needed by this forcefield for the molecule with number 'molnum' */
bool InterFF::_pvt_wouldChangeProperties(SireMol::MolNum molnum,
                                         const SireBase::PropertyMap &map) const
{
    if (map.isDefault())
    {
        return d->maps_for_mol.contains(molnum);
    }
    else
    {
        return not d->maps_for_mol.contains(molnum);
    }
}

/** Return whether or not this forcefield is using a temporary workspace that 
    needs to be accepted */
bool InterFF::needsAccepting() const
{
    return needs_reboxing or G1FF::needsAccepting();
}

/** Tell the forcefield that the last move was accepted. This tells the
    forcefield to make permanent any temporary changes that were used a workspace
    to avoid memory allocation during a move */
void InterFF::accept()
{
    if (needs_reboxing)
        this->reboxChangedAtoms();
    
    G1FF::accept();
}
