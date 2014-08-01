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

#include "cljextractor.h"

#include "SireMol/atomcoords.h"
#include "SireMol/atomcharges.h"
#include "SireMol/residue.h"
#include "SireMM/atomljs.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<CLJExtractor> r_cljext( NO_ROOT );

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJExtractor &cljext)
{
    writeHeader(ds, r_cljext, 1);
    
    SharedDataStream sds(ds);
    
    sds << cljext.mol << cljext.newmol
        << cljext.selected_atoms << cljext.new_selected_atoms
        << cljext.props << cljext.cljidxs
        << cljext.cljdeltas
        << cljext.id_source
        << cljext.extract_by_residue;

    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJExtractor &cljext)
{
    VersionID v = readHeader(ds, r_cljext);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> cljext.mol >> cljext.newmol
            >> cljext.selected_atoms >> cljext.new_selected_atoms
            >> cljext.props >> cljext.cljidxs
            >> cljext.cljdeltas
            >> cljext.id_source
            >> cljext.extract_by_residue;
    }
    else
        throw version_error(v, "1", r_cljext, CODELOC);

    return ds;
}

/** Null constructor */
CLJExtractor::CLJExtractor() : id_source(CLJAtoms::USE_MOLNUM), extract_by_residue(true)
{}

/** Construct to extract the CLJ properties from the passed molecule, extracting
    information per-residue, and using the supplied property map to find the 
    correct properties */
CLJExtractor::CLJExtractor(const MoleculeView &molecule, const PropertyMap &map)
             : props(map), id_source(CLJAtoms::USE_MOLNUM), extract_by_residue(true)
{
    mol = molecule.molecule();
    newmol = mol;

    if (not mol.selectedAll())
    {
        selected_atoms = mol.selection();
        new_selected_atoms = selected_atoms;
    }
}

/** Construct to extract the CLJ properties from the passed molecule, specifying
    wheter or not to extract information per-residue or per-molecule, 
    and using the supplied property map to find the
    correct properties */
CLJExtractor::CLJExtractor(const MoleculeView &molecule, bool split_by_residue,
                           const PropertyMap &map)
             : props(map), id_source(CLJAtoms::USE_MOLNUM), extract_by_residue(split_by_residue)
{
    mol = molecule.molecule();
    newmol = mol;

    if (not mol.selectedAll())
    {
        selected_atoms = mol.selection();
        new_selected_atoms = selected_atoms;
    }
}

/** Construct to extract the CLJ properties from the passed molecule, extracting
    information per-residue, and using the supplied property map to find the 
    correct properties */
CLJExtractor::CLJExtractor(const MoleculeView &molecule, CLJAtoms::ID_SOURCE id,
                           const PropertyMap &map)
             : props(map), id_source(id), extract_by_residue(true)
{
    mol = molecule.molecule();
    
    if (not mol.selectedAll())
    {
        selected_atoms = mol.selection();
    }
}

/** Construct to extract the CLJ properties from the passed molecule, specifying
    wheter or not to extract information per-residue or per-molecule, 
    and using the supplied property map to find the
    correct properties */
CLJExtractor::CLJExtractor(const MoleculeView &molecule, CLJAtoms::ID_SOURCE id,
                           bool split_by_residue, const PropertyMap &map)
             : props(map), id_source(id), extract_by_residue(split_by_residue)
{
    mol = molecule.molecule();
    newmol = mol;
    
    if (not mol.selectedAll())
    {
        selected_atoms = mol.selection();
        new_selected_atoms = selected_atoms;
    }
}

/** Copy constructor */
CLJExtractor::CLJExtractor(const CLJExtractor &other)
             : mol(other.mol), selected_atoms(other.selected_atoms),
               newmol(other.newmol), new_selected_atoms(other.new_selected_atoms),
               props(other.props),
               cljidxs(other.cljidxs), cljdeltas(other.cljdeltas),
               id_source(other.id_source), extract_by_residue(other.extract_by_residue)
{}

/** Destructor */
CLJExtractor::~CLJExtractor()
{}

/** Copy assignment operator */
CLJExtractor& CLJExtractor::operator=(const CLJExtractor &other)
{
    if (this != &other)
    {
        mol = other.mol;
        selected_atoms = other.selected_atoms;
        newmol = other.newmol;
        new_selected_atoms = other.new_selected_atoms;
        props = other.props;
        cljidxs = other.cljidxs;
        cljdeltas = other.cljdeltas;
        id_source = other.id_source;
        extract_by_residue = other.extract_by_residue;
    }
    
    return *this;
}

/** Comparison operator */
bool CLJExtractor::operator==(const CLJExtractor &other) const
{
    return (this == &other) or
           (mol.number() == other.mol.number() and mol.version() == other.mol.version() and
            newmol.number() == other.newmol.number() and
            newmol.version() == other.newmol.version() and
            selected_atoms == other.selected_atoms and
            new_selected_atoms == other.new_selected_atoms and
            props == other.props and
            cljidxs == other.cljidxs and cljdeltas == other.cljdeltas and
            id_source == other.id_source and extract_by_residue == other.extract_by_residue);
}

/** Comparison operator */
bool CLJExtractor::operator!=(const CLJExtractor &other) const
{
    return not operator==(other);
}

/** Return whether or not this molecule has been changed during the move */
bool CLJExtractor::changed() const
{
    return newmol.version() != mol.version() or
           new_selected_atoms != selected_atoms;
}

/** Return whether or not this extractor needs to be committed (i.e. whether
    or not the molecule has changed in any way) */
bool CLJExtractor::needsCommitting() const
{
    return changed();
}

/** Return whether or not this extractor is empty (contains no atoms) */
bool CLJExtractor::isEmpty() const
{
    return selected_atoms.selectedNone() and new_selected_atoms.selectedNone();
}

/** Return whether or not this extractor is null (contains no molecule information) */
bool CLJExtractor::isNull() const
{
    return mol.isNull();
}

/** Return the molecule as it exists before any changes were made */
PartialMolecule CLJExtractor::oldmolecule() const
{
    if (selected_atoms.isNull())
        return mol;
    else
        return PartialMolecule(mol, selected_atoms);
}

/** Return the molecule as it exists after the changes have been made */
PartialMolecule CLJExtractor::newMolecule() const
{
    if (new_selected_atoms.isNull())
        return newmol;
    else
        return PartialMolecule(newmol, new_selected_atoms);
}

/** Return the property map used to get the names of the coordinates,
    charge and LJ properties from the molecule */
PropertyMap CLJExtractor::propertyMap() const
{
    return props;
}

/** Return the property used to find the coordinates */
PropertyName CLJExtractor::coordinatesProperty() const
{
    return props["coordinates"];
}

/** Return the property used to find the charges */
PropertyName CLJExtractor::chargeProperty() const
{
    return props["charge"];
}

/** Return the property used to find the LJ parameters */
PropertyName CLJExtractor::ljProperty() const
{
    return props["LJ"];
}

/** Return whether or not atoms are extracted by residue */
bool CLJExtractor::extractingByResidue() const
{
    return extract_by_residue;
}

/** Return whether or not atoms are extracted by molecule */
bool CLJExtractor::extractingByMolecule() const
{
    return not extract_by_residue;
}

/** Return the source of the ID property for each CLJAtom */
CLJAtoms::ID_SOURCE CLJExtractor::idSource() const
{
    return id_source;
}

/** Add the extra atoms in 'new_molecule' to the molecule */
void CLJExtractor::add(const MoleculeView &new_molecule, CLJBoxes &boxes, CLJWorkspace &workspace)
{
    this->update(new_molecule, boxes, workspace);
    this->add(new_molecule.selection(), boxes, workspace);
}

/** Add the extra atoms in 'new_selection' to the molecule */
void CLJExtractor::add(const AtomSelection &new_selection,
                       CLJBoxes &boxes, CLJWorkspace &workspace)
{
    if (new_selected_atoms.isNull() or new_selected_atoms.selectedAll() or
           new_selected_atoms.contains(new_selection))
    {
        //we are not adding any more atoms
        return;
    }
    else
    {
        AtomSelection news = new_selected_atoms;
        news = news.unite(new_selection);
    
        this->updateSelection(news, boxes, workspace);
    }
}

/** Remove the atoms in 'new_molecule' from the molecule */
void CLJExtractor::remove(const MoleculeView &new_molecule,
                          CLJBoxes &boxes, CLJWorkspace &workspace)
{
    this->update(new_molecule, boxes, workspace);
    this->remove(new_molecule.selection(), boxes, workspace);
}

/** Remove the atoms in 'new_selection' from the molecule */
void CLJExtractor::remove(const AtomSelection &new_selection,
                          CLJBoxes &boxes, CLJWorkspace &workspace)
{
    if (new_selected_atoms.isNull() or not new_selected_atoms.selectedNone())
    {
        //we have atoms to remove
        AtomSelection news = new_selected_atoms;
        news = news.subtract(new_selection);
        
        this->updateSelection(news, boxes, workspace);
    }
    else
    {
        //we have selected no atoms, so cannot remove any more
        return;
    }
}

/** This function is used to update the molecule to use only the passed
    selection as the selected atoms */
void CLJExtractor::updateSelection(const AtomSelection &new_selection,
                                   CLJBoxes &boxes, CLJWorkspace &workspace)
{
    if (cljidxs.isEmpty())
        this->initialise(boxes, workspace);

    //have to use 'old_selection' as we use a null AtomSelection to indicate
    //that the entire molecule has been selected
    AtomSelection old_selection = new_selected_atoms;

    if (new_selected_atoms.isNull())
        old_selection = mol.selection();

    old_selection.assertCompatibleWith(new_selection);

    if (old_selection == new_selection)
        //nothing to do
        return;
    
    if (new_selection.selectedNone())
    {
        //all atoms have been removed
        for (int i=0; i<cljidxs.count(); ++i)
        {
            cljdeltas[i] = workspace.push(boxes, cljidxs.at(i), CLJAtoms(), cljdeltas[i]);
        }
    }
    else if (new_selection.selectedAll())
    {
        if (extract_by_residue)
        {
            for (ResIdx i(0); i<mol.nResidues(); ++i)
            {
                if (not old_selection.selectedAll(i))
                {
                    cljdeltas[i] = workspace.push(boxes, cljidxs.at(i),
                                                  CLJAtoms(newmol.residue(i), id_source, props),
                                                  cljdeltas[i]);
                }
            }
        }
        else
        {
            cljdeltas[0] = workspace.push(boxes, cljidxs.at(0),
                                          CLJAtoms(newmol, id_source, props),
                                          cljdeltas[0]);
        }
    }
    else
    {
        //something in between...
        if (extract_by_residue)
        {
            for (ResIdx i(0); i<mol.nResidues(); ++i)
            {
                AtomSelection old_res_selection = old_selection;
                AtomSelection new_res_selection = new_selection;
                
                old_res_selection = old_res_selection.mask(i);
                new_res_selection = new_res_selection.mask(i);
            
                if (old_res_selection != new_res_selection)
                {
                    cljdeltas[i] = workspace.push(boxes, cljidxs.at(i),
                                            CLJAtoms( PartialMolecule(newmol, new_res_selection),
                                                      id_source, props ), cljdeltas[i]);
                }
            }
        }
        else
        {
            cljdeltas[0] = workspace.push(boxes, cljidxs.at(0),
                                          CLJAtoms( PartialMolecule(newmol,new_selection),
                                                    id_source, props ), cljdeltas[0]);
        }
    }
    
    if (new_selection.selectedAll())
    {
        new_selected_atoms = AtomSelection();
    }
    else
    {
        new_selected_atoms = new_selection;
    }
}

/** Update the molecule, calculating the change in CLJAtoms as a CLJDelta that is
    added to the passed CLJWorkspace. Any atoms that have changed are removed
    from the passed CLJBoxes */
void CLJExtractor::update(const MoleculeView &new_molecule,
                          CLJBoxes &boxes, CLJWorkspace &workspace)
{
    if (new_molecule.data().number() != newmol.number())
        throw SireError::incompatible_error( QObject::tr(
                "You cannot update one molecule (%1) with data from another molecule!")
                    .arg(mol.toString()).arg(new_molecule.molecule().toString()),
                        CODELOC );
    
    if (new_molecule.data().version() == newmol.version())
        //there is nothing to update
        return;

    else if (cljidxs.isEmpty())
    {
        //we haven't been added to the boxes or workspace yet. We need to completely
        //re-add the molecule
        this->initialise(boxes, workspace);
        return;
    }
    else
    {
        //we are updating the molecule - see if we need to update the
        //coordinates, charge or LJ properties...
        const PropertyName coords_property = coordinatesProperty();
        const PropertyName charge_property = chargeProperty();
        const PropertyName lj_property = ljProperty();

        bool changed_coords = newmol.version(coords_property) !=
                                    new_molecule.data().version(coords_property);

        bool changed_charge = newmol.version(charge_property) !=
                                    new_molecule.data().version(charge_property);
        
        bool changed_lj = newmol.version(lj_property) !=
                                    new_molecule.data().version(lj_property);
        
        if (not (changed_coords or changed_charge or changed_lj))
        {
            //nothing important has changed :-)
            newmol = new_molecule.molecule();
            return;
        }
        
        //do we have multiple CLJAtoms groups to extract?
        if (extract_by_residue and newmol.nResidues() > 1)
        {
            //generate a list of changed residue indicies
            QSet<qint32> changed_residues;
            
            if (changed_coords)
            {
                const AtomCoords &old_coords = newmol.property(coords_property).asA<AtomCoords>();
                const AtomCoords &new_coords = new_molecule.data()
                                                    .property(coords_property).asA<AtomCoords>();
        
                for (CGIdx i(0); i<newmol.nCutGroups(); ++i)
                {
                    const Vector *oldc = old_coords.constData(i);
                    const Vector *newc = new_coords.constData(i);
                    
                    if (oldc != newc)
                    {
                        for (Index j(0); j<mol.data().info().nAtoms(i); ++j)
                        {
                            if (oldc[j] != newc[j])
                            {
                                //this atom has moved
                                changed_residues.insert(
                                    newmol.data().info().parentResidue( CGAtomIdx(i,j) ).value() );
                            }
                        }
                    }
                }
            }
            
            if (changed_charge and changed_residues.count() < newmol.nResidues())
            {
                const AtomCharges &old_chgs = newmol.property(charge_property).asA<AtomCharges>();
                const AtomCharges &new_chgs = new_molecule.data()
                                                    .property(charge_property).asA<AtomCharges>();
            
                for (CGIdx i(0); i<newmol.nCutGroups(); ++i)
                {
                    const Charge *oldc = old_chgs.constData(i);
                    const Charge *newc = new_chgs.constData(i);
                    
                    if (oldc != newc)
                    {
                        for (Index j(0); j<mol.data().info().nAtoms(i); ++j)
                        {
                            if (oldc[j] != newc[j])
                            {
                                //the charge on this atom has changed
                                changed_residues.insert(
                                    newmol.data().info().parentResidue( CGAtomIdx(i,j) ).value() );
                            }
                        }
                    }
                }
            }
            
            if (changed_lj and changed_residues.count() < newmol.nResidues())
            {
                const AtomLJs &old_ljs = newmol.property(lj_property).asA<AtomLJs>();
                const AtomLJs &new_ljs = new_molecule.data().property(lj_property).asA<AtomLJs>();

                for (CGIdx i(0); i<newmol.nCutGroups(); ++i)
                {
                    const LJParameter *oldlj = old_ljs.constData(i);
                    const LJParameter *newlj = new_ljs.constData(i);
                    
                    if (oldlj != newlj)
                    {
                        for (Index j(0); j<mol.data().info().nAtoms(i); ++j)
                        {
                            if (oldlj[j] != newlj[j])
                            {
                                //the LJ parameter of this atom has changed
                                changed_residues.insert(
                                    newmol.data().info().parentResidue( CGAtomIdx(i,j) ).value() );
                            }
                        }
                    }
                }
            }
        
            newmol = new_molecule.molecule();
        
            //loop over all of the changed residues
            foreach (int changed_residue, changed_residues)
            {
                ResIdx i(changed_residue);
            
                if (new_selected_atoms.selectedAll() or new_selected_atoms.selectedAll(i))
                {
                    cljdeltas[i] = workspace.push(boxes, cljidxs.at(i),
                                                  CLJAtoms(newmol.residue(i), id_source, props),
                                                  cljdeltas[i]);
                }
                else
                {
                    AtomSelection selected_resatoms = new_selected_atoms;
                    selected_resatoms = selected_resatoms.intersect(i);
                
                    cljdeltas[i] = workspace.push(boxes, cljidxs.at(i),
                                     CLJAtoms( PartialMolecule(newmol, selected_resatoms),
                                               id_source, props ), cljdeltas[i]);
                }
            }
        }
        else
        {
            newmol = new_molecule.molecule();

            //there is only one residue, or we are extracting as a single molecule
            if (new_selected_atoms.isNull())
            {
                //we have selected all atoms
                cljdeltas[0] = workspace.push(boxes, cljidxs.at(0),
                                              CLJAtoms(newmol, id_source, props),
                                              cljdeltas[0]);
            }
            else
            {
                cljdeltas[0] = workspace.push(boxes, cljidxs.at(0),
                                    CLJAtoms( PartialMolecule(newmol.data(),new_selected_atoms),
                                              id_source, props ), cljdeltas[0] );
            }
        }
    }
}

/** Commit the changes */
void CLJExtractor::commit(CLJBoxes &boxes, CLJWorkspace &workspace)
{
    mol = newmol;
    selected_atoms = new_selected_atoms;
    
    for (int i=0; i<cljidxs.count(); ++i)
    {
        cljidxs[i] = workspace.commit(boxes, cljdeltas[i]);
        cljdeltas[i] = CLJDelta();
    }
}

/** Revert the changes */
void CLJExtractor::revert(CLJBoxes &boxes, CLJWorkspace &workspace)
{
    newmol = mol;
    new_selected_atoms = selected_atoms;
    
    //have we only initialised this extractor?
    bool initialised_only = true;
    
    for (int i=0; i<cljidxs.count(); ++i)
    {
        if (not cljidxs.at(i).isEmpty())
        {
            initialised_only = false;
            break;
        }
    }
    
    if (initialised_only)
    {
        cljidxs.clear();
        cljdeltas.clear();
    }
    else
    {
        for (int i=0; i<cljidxs.count(); ++i)
        {
            cljidxs[i] = workspace.revert(boxes, cljdeltas[i]);
            cljdeltas[i] = CLJDelta();
        }
    }
}

/** Initialise this molecule with the passed CLJBoxes / CLJWorkspace */
void CLJExtractor::initialise(CLJBoxes &boxes, CLJWorkspace &workspace)
{
    if (not cljidxs.isEmpty())
        throw SireError::program_bug( QObject::tr(
                "It is a mistake to initialise a CLJExtractor more than once!"),
                    CODELOC );
    
    if (extract_by_residue)
    {
        cljdeltas = QVector<CLJDelta>( newmol.nResidues(), CLJDelta() );
        cljidxs = QVector< QVector<CLJBoxIndex> >( newmol.nResidues(), QVector<CLJBoxIndex>() );

        for (ResIdx i(0); i<newmol.nResidues(); ++i)
        {
            if (new_selected_atoms.isNull() or new_selected_atoms.selectedAll(i))
            {
                cljdeltas[i] = workspace.push(boxes, QVector<CLJBoxIndex>(),
                                              CLJAtoms(newmol.residue(i), id_source, props),
                                              CLJDelta());
            }
            else
            {
                AtomSelection res_selection = new_selected_atoms;
                res_selection = res_selection.mask(i);
                
                cljdeltas[i] = workspace.push(boxes, QVector<CLJBoxIndex>(),
                                              CLJAtoms( PartialMolecule(newmol,res_selection),
                                                        id_source, props ), CLJDelta());
            }
        }
    }
    else
    {
        cljdeltas = QVector<CLJDelta>( 1, CLJDelta() );
        cljidxs = QVector< QVector<CLJBoxIndex> >( 1, QVector<CLJBoxIndex>() );
        
        if (new_selected_atoms.isNull())
        {
            cljdeltas[0] = workspace.push(boxes, QVector<CLJBoxIndex>(),
                                          CLJAtoms(newmol, id_source, props), CLJDelta());
        }
        else
        {
            cljdeltas[0] = workspace.push(boxes, QVector<CLJBoxIndex>(),
                                          CLJAtoms( PartialMolecule(newmol,new_selected_atoms),
                                                    id_source, props ), CLJDelta());
        }
    }
}
