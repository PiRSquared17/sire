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

#include "moleditor.h"

#include "segeditor.h"
#include "chaineditor.h"
#include "reseditor.h"
#include "cgeditor.h"
#include "atomeditor.h"

#include "segment.h"
#include "chain.h"
#include "residue.h"
#include "cutgroup.h"
#include "atom.h"

#include "selector.hpp"
#include "mover.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireStream;

//instantiate the underlying Editor class
template class Editor<MolEditor, Molecule>;

//////////
////////// Implementation of MolEditor
//////////

static const RegisterMetaType<MolEditor> r_moleditor;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const MolEditor &moleditor)
{
    writeHeader(ds, r_moleditor, 1);
    
    ds << static_cast<const Editor<MolEditor, Molecule>&>(moleditor);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       MolEditor &moleditor)
{
    VersionID v = readHeader(ds, r_moleditor);
    
    if (v == 1)
    {
        ds >> static_cast<Editor<MolEditor, Molecule>&>(moleditor);
    }
    else
        throw version_error(v, "1", r_moleditor, CODELOC);
        
    return ds;
}

/** Null constructor */
MolEditor::MolEditor() : Editor<MolEditor, Molecule>()
{}

/** Construct an editor to edit a copy of 'molecule' */
MolEditor::MolEditor(const Molecule &molecule)
          : Editor<MolEditor, Molecule>(molecule)
{}

/** Copy constructor */
MolEditor::MolEditor(const MolEditor &other)
          : Editor<MolEditor, Molecule>(other)
{}

/** Destructor */
MolEditor::~MolEditor()
{}

/** Assign so that this edits a copy of 'molecule' */
MolEditor& MolEditor::operator=(const Molecule &molecule)
{
    Editor<MolEditor, Molecule>::operator=(molecule);
    return *this;
}

/** Copy assignment operator */
MolEditor& MolEditor::operator=(const MolEditor &other)
{
    Editor<MolEditor, Molecule>::operator=(other);
    return *this;
}

/** Return an editor for the atom at ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomEditor MolEditor::select(const AtomID &atomid) const
{
    return AtomEditor( Molecule::select(atomid) );
}

/** Return an editor for the CutGroup at ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CGEditor MolEditor::select(const CGID &cgid) const
{
    return CGEditor( Molecule::select(cgid) );
}

/** Return an editor for the residue at ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResEditor MolEditor::select(const ResID &resid) const
{
    return ResEditor( Molecule::select(resid) );
}

/** Return an editor for the chain at ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
ChainEditor MolEditor::select(const ChainID &chainid) const
{
    return ChainEditor( Molecule::select(chainid) );
}

/** Return an editor for the segment at ID 'segid'

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
SegEditor MolEditor::select(const SegID &segid) const
{
    return SegEditor( Molecule::select(segid) );
}

/** Return an editor for the atom at ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomEditor MolEditor::atom(const AtomID &atomid) const
{
    return AtomEditor( Molecule::atom(atomid) );
}

/** Return an editor for the CutGroup at ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CGEditor MolEditor::cutGroup(const CGID &cgid) const
{
    return CGEditor( Molecule::cutGroup(cgid) );
}

/** Return an editor for the residue at ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResEditor MolEditor::residue(const ResID &resid) const
{
    return ResEditor( Molecule::residue(resid) );
}

/** Return an editor for the chain at ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
ChainEditor MolEditor::chain(const ChainID &chainid) const
{
    return ChainEditor( Molecule::chain(chainid) );
}

/** Return an editor for the segment at ID 'segid'

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
SegEditor MolEditor::segment(const SegID &segid) const
{
    return SegEditor( Molecule::segment(segid) );
}

/** Rename this molecule to 'newname' */
MolEditor& MolEditor::rename(const MolName &newname)
{
    if (newname == this->name())
        //nothing needs doing
        return *this;
        
    d->rename(newname);
    
    return *this;
}

/** Give this molecule a new, unique ID number */
MolEditor& MolEditor::renumber()
{
    d->renumber();
    return *this;
}

/** Give this molecule the ID number 'newnum' */
MolEditor& MolEditor::renumber(MolNum newnum)
{
    if (newnum == this->number())
        //nothing needs doing
        return *this;
        
    d->renumber(newnum);
    
    return *this;
}

/** Add an atom called 'name' and return an editor that can
    be used to edit it */
AtomStructureEditor MolEditor::add(const AtomName &name) const
{
    MolStructureEditor editor(*this);
    return editor.add(name);
}

/** Add an atom with number 'number' and return an editor that
    can be used to edit it */
AtomStructureEditor MolEditor::add(const AtomNum &number) const
{
    MolStructureEditor editor(*this);
    return editor.add(number);
}

/** Add a residue called 'name' and return an editor that can
    be used to edit it */
ResStructureEditor MolEditor::add(const ResName &name) const
{
    MolStructureEditor editor(*this);
    return editor.add(name);
}

/** Add a residue with number 'number' and return an editor that
    can be used to edit it */
ResStructureEditor MolEditor::add(const ResNum &number) const
{
    MolStructureEditor editor(*this);
    return editor.add(number);
}

/** Add a CutGroup called 'name' and return an editor that can
    be used to edit it */
CGStructureEditor MolEditor::add(const CGName &name) const
{
    MolStructureEditor editor(*this);
    return editor.add(name);
}

/** Add a chain called 'name' and return an editor that can
    be used to edit it */
ChainStructureEditor MolEditor::add(const ChainName &name) const
{
    MolStructureEditor editor(*this);
    return editor.add(name);
}

/** Add a segment called 'name' and return an editor that can
    be used to edit it */
SegStructureEditor MolEditor::add(const SegName &name) const
{
    MolStructureEditor editor(*this);
    return editor.add(name);
}

/** Remove all atoms from this molecule that have the ID 'atomid'.
    This does nothing if there are no atoms with this ID.
    This returns an editor that can be used to further edit
    the molecule */
MolStructureEditor MolEditor::remove(const AtomID &atomid) const
{
    MolStructureEditor editor(*this);
    return editor.remove(atomid);
}

/** Remove all CutGroups from this molecule that have the ID 'atomid'.
    This does nothing if there are no CutGroups with this ID.
    This returns an editor that can be used to further edit
    the molecule */
MolStructureEditor MolEditor::remove(const CGID &cgid) const
{
    MolStructureEditor editor(*this);
    return editor.remove(cgid);
}

/** Remove all residues from this molecule that have the ID 'atomid'.
    This does nothing if there are no residues with this ID.
    This returns an editor that can be used to further edit
    the molecule */
MolStructureEditor MolEditor::remove(const ResID &resid) const
{
    MolStructureEditor editor(*this);
    return editor.remove(resid);
}

/** Remove all chains from this molecule that have the ID 'atomid'.
    This does nothing if there are no chains with this ID.
    This returns an editor that can be used to further edit
    the molecule */
MolStructureEditor MolEditor::remove(const ChainID &chainid) const
{
    MolStructureEditor editor(*this);
    return editor.remove(chainid);
}

/** Remove all segments from this molecule that have the ID 'atomid'.
    This does nothing if there are no segments with this ID.
    This returns an editor that can be used to further edit
    the molecule */
MolStructureEditor MolEditor::remove(const SegID &segid) const
{
    MolStructureEditor editor(*this);
    return editor.remove(segid);
}

/** Commit these changes and return a copy of the
    edited molecule */
Molecule MolEditor::commit() const
{
    return *this;
}

//////////
////////// Implementation of MolStructureEditor
//////////

static const RegisterMetaType<MolStructureEditor> r_molstructeditor;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const MolStructureEditor &moleditor)
{
    writeHeader(ds, r_molstructeditor, 1);
    
    ds << static_cast<const StructureEditor&>(moleditor);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       MolStructureEditor &moleditor)
{
    VersionID v = readHeader(ds, r_molstructeditor);
    
    if (v == 1)
    {
        ds >> static_cast<MolStructureEditor&>(moleditor);
    }
    else
        throw version_error( v, "1", r_molstructeditor, CODELOC );
        
    return ds;
}

/** Null constructor */
MolStructureEditor::MolStructureEditor() : StructureEditor()
{}

/** Construct a structure editor of the molecule viewed by 'molview' */
MolStructureEditor::MolStructureEditor(const MoleculeView &molview)
                   : StructureEditor(molview.data())
{}

/** Construct a structure editor that edits the structure of the molecule
    being edited by 'other' */
MolStructureEditor::MolStructureEditor(const StructureEditor &other)
                   : StructureEditor(other)
{}

/** Copy constructor */
MolStructureEditor::MolStructureEditor(const MolStructureEditor &other)
                   : StructureEditor(other)
{}

/** Destructor */
MolStructureEditor::~MolStructureEditor()
{}

/** Assign to edit the structure of a copy of the molecule viewed in 'molview' */
MolStructureEditor& MolStructureEditor::operator=(const MoleculeView &molview)
{
    StructureEditor::operator=(molview.data());
    return *this;
}

/** Copy assignment from another editor */
MolStructureEditor& MolStructureEditor::operator=(const StructureEditor &other)
{
    StructureEditor::operator=(other);
    return *this;
}

/** Copy assignment from another editor */
MolStructureEditor& MolStructureEditor::operator=(const MolStructureEditor &other)
{
    StructureEditor::operator=(other);
    return *this;
}

/** Return an editor for the atom at ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomStructureEditor MolStructureEditor::select(const AtomID &atomid)
{
    return AtomStructureEditor(*this, this->atomIdx(atomid));
}

/** Return an editor for the CutGroup at ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CGStructureEditor MolStructureEditor::select(const CGID &cgid)
{
    return CGStructureEditor(*this, this->cgIdx(cgid));
}

/** Return an editor for the residue at ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResStructureEditor MolStructureEditor::select(const ResID &resid)
{
    return ResStructureEditor(*this, this->resIdx(resid));
}

/** Return an editor for the chain at ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
ChainStructureEditor MolStructureEditor::select(const ChainID &chainid)
{
    return ChainStructureEditor(*this, this->chainIdx(chainid));
}

/** Return an editor for the segment at ID 'segid'

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
SegStructureEditor MolStructureEditor::select(const SegID &segid)
{
    return SegStructureEditor(*this, this->segIdx(segid));
}

/** Return an editor for the atom at ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomStructureEditor MolStructureEditor::atom(const AtomID &atomid)
{
    return AtomStructureEditor(*this, this->atomIdx(atomid));
}

/** Return an editor for the CutGroup at ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CGStructureEditor MolStructureEditor::cutGroup(const CGID &cgid)
{
    return CGStructureEditor(*this, this->cgIdx(cgid));
}

/** Return an editor for the residue at ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResStructureEditor MolStructureEditor::residue(const ResID &resid)
{
    return ResStructureEditor(*this, this->resIdx(resid));
}

/** Return an editor for the chain at ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
ChainStructureEditor MolStructureEditor::chain(const ChainID &chainid)
{
    return ChainStructureEditor(*this, this->chainIdx(chainid));
}

/** Return an editor for the segment at ID 'segid'

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
SegStructureEditor MolStructureEditor::segment(const SegID &segid)
{
    return SegStructureEditor(*this, this->segIdx(segid));
}

/** Return the name of this molecule */
const MolName& MolStructureEditor::name() const
{
    return this->molName();
}

/** Return this molecule's ID number */
MolNum MolStructureEditor::number() const
{
    return this->molNum();
}

/** Return the number of atoms in this molecule (may be zero!) */
int MolStructureEditor::nAtoms() const
{
    return this->nAtomsInMolecule();
}

/** Return the number of CutGroups in this molecule (may be zero!) */
int MolStructureEditor::nCutGroups() const
{
    return this->nCutGroupsInMolecule();
}

/** Return the number of residues in this molecule (may be zero!) */
int MolStructureEditor::nResidues() const
{
    return this->nResiduesInMolecule();
}

/** Return the number of chains in this molecule (may be zero!) */
int MolStructureEditor::nChains() const
{
    return this->nChainsInMolecule();
}

/** Return the number of segments in this molecule (may be zero!) */
int MolStructureEditor::nSegments() const
{
    return this->nSegmentsInMolecule();
}

/** Rename this molecule to 'newname' */
MolStructureEditor& MolStructureEditor::rename(const MolName &newname)
{
    this->renameMolecule(newname);
    return *this;
}

/** Give this molecule a new, unique ID number */
MolStructureEditor& MolStructureEditor::renumber()
{
    this->renumberMolecule();
    return *this;
}

/** Renumber this molecule to have the number 'newnum' */
MolStructureEditor& MolStructureEditor::renumber(MolNum newnum)
{
    this->renumberMolecule(newnum);
    return *this;
}

/** Add an atom called 'name' to this molecule and return an
    editor for this atom */
AtomStructureEditor MolStructureEditor::add(const AtomName &name)
{
    AtomStructureEditor editor = this->addAtom();
    editor.rename(name);
    return editor;
}

/** Add an atom with number 'number' to this molecule and 
    return an editor for this atom */
AtomStructureEditor MolStructureEditor::add(const AtomNum &number)
{
    AtomStructureEditor editor = this->addAtom();
    editor.renumber(number);
    return editor;
}

/** Add a residue called 'name' to this molecule and return an
    editor for this residue */
ResStructureEditor MolStructureEditor::add(const ResName &name)
{
    ResStructureEditor editor = this->addResidue();
    editor.rename(name);
    return editor;
}

/** Add a residue with number 'number' to this molecule and 
    return an editor for this residue */
ResStructureEditor MolStructureEditor::add(const ResNum &number)
{
    ResStructureEditor editor = this->addResidue();
    editor.renumber(number);
    return editor;
}

/** Add a CutGroup called 'name' to this molecule and return an
    editor for this CutGroup */
CGStructureEditor MolStructureEditor::add(const CGName &name)
{
    CGStructureEditor editor = this->addCutGroup();
    editor.rename(name);
    return editor;
}

/** Add a chain called 'name' to this molecule and return an
    editor for this chain */
ChainStructureEditor MolStructureEditor::add(const ChainName &name)
{
    ChainStructureEditor editor = this->addChain();
    editor.rename(name);
    return editor;
}

/** Add a segment called 'name' to this molecule and return an
    editor for this segment */
SegStructureEditor MolStructureEditor::add(const SegName &name)
{
    SegStructureEditor editor = this->addSegment();
    editor.rename(name);
    return editor;
}

/** Remove all atoms that match the ID 'atomid' - this does
    nothing if there are no atoms that match this ID */
MolStructureEditor& MolStructureEditor::remove(const AtomID &atomid)
{
    this->removeAtoms(atomid);
    return *this;
}

/** Remove all CutGroups that match the ID 'cgid' - this does
    nothing if there are no CutGroups that match this ID */
MolStructureEditor& MolStructureEditor::remove(const CGID &cgid)
{
    this->removeCutGroups(cgid);
    return *this;
}

/** Remove all residues that match the ID 'resid' - this does
    nothing if there are no residues that match this ID */
MolStructureEditor& MolStructureEditor::remove(const ResID &resid)
{
    this->removeResidues(resid);
    return *this;
}

/** Remove all chains that match the ID 'chainid' - this does
    nothing if there are no chains that match this ID */
MolStructureEditor& MolStructureEditor::remove(const ChainID &chainid)
{
    this->removeChains(chainid);
    return *this;
}

/** Remove all segments that match the ID 'segid' - this does
    nothing if there are no segments that match this ID */
MolStructureEditor& MolStructureEditor::remove(const SegID &segid)
{
    this->removeSegments(segid);
    return *this;
}

/** Commit the changes and return a Molecule copy */
Molecule MolStructureEditor::commit() const
{
    return Molecule( this->commitChanges() );
}

/** Automatically convert this editor to a Molecule() */
MolStructureEditor::operator Molecule() const
{
    return this->commit();
}
