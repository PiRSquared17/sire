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

#include "cljdelta.h"

#include "SireMol/moleculeview.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<CLJDelta> r_delta(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJDelta &delta)
{
    writeHeader(ds, r_delta, 1);
    
    SharedDataStream sds(ds);
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJDelta &delta)
{
    VersionID v = readHeader(ds, r_delta);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
    }
    else
        throw version_error(v, "1", r_delta, CODELOC);
    
    return ds;
}

/** Construct from the passed data, but re-using the memory already allocated
    for this delta */
void CLJDelta::buildFrom(const CLJBoxes &boxes, const QVector<CLJBoxIndex> &indicies,
                         const MoleculeView &new_molecule, CLJAtoms::ID_SOURCE source,
                         const PropertyMap &map)
{
    CLJBoxIndex min_box;
    CLJBoxIndex max_box;

    //how many changed atoms are there?
    const int nold = CLJBoxIndex::countNonDummies(indicies);
    
    //construct the new molecule
    new_atoms.reconstruct(new_molecule, source, map);

    const int nnew = new_atoms.nAtoms();

    if (nnew > 0)
    {
        //box up each of the new atoms
        const MultiFloat *x = new_atoms.x().constData();
        const MultiFloat *y = new_atoms.y().constData();
        const MultiFloat *z = new_atoms.z().constData();
        const MultiInt *id = new_atoms.ID().constData();
        
        const qint32 id_of_dummy = CLJAtoms::idOfDummy()[0];
        
        const float inv_length = 1.0f / boxes.length().value();
        
        for (int i=0; i<new_atoms.x().count(); ++i)
        {
            for (int j=0; j<MultiFloat::count(); ++j)
            {
                if (id[i][j] != id_of_dummy)
                {
                    CLJBoxIndex idx = CLJBoxIndex::createWithInverseBoxLength(
                                                x[i][j], y[i][j], z[i][j], inv_length);
                    
                    if (min_box.isNull())
                    {
                        min_box = idx;
                        max_box = idx;
                    }
                    else
                    {
                        min_box = min_box.min(idx);
                        max_box = max_box.max(idx);
                    }
                }
            }
        }
    }

    changed_atoms.resize(nold + nnew);
    changed_atoms.copyIn(new_atoms);

    //now grab a copy of all of the old atoms, putting the negative
    //into changed_atoms and working out the range of boxes covered by
    //these atoms
    if (not indicies.isEmpty())
    {
        int old_idx = nnew;
    
        //find the range of boxes occupied by the atoms and copy out the data
        const CLJBoxIndex *indicies_array = indicies.constData();

        for (int i=0; i<indicies.count(); ++i)
        {
            const CLJBoxIndex &index = indicies_array[i];
            
            if (not index.isNull())
            {
            
                if (min_box.isNull())
                {
                    min_box = index;
                    max_box = min_box;
                }
                else
                {
                    min_box = min_box.min(index);
                    max_box = max_box.max(index);
                }
                
                //copy in the data from the ith old atom
                CLJAtom old_atom = boxes[index];
                
                if (not old_atom.isDummy())
                {
                    changed_atoms.set(old_idx, old_atom.negate());
                    old_idx += 1;
                }
            }
        }
    }

    old_indicies = indicies;
    
    //now work out the range of boxes covered by these atoms
    box_index = min_box.boxOnly();
    nbox_x = max_box.i() - min_box.i() + 1;
    nbox_y = max_box.j() - min_box.j() + 1;
    nbox_z = max_box.k() - min_box.k() + 1;
    is_single_box = (max_box.i() == min_box.i()) and
                    (max_box.j() == min_box.j()) and
                    (max_box.k() == min_box.k());
    
    //save the box length
    box_length = boxes.length().value();
}

/** Null constructor */
CLJDelta::CLJDelta() : nbox_x(0), nbox_y(0), nbox_z(0), is_single_box(0), idnum(0), box_length(0)
{}

/** Construct the delta that changes from the atoms at indicies 'old_atoms' in the 
    CLJBoxes 'boxes' to the atoms in view 'new_atoms', using the optionally supplied
    property map to extract the necessary properties */
CLJDelta::CLJDelta(quint32 id, const CLJBoxes &boxes, const QVector<CLJBoxIndex> &old_atoms,
                   const MoleculeView &new_atoms, const PropertyMap &map)
         : nbox_x(0), nbox_y(0), nbox_z(0), is_single_box(0), idnum(id), box_length(0)
{
    this->buildFrom(boxes, old_atoms, new_atoms, CLJAtoms::USE_MOLNUM, map);
}

/** Construct the delta that changes from the atoms at indicies 'old_atoms' in the 
    CLJBoxes 'boxes' to the atoms in view 'new_atoms', using the optionally supplied
    property map to extract the necessary properties */
void CLJDelta::reconstruct(quint32 id, const CLJBoxes &boxes, const QVector<CLJBoxIndex> &old_atoms,
                           const MoleculeView &new_atoms, const PropertyMap &map)
{
    nbox_x = 0;
    nbox_y = 0;
    nbox_z = 0;
    is_single_box = 0;
    idnum = id;
    box_length = 0;
    this->buildFrom(boxes, old_atoms, new_atoms, CLJAtoms::USE_MOLNUM, map);
}

/** Construct the delta that changes from the atoms at indicies 'old_atoms' in the 
    CLJBoxes 'boxes' to the atoms in view 'new_atoms', using the optionally supplied
    property map to extract the necessary properties */
CLJDelta::CLJDelta(quint32 id, const CLJBoxes &boxes, const QVector<CLJBoxIndex> &old_atoms,
                   const MoleculeView &new_atoms, CLJAtoms::ID_SOURCE source,
                   const PropertyMap &map)
         : nbox_x(0), nbox_y(0), nbox_z(0), is_single_box(0), idnum(id), box_length(0)
{
    this->buildFrom(boxes, old_atoms, new_atoms, source, map);
}

/** Construct the delta that changes from the atoms at indicies 'old_atoms' in the 
    CLJBoxes 'boxes' to the atoms in view 'new_atoms', using the optionally supplied
    property map to extract the necessary properties */
void CLJDelta::reconstruct(quint32 id, const CLJBoxes &boxes, const QVector<CLJBoxIndex> &old_atoms,
                           const MoleculeView &new_atoms, CLJAtoms::ID_SOURCE source,
                           const PropertyMap &map)
{
    nbox_x = 0;
    nbox_y = 0;
    nbox_z = 0;
    is_single_box = 0;
    idnum = id;
    box_length = 0;
    this->buildFrom(boxes, old_atoms, new_atoms, source, map);
}

/** Copy constructor */
CLJDelta::CLJDelta(const CLJDelta &other)
         : new_atoms(other.new_atoms), changed_atoms(other.changed_atoms),
           box_index(other.box_index), old_indicies(other.old_indicies),
           nbox_x(other.nbox_x), nbox_y(other.nbox_y), nbox_z(other.nbox_z),
           is_single_box(other.is_single_box), idnum(other.idnum), box_length(other.box_length)
{}

/** Destructor */
CLJDelta::~CLJDelta()
{}

/** Copy assignment operator */
CLJDelta& CLJDelta::operator=(const CLJDelta &other)
{
    if (this != &other)
    {
        new_atoms = other.new_atoms;
        changed_atoms = other.changed_atoms;
        box_index = other.box_index;
        old_indicies = other.old_indicies;
        nbox_x = other.nbox_x;
        nbox_y = other.nbox_y;
        nbox_z = other.nbox_z;
        is_single_box = other.is_single_box;
        idnum = other.idnum;
        box_length = other.box_length;
    }
    
    return *this;
}

/** Comparison operator */
bool CLJDelta::operator==(const CLJDelta &other) const
{
    return this == &other or
           (new_atoms == other.new_atoms and changed_atoms == other.changed_atoms and
            box_index == other.box_index and old_indicies == other.old_indicies and
            nbox_x == other.nbox_x and nbox_y == other.nbox_y and nbox_z == other.nbox_z and
            is_single_box == other.is_single_box and idnum == other.idnum and
            box_length == other.box_length);
}

/** Comparison operator */
bool CLJDelta::operator!=(const CLJDelta &other) const
{
    return not operator==(other);
}

const char* CLJDelta::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJDelta>() );
}

const char* CLJDelta::what() const
{
    return CLJDelta::typeName();
}

QString CLJDelta::toString() const
{
    return QObject::tr("CLJDelta( nChanged() = %1, nBoxes() = %2 )")
                .arg(changedAtoms().count()).arg(nBoxes());
}

/** Return the total number of boxes covered by these atoms */
int CLJDelta::nBoxes() const
{
    return int(nbox_x) * int(nbox_y) * int(nbox_z);
}

/** Return the index of the minimum box occupied by atoms of this delta */
CLJBoxIndex CLJDelta::minBox() const
{
    return box_index.boxOnly();
}

/** Return the index of the maximum box occupied by atoms of this delta */
CLJBoxIndex CLJDelta::maxBox() const
{
    return CLJBoxIndex( box_index.i() + nbox_x - 1,
                        box_index.j() + nbox_y - 1,
                        box_index.k() + nbox_z - 1 );
}

/** Return the old version of the changed atoms */
CLJAtoms CLJDelta::oldAtoms() const
{
    //the old atoms are the ones after the new atoms in 'changed_atoms'
    if (changed_atoms.isEmpty())
        return QVector<CLJAtom>();
    
    const int nold = CLJBoxIndex::countNonDummies(old_indicies);
    
    if (nold == 0)
        return QVector<CLJAtom>();
    
    const int nnew = new_atoms.nAtoms();
    
    if (changed_atoms.nAtoms() - nnew != nold)
        throw SireError::program_bug( QObject::tr(
                "Something wrong? %1 - %2 != %3")
                    .arg(changed_atoms.nAtoms())
                    .arg(nnew)
                    .arg(nold), CODELOC );
    
    CLJAtoms old_atoms;
    old_atoms.resize(nold);
    
    for (int i=0; i<nold; ++i)
    {
        old_atoms.set(i, changed_atoms[nnew + i].negate());
    }

    return old_atoms;
}

/** Merge together 'n' deltas from the passed array into a single delta.
    If 'changes_only' is true, then this will only merge together the changed
    atoms with box and will ignore the bookkeeping data. The resulting delta
    would thus only be useful for delta energy calculations */
CLJDelta CLJDelta::merge(const CLJDelta *deltas, int n, bool changes_only)
{
    if (n == 0)
        return CLJDelta();
    
    else if (n == 1)
        return deltas[0];

    else if (changes_only)
    {
        CLJDelta ret = deltas[0];
        ret.new_atoms = CLJAtoms();
        ret.old_indicies = QVector<CLJBoxIndex>();
        ret.idnum = 0;
        
        CLJBoxIndex min_box = ret.minBox();
        CLJBoxIndex max_box = ret.maxBox();
        
        for (int i=1; i<n; ++i)
        {
            const CLJDelta &delta = deltas[i];
            
            if (delta.box_length != ret.box_length)
                throw SireError::incompatible_error( QObject::tr(
                        "You cannot merge together deltas that have different box lengths! "
                        "%1 A vs. %2 A").arg(ret.box_length).arg(delta.box_length), CODELOC );
            
            ret.changed_atoms += delta.changed_atoms;
            min_box = min_box.min( delta.minBox() );
            max_box = max_box.max( delta.maxBox() );
        }
        
        ret.box_index = min_box.boxOnly();
        ret.nbox_x = max_box.i() - min_box.i() + 1;
        ret.nbox_y = max_box.j() - min_box.j() + 1;
        ret.nbox_z = max_box.k() - min_box.k() + 1;
        ret.is_single_box = (max_box.i() == min_box.i()) and
                            (max_box.j() == min_box.j()) and
                            (max_box.k() == min_box.k());
        
        return ret;
    }
    else
    {
        CLJDelta ret = deltas[0];
        
        CLJBoxIndex min_box = ret.minBox();
        CLJBoxIndex max_box = ret.maxBox();
        
        for (int i=1; i<n; ++i)
        {
            const CLJDelta &delta = deltas[i];
            
            if (delta.box_length != ret.box_length)
                throw SireError::incompatible_error( QObject::tr(
                        "You cannot merge together deltas that have different box lengths! "
                        "%1 A vs. %2 A").arg(ret.box_length).arg(delta.box_length), CODELOC );

            if (delta.ID() != ret.idnum)
                throw SireError::incompatible_error( QObject::tr(
                        "You cannot merge together deltas that have different ID numbers! "
                        "%1 vs. %2.").arg(ret.ID()).arg(delta.ID()), CODELOC );
            
            ret.new_atoms += delta.new_atoms;
            ret.changed_atoms += delta.changed_atoms;
            ret.old_indicies += delta.old_indicies;
            
            min_box = min_box.min( delta.minBox() );
            max_box = max_box.max( delta.maxBox() );
        }
        
        ret.box_index = min_box.boxOnly();
        ret.nbox_x = max_box.i() - min_box.i() + 1;
        ret.nbox_y = max_box.j() - min_box.j() + 1;
        ret.nbox_z = max_box.k() - min_box.k() + 1;
        ret.is_single_box = (max_box.i() == min_box.i()) and
                            (max_box.j() == min_box.j()) and
                            (max_box.k() == min_box.k());
        
        return ret;
    }
}
