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

#ifndef SIREMM_CLJDELTA_H
#define SIREMM_CLJDELTA_H

#include "cljatoms.h"
#include "cljboxes.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class CLJDelta;
}

QDataStream& operator<<(QDataStream&, const SireMM::CLJDelta&);
QDataStream& operator>>(QDataStream&, SireMM::CLJDelta&);

namespace SireMM
{

class CLJBoxes;

/** This class is used to hold the change in coordinates etc. of a set of atoms caused
    by e.g. a Monte Carlo move
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJDelta
{

friend QDataStream& ::operator<<(QDataStream&, const CLJDelta&);
friend QDataStream& ::operator>>(QDataStream&, CLJDelta&);

public:
    CLJDelta();
    CLJDelta(quint32 idnum, const CLJBoxes &boxes, const QVector<CLJBoxIndex> &old_atoms,
             const MoleculeView &new_atoms, const PropertyMap &map = PropertyMap());
    CLJDelta(quint32 idnum, const CLJBoxes &boxes, const QVector<CLJBoxIndex> &old_atoms,
             const MoleculeView &new_atoms, CLJAtoms::ID_SOURCE source,
             const PropertyMap &map = PropertyMap());
    
    CLJDelta(const CLJDelta &other);
    
    ~CLJDelta();
    
    CLJDelta& operator=(const CLJDelta &other);
    
    bool operator==(const CLJDelta &other) const;
    bool operator!=(const CLJDelta &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    void reconstruct(quint32 idnum, const CLJBoxes &boxes, const QVector<CLJBoxIndex> &old_atoms,
                     const MoleculeView &new_atoms, const PropertyMap &map = PropertyMap());
    
    void reconstruct(quint32 idnum, const CLJBoxes &boxes, const QVector<CLJBoxIndex> &old_atoms,
                     const MoleculeView &new_atoms, CLJAtoms::ID_SOURCE source,
                     const PropertyMap &map = PropertyMap());

    QString toString() const;
    
    bool isSingleBox() const;
    
    const CLJBoxIndex& boxIndex() const;
    
    CLJBoxIndex minBox() const;
    CLJBoxIndex maxBox() const;
    
    bool isNull() const;
    
    quint32 ID() const;
    
    quint8 nBoxX() const;
    quint8 nBoxY() const;
    quint8 nBoxZ() const;
    
    int nBoxes() const;
    
    const QVector<CLJBoxIndex>& oldIndicies() const;
    
    CLJAtoms oldAtoms() const;
    
    const CLJAtoms& changedAtoms() const;
    
    const CLJAtoms& newAtoms() const;
    
    float boxLength() const;
    
    static CLJDelta merge(const CLJDelta *deltas, int count, bool changes_only=false);
    
private:
    friend class CLJBoxes;
    void buildFrom(const CLJBoxes &boxes, const QVector<CLJBoxIndex> &indicies,
                   const MoleculeView &new_molecule, CLJAtoms::ID_SOURCE source,
                   const PropertyMap &map);

    /** The new atoms */
    CLJAtoms new_atoms;

    /** The combined old and new atoms */
    CLJAtoms changed_atoms;
    
    /** The CLJBoxIndex of the lowest box of the changed atoms */
    CLJBoxIndex box_index;
    
    /** The indicies of the old atoms in the CLJBoxes object */
    QVector<CLJBoxIndex> old_indicies;
    
    /** The number of boxes in the x, y, and z dimensions
        that these changed atoms encompass */
    quint8 nbox_x, nbox_y, nbox_z;
    
    /** Whether or not the changed atoms fit into a single box */
    quint8 is_single_box;
    
    /** The ID number of the molecule or part of molecule that this
        delta represents. This is used for book-keeping */
    quint32 idnum;
    
    /** The box length used by this delta */
    float box_length;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return whether or not the CLJAtoms are all contained in a single box */
inline bool CLJDelta::isSingleBox() const
{
    return is_single_box;
}

/** Return the index of the box that holds the lowest index of the atoms */
inline const CLJBoxIndex& CLJDelta::boxIndex() const
{
    return box_index;
}

/** Return whether or not this delta is null */
inline bool CLJDelta::isNull() const
{
    return nbox_x == 0;
}

/** Return the ID number of this delta. This can be used for book-keeping
    by the object that created this delta */
inline quint32 CLJDelta::ID() const
{
    return idnum;
}

/** Return the number of boxes along the X dimension that the changed
    atoms cover */
inline quint8 CLJDelta::nBoxX() const
{
    return nbox_x;
}

/** Return the number of boxes along the X dimension that the changed
    atoms cover */
inline quint8 CLJDelta::nBoxY() const
{
    return nbox_y;
}

/** Return the number of boxes along the X dimension that the changed
    atoms cover */
inline quint8 CLJDelta::nBoxZ() const
{
    return nbox_z;
}

/** Return the of the old atoms in the CLJBoxes used to construct this object */
inline const QVector<CLJBoxIndex>& CLJDelta::oldIndicies() const
{
    return old_indicies;
}

/** Return the set of changed atoms (has a negated version of the old
    atoms so that we can calculate energy deltas directly) */
inline const CLJAtoms& CLJDelta::changedAtoms() const
{
    return changed_atoms;
}

/** Return the new version of the atoms in a format that should make
    them each to rebox */
inline const CLJAtoms& CLJDelta::newAtoms() const
{
    return new_atoms;
}

/** Return the length of the box (in angstroms) */
inline float CLJDelta::boxLength() const
{
    return box_length;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMM::CLJDelta )

SIRE_EXPOSE_CLASS( SireMM::CLJDelta )

SIRE_END_HEADER

#endif
