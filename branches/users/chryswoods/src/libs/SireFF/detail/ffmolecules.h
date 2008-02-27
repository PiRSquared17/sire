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

#ifndef SIREFF_DETAIL_FFMOLECULES_H
#define SIREFF_DETAIL_FFMOLECULES_H

#include "partialmolecule.h"

SIRE_BEGIN_HEADER

namespace SireFF
{

namespace detail
{

/** This is the base class of most (all?) of the forcefield
    specialised molecules. Forcefields used specialised molecule
    classes as a way of storing important data in a way that
    is best for the forcefield itself. As these specialised
    classes are only used within a forcefield, they do not
    form part of the public API of Sire, and should not be
    used in any other code.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFMoleculeBase
{

friend QDataStream& ::operator<<(QDataStream&, const FFMolecule&);
friend QDataStream& ::operator>>(QDataStream&, FFMolecule&);

public:
    FFMolecule();
    FFMolecule(const PartialMolecule &molview);
    
    ~FFMolecule();
    
    const PartialMolecule& molecule() const;
    
    void change(const PartialMolecule &molecule);
    
private:
    /** Copy of the view that this object represents */
    PartialMolecule molview;
};

/** This is the forcefield specific part of the FFMolecule class
    that holds the forcefield parameters for this part of the molecule
    
    @author Christopher Woods
*/
template<class FFIELD>
class FFMolecule : public FFMoleculeBase
{
public:
    typedef FFIELD::Parameters Parameters;
    typedef FFIELD::ParameterNames ParameterNames;

    FFMolecule();
    
    FFMolecule(const PartialMolecule &molecule,
               FFIELD &forcefield,
               const ParameterNames &parameters);

    FFMolecule(const FFMolecule<FFIELD> &other);

    ~FFMolecule();

    FFMolecule<FFIELD>& operator=(const FFMolecule<FFIELD> &other);
    
    bool operator==(const FFMolecule<FFIELD> &other) const;
    bool operator!=(const FFMolecule<FFIELD> &other) const;

    const Parameters& parameters() const;
    
    void change(const PartialMolecule &molecule,
                FFIELD &forcefield,
                const ParameterNames &parameters);
private:
    /** The parameters for this molecule */
    Parameters molparams;
};

/** This is the base class of most (all?) of the 'Molecules'
    classes. These classes hold groups of molecules together
    in a way that allows very rapid iteration over them and 
    their forcefield parameters. By its nature, this is a 
    template class, so this provides the template independent
    parts.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFMoleculesBase
{
public:
    FFMoleculesBase();

    const QVector<MolNum> molNumsByIndex() const;
    const QHash<MolNum,quint32>& indexesByMolNum() const;

    quint32 indexOf(MolNum molnum) const;
    bool contains(MolNum molnum) const;

protected:
    void _pvt_removed(MolNum molnum);
    void _pvt_reindex();
    
private:
    /** The number of each molecule in the group, 
        in the order they appear in the arrays */
    QVector<MolNum> molnums_by_idx;
    
    /** The indicies of each molecule in the array
        indexed by molecule number */
    QHash<MolNum,quint32> idxs_by_molnum;
};

/** This class provides an array of FF specialised molecules, 
    arranged in such a way as to speed up indexing over them.
    
    The template parameter is the forcefield type for which 
    this group is used. This is necessary, as different
    forcefield types use different parameterisation methods,
    and can supply different specialised molecule types.
    
    @author Christopher Woods
*/
template<class FFIELD>
class FFMolecules : public FFMoleculesBase
{
public:
    typedef typename FFIELD::Molecule Molecule;
    typedef typename FFIELD::ChangedMolecule Molecule;

    typedef typename FFIELD::ParameterNames ParameterNames;

    FFMolecules();

    ChangedMolecule change(const PartialMolecule &molecule,
                           T &forcefield);
                           
    ChangedMolecule add(const PartialMolecule &molecule,
                        T &forcefield);
                        
    ChangedMolecule remove(const PartialMolecule &molecule,
                           T &forcefield);
    
    const QVector<Molecule>& moleculesByIndex() const;
    const QVector<ParameterNames>& parameterNamesByIndex() const;
    
protected:
    /** The array of forcefield-specialised molecules in this group */
    QVector<Molecule> mols_by_idx;
    
    /** The current parameter locations used for each molecule */
    QVector<ParameterNames> parameter_names;
};

/** This implementation class holds the information about 
    a change in a molecule. This can be used by a
    forcefield to simplify the calculation of the change in
    energy associated with a change in the molecule(s)
    (e.g. if only changes in energy associated with the 
    changed molecules needs to be evaluated) 
    
    The template parameter 'FFMOL' is the specialised Molecule
    type used by the forcefield. 
    
    @author Christopher Woods
*/
template<class FFMOL>
class ChangedMolecule
{
public:
    typedef FFMOL Molecule;

    ChangedMolecule();
    
    ChangedMolecule(const FFMOL &molecule);
    ChangedMolecule(const FFMOL &oldmol, const FFMOL &newmol);
    
    ChangedMolecule(const ChangedMolecule<FFMOL> &other);
    
    ~ChangedMolecule();
    
    ChangedMolecule<FFMOL>& operator=(const ChangedMolecule<FFMOL> &other);
    
    bool operator==(const ChangedMolecule<FFMOL> &other) const;
    bool operator!=(const ChangedMolecule<FFMOL> &other) const;
    
    bool isEmpty() const;
    
    bool changedAll() const;
    bool nothingChanged() const;
    
    const FFMOL& oldMolecule() const;
    const FFMOL& newMolecule() const;
    
    const FFMOL& oldParts() const;
    const FFMOL& newParts() const;
    
    ChangedMolecule<FFMOL>& change(const FFMOL &new_molecule);
    
private:
    /** The complete, old molecule */
    FFMOL old_molecule;
    
    /** The complete, new molecule */
    FFMOL new_molecule;
    
    /** The parts of the old molecule that have changed */
    FFMOL old_parts;
    
    /** The corresponding changed parts of the new molecule */
    FFMOL new_parts;
};

}

}

SIRE_END_HEADER

#endif
