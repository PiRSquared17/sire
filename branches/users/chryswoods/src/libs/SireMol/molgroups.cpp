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

#include "molgroups.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireStream;

/////////////
///////////// Implementation of MolGroupsBase
/////////////

RegisterMetaType<MolGroupsBase> r_molgroupsbase(MAGIC_ONLY,
                                                "SireMol::MolGroupsBase");
                                                
/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const MolGroupsBase &molgroupsbase)
{
    writeHeader(ds, r_molgroupsbase, 1);
    
    SharedDataStream sds(ds);
    
    sds << molgroupsbase.mgidx_to_num
        << molgroupsbase.mgname_to_mgnum
        << molgroupsbase.molnum_to_mgnum
        << static_cast<const PropertyBase&>(molgroupsbase);
        
    return ds;
}

/** Read from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       MolGroupsBase &molgroupsbase)
{
    VersionID v = readHeader(ds, r_molgroupsbase);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> molgroupsbase.mgidx_to_num
            >> molgroupsbase.mgname_to_mgnum
            >> molgroupsbase.molnum_to_mgnum
            >> static_cast<PropertyBase&>(molgroupsbase);
    }
    else
        throw version_error(r_molgroupsbase, "1", CODELOC);
        
    return ds;
}

/** Null constructor */
MolGroupsBase::MolGroupsBase() : PropertyBase()
{}

/** Copy constructor */
MolGroupsBase::MolGroupsBase(const MolGroupsBase &other)
              : PropertyBase(other),
                mgidx_to_num(other.mgidx_to_num),
                mgname_to_mgnum(other.mgname_to_mgnum),
                molnum_to_mgnum(other.molnum_to_mgnum)
{}

MolGroupsBase::~MolGroupsBase();

/** Copy assignment operator */
MolGroupsBase& MolGroupsBase::operator=(const MolGroupsBase &other)
{
    if (this != &other)
    {
        mgidx_to_num = other.mgidx_to_num;
        mgname_to_mgnum = other.mgname_to_mgnum;
        molnum_to_mgnum = other.molnum_to_mgnum;
    }
    
    return *this;
}

/** Return a const reference to the molecule group with number 'mgnum'

    \throw SireMol::missing_group
*/
const MolGroup& MolGroupsBase::operator[](MGNum mgnum) const
{
    return this->at(mgnum);
}

/** Return a const reference to the molecule group with name 'mgname'

    \throw SireMol::missing_group
*/
const MolGroup& MolGroupsBase::operator[](const MGName &mgname) const
{
    return this->at(mgname);
}

/** Return a const reference to the molecule group that is 
    identified by 'mgid'
    
    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
*/
const MolGroup& MolGroupsBase::operator[](const MGID &mgid) const
{
    return this->at(mgid);
}

/** Return all of the views of the molecule with number 'molnum'
    that are contained in this set of groups. Note that if the 
    same view appears in multiple groups, then it will be returned
    multiple times in the returned set of views
    
    \throw SireMol::missing_molecule
*/
ViewsOfMol MolGroupsBase::operator[](MolNum molnum) const
{
    return this->at(molnum);
}

/** Return all of the views of the molecule identified by 'molid'
    that are contained in this set of groups. Note that if the 
    same view appears in multiple groups, then it will be returned
    multiple times in the returned set of views
    
    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
ViewsOfMol MolGroupsBase::operator[](const MolID &molid) const
{
    return this->at(molid);
}

Segment MolGroupsBase::operator[](const SegID &segid) const
{
    return this->at(segid);
}

Chain MolGroupsBase::operator[](const ChainID &chainid) const
{
    return this->at(chainid);
}

Residue MolGroupsBase::operator[](const ResID &resid) const
{
    return this->at(resid);
}

CutGroup MolGroupsBase::operator[](const CGID &cgid) const
{
    return this->at(cgid);
}

Atom MolGroupsBase::operator[](const AtomID &atomid) const
{
    return this->at(atomid);
}

/** Get the number of the molecule group whose number is 'mgnum'.
    This is an obvious function, only provided as a shortcut
    to prevent the MGID function being called if an MGNum is passed.
    
    \throw SireMol::missing_group
*/
MGNum MolGroupsBase::getGroupNumber(MGNum mgnum) const
{
    this->assertContains(mgnum);
    return mgnum;
}

/** Return the number of the molecule group that is called 'mgname'.

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
*/
MGNum MolGroupsBase::getGroupNumber(const MGName &mgname) const
{
    QList<MGNum> mgnums = this->map(mgname);
    
    if (mgnums.count() > 1)
        throw SireMol::duplicate_group( QObject::tr(
            "There are multiple molecule groups called \"%1\" "
            "in this set - with group numbers %2.")
                .arg(mgname).arg(Sire::toString(mgnums)), CODELOC );
                
    return mgnums.first();
}

/** Return the number of the group at index 'mgidx'

    \throw SireError::invalid_index
*/
MGNum MolGroupsBase::getGroupNumber(MGIdx mgidx) const
{
    return mgidx_to_num.constData()[ mgidx.map(mgidx_to_num.count()) ];
}

/** Return the number of the groups that matches 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
MGNum MolGroupsBase::getGroupNumber(const MGID &mgid) const
{
    QList<MGNum> mgnums = this->map(mgid);
    
    if (mgnums.count() > 1)
        throw SireMol::duplicate_group( QObject::tr(
            "There is more than one molecule group that matches "
            "the ID \"%1\". Matching groups have numbers %2.")
                .arg(mgid.toString()).arg(Sire::toString(mgnums)), 
                    CODELOC );
                    
    return mgnums.first();
}

/** Return the numbers of all groups in this set that are called
    'mgname'
    
    \throw SireMol::missing_group
*/
QList<MGNum> MolGroups::map(const MGName &mgname) const
{
    QHash< MGName,QList<MGNum> >::const_iterator it = mgname_to_mgnum.find(mgname);
    
    if (it == mgname_to_mgnum.end())
        throw SireMol::missing_group( QObject::tr(
            "There are no molecule groups called \"%1\" in this set. "
            "Available groups are %2.")
                .arg(mgname).arg(Sire::toString(mgname_to_mgnum.keys())), 
                    CODELOC );
                    
    return *it;
}

/** Return the list of numbers of groups that have the number 'mgnum'.
    This is a simple and obvious function that acts as a shortcut
    preventing map(const MGID&) being called for an MGNum
    
    \throw SireMol::missing_group
*/
QList<MGNum> MolGroupsBase::map(MGNum mgnum) const
{
    this->assertContains(mgnum);
    
    QList<MGNum> mgnums;
    mgnums.append(mgnum);
    
    return mgnums;
}

/** Return the list (of only one) molecule group that is at  
    index 'mgidx'
    
    \throw SireError::invalid_index
*/
QList<MGNum> MolGroupsBase::map(MGIdx mgidx) const
{
    QList<MGNum> mgnums;
    
    mgnums.append( this->getGroupNumber(mgidx) );
    
    return mgnums;
}

/** Map the molecule group ID 'mgid' to the list of molecule
    group numbers of the groups that match this ID in this set.
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
QList<MGNum> MolGroupsBase::map(const MGID &mgid) const
{
    return mgid.map(*this);
}

/** Simple function that just checks if a molecule with number 
    'molnum' is in the set, and returns it. This shortcuts
    the getMoleculeNumber(const MolID&) function in the case
    of MolNums
    
    \throw SireMol::missing_molecule
*/
MolNum MolGroupsBase::getMoleculeNumber(MolNum molnum) const
{
    this->assertContains(molnum);
    return molnum;
}

/** Return the number of the molecule at index 'molidx' in 
    this set
    
    \throw SireError::invalid_index
*/
MolNum MolGroupsBase::getMoleculeNumber(MolIdx molidx) const
{
    int n = molidx.map( molnum_to_mgnum.count() );
    
    //find the nth molnum in the hash - this is potentially very slow!
    QHash< MolNum,QList<MGNum> >::const_iterator it = molnum_to_mgnum.begin();
    
    for (int i=0; i<n; ++i)
    {
        ++it;
    }
    
    return it.key();
}

/** Return the number of the molecule called 'molname' from this set.
    
    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
MolNum MolGroupsBase::getMoleculeNumber(const MolName &molname) const
{
    QList<MolNum> molnums = this->map(molname);
    
    if (molnums.count() > 1)
        throw SireMol::duplicate_molecule( QObject::tr(
            "There is more than one molecule with the name \"%1\" "
            "in this set of molecule groups. Molecules with this "
            "name have numbers %2.")
                .arg(molname).arg(Sire::toString(molnums)), CODELOC );
                
    return molnums.first();
}

/** Return the number of the molecule that matches the ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
MolNum MolGroupsBase::getMoleculeNumber(const MolID &molid) const
{
    QList<MolNum> molnums = this->map(molid);
    
    if (molnums.count() > 1)
        throw SireMol::duplicate_molecule( QObject::tr(
            "There is more than one molecule that matches the "
            "ID \"%1\". Molecules that match have number %2.")
                .arg(molid.toString()).arg(Sire::toString(molnums)),
                    CODELOC );
                    
    return molnums.first();
}

/** Simple function that provides a shortcut for map(const MolID&)

    \throw SireMol::missing_molecule
*/
QList<MolNum> MolGroupsBase::map(MolNum molnum) const
{
    this->assertContains(molnum);
    
    QList<MolNum> molnums;
    molnums.append(molnum);
    
    return molnums;
}

/** Return the number of the molecule at index 'molidx'

    \throw SireError::invalid_index
*/
QList<MolNum> MolGroupsBase::map(MolIdx molidx) const
{
    QList<MolNum> molnums;
    molnums.append( this->getMoleculeNumber(molidx) );
    
    return molnums;
}

/** Return the numbers of all of the molecules that have the
    name 'molname'
    
    \throw SireMol::missing_molecule
*/
QList<MolNum> MolGroupsBase::map(const MolName &molname) const
{
    QList<MolNum> molnums;
    
    //loop over all of the groups in this set
    const QHash<MGNum,const MolGroup*> molgroups = this->getGroups();
    
    for (QHash<MGNum,const MolGroup*>::const_iterator it = molgroups.constBegin();
         it != molgroups.constEnd();
         ++it)
    {
        try
        {
            molnums.append( (*it)->map(molname) );
        }
        catch(...)
        {}
    }
    
    if (molnums.isEmpty())
        throw SireMol::missing_molecule( QObject::tr(
            "There are no molecules called \"%1\" in any of the "
            "molecule groups in this set.")
                .arg(molname), CODELOC );
    
    else if (molnums.count() > 1)
        //remove any duplicates
        molnums = molnums.toSet().toList();

    return molnums;
}

/** Return the numbers of all molecules that match the ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
QList<MolNum> MolGroupsBase::map(const MolID &molid) const
{
    return molid.map(*this);
}

/** Return a const reference to the group at index 'mgidx'

    \throw SireError::invalid_index
*/
const MolGroup& MolGroupsBase::at(MGIdx mgidx) const
{
    return this->at( this->getGroupNumber(mgidx) );
}

/** Return a const reference to the group in this set 
    called 'mgname'
    
    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
*/
const MolGroup& MolGroupsBase::at(const MGName &mgname) const
{
    return this->at( this->getGroupNumber(mgname) );
}

/** Return a const reference to the group in this set that
    is identified by 'mgid' 
    
    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
    \throw SireError::invalid_index
*/
const MolGroup& MolGroupsBase::at(const MGID &mgid) const
{
    return this->at( this->getGroupNumber(mgid) );
}

/** Return all of the views of the molecule with number 'molnum'
    that are contained in this set of groups. Note that if the 
    same view appears in multiple groups, then it will be returned
    multiple times in the returned set of views
    
    \throw SireMol::missing_molecule
*/
ViewsOfMol MolGroupsBase::at(MolNum molnum) const
{
    //get the list of groups that contain this molecule
    const QList<MGNum> &mgnums = this->groupsContaining(molnum);
    
    if (mgnums.count() == 1)
        return this->at( mgnums.first() ).molecule(molnum);
    else
    {
        ViewsOfMol mol = this->at( mgnums.first() ).molecule(molnum);
        
        for (int i=1; i<mgnums.count(); ++i)ß
        {
            mol.add( this->at( mgnums.at(i) ).molecule(molnum) );
        }
        
        return mol;
    }
}

/** Return all of the views of the molecule identified by 'molid'
    that are contained in this set of groups. Note that if the 
    same view appears in multiple groups, then it will be returned
    multiple times in the returned set of views
    
    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
ViewsOfMol MolGroupsBase::at(const MolID &molid) const
{
    QList<MolNum> molnums = molid.map(*this);
    
    if (molnums.count() > 1)
        throw SireMol::duplicate_molecule( QObject::tr(
            "There is more than one molecule that matches the ID "
            "\"%1\". Matching molecules have numbers %2.")
                .arg(molid.toString())
                .arg(Sire::toString(molnums)), CODELOC );
                
    return this->at(molnums.first());
}

/** Return the segment from this set that matches the ID 'segid'.
    This segment must be wholly contained by one of the groups
    in this set
    
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
Segment MolGroupsBase::at(const SegID &segid) const
{
    return segid.selectFrom(*this);
}

/** Return the chain from this set that matches the ID 'chainid'.
    This chain must be wholly contained by one of the groups
    in this set
    
    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
Chain MolGroupsBase::at(const ChainID &chainid) const
{
    return chainid.selectFrom(*this);
}

/** Return the residue from this set that matches the ID 'resid'.
    This residue must be wholly contained by one of the groups
    in this set
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
Residue MolGroupsBase::at(const ResID &resid) const
{
    return resid.selectFrom(*this);
}

/** Return the CutGroup from this set that matches the ID 'cgid'.
    This CutGroup must be wholly contained by one of the groups
    in this set
    
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CutGroup MolGroupsBase::at(const CGID &cgid) const
{
    return cgid.selectFrom(*this);
}

/** Return the atom from this set that matches the ID 'atomid'.
    This atom must be contained in one of the groups in this set.
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Atom MolGroupsBase::at(const AtomID &atomid) const
{
    return atomid.selectFrom(*this);
}

/** Return the MoleculeGroup that matches the ID 'mgid'

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
    \throw SireError::invalid_index
*/
MoleculeGroup MolGroupsBase::select(const MGID &mgid) const
{
    return this->at(mgid);
}

/** Return all of the views of the molecule with number 'molnum'
    that are contained in this set of groups. Note that if the 
    same view appears in multiple groups, then it will be returned
    multiple times in the returned set of views
    
    \throw SireMol::missing_molecule
*/
ViewsOfMol MolGroupsBase::select(const MolID &molid) const
{
    return this->at(molid);
}

/** Return the segment from this set that matches the ID 'segid'.
    This segment must be wholly contained by one of the groups
    in this set
    
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
Segment MolGroupsBase::select(const SegID &segid) const
{
    return this->at(segid);
}

/** Return the chain from this set that matches the ID 'chainid'.
    This chain must be wholly contained by one of the groups
    in this set
    
    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
Chain MolGroupsBase::select(const ChainID &chainid) const
{
    return this->at(chainid);
}

/** Return the residue from this set that matches the ID 'resid'.
    This residue must be wholly contained by one of the groups
    in this set
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
Residue MolGroupsBase::select(const ResID &resid) const
{
    return this->at(resid);
}

/** Return the CutGroup from this set that matches the ID 'cgid'.
    This CutGroup must be wholly contained by one of the groups
    in this set
    
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CutGroup MolGroupsBase::select(const CGID &cgid) const
{
    return this->at(cgid);
}

/** Return the atom from this set that matches the ID 'atomid'.
    This atom must be contained in one of the groups in this set.
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Atom MolGroupsBase::select(const AtomID &atomid) const
{
    return this->at(atomid);
}

/** Obvious shortcut for select(const MGID&) 
    
    \throw SireMol::missing_group
*/
QList<MoleculeGroup> MolGroupsBase::selectAll(MGNum mgnum) const
{
    QList<MoleculeGroup> molgroups;
    molgroups.append( this->at(mgnum) );
    
    return molgroups;
}

/** Obvious shortcut for select(const MGID&) 

    \throw SireError::invalid_index
*/
QList<MoleculeGroup> MolGroupsBase::selectAll(MGIdx mgidx) const
{
    QList<MoleculeGroup> molgroups;
    molgroups.append( this->at(mgidx) );
    
    return molgroups;
}

/** Return all of the molecule groups that are called 'mgname'

    \throw SireMol::missing_group
*/
QList<MoleculeGroup> MolGroupsBase::selectAll(const MGName &mgname) const
{
    QHash< MGName,QList<MGNum> >::const_iterator it = mgname_to_mgnum.find(mgname);
    
    if (it == mgname_to_mgnum.end())
        throw SireMol::missing_group( QObject::tr(
            "There are no groups in this set called \"%1\". "
            "The groups in this set are called %2.")
                .arg(mgname).arg(Sire::toString(mgname_to_mgnum.keys())),
                    CODELOC );
          
    //now get all of the groups          
    QVarLengthArray<const MolGroup*,10> groups;
    this->getGroups(*it, groups);
    
    QList<MoleculeGroup> molgroups;
    
    foreach (const MolGroup *group, groups)
    {
        molgroups.append( MoleculeGroup(*group) );
    }
    
    return molgroups;
}

/** Return all of the molecule groups that match the ID 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
QList<MoleculeGroup> MolGroupsBase::selectAll(const MGID &mgid) const
{
    //get the list of numbers that match this ID
    QList<MGNum> mgnums = this->map(mgid);
    
    //now pick up those groups...
    QVarLengthArray<const MolGroup*,10> groups;
    this->getGroups(mgnums, groups);

    QList<MoleculeGroup> molgroups;
    
    foreach (const MolGroup *group, groups)
    {
        molgroups.append( MoleculeGroup(*group) );
    }
    
    return molgroups;
}

/** Return the views of the molecule(s) that match the molecule ID
    'molid'. This returns all views of the molecule in the groups,
    and if a view is contained multiple times, then multiple copies
    of that view will be returned.
    
    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
QList<ViewsOfMol> MolGroupsBase::selectAll(const MolID &molid) const
{
    //get the numbers of molecules that match this ID
    QList<MolNum> molnums = this->map(molid);
    
    QList<ViewsOfMol> molviews;
    
    foreach (MolNum molnum, molnums)
    {
        molviews.append(this->at(molnum));
    }
    
    return molviews;
}

/** Return all of the segments from this set that match the ID 'segid'.
    The returned segments are arranged by molecule, and only one copy
    of each segment is returned, regardless of how many times it appears
    in this set.
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Segment> > MolGroupsBase::selectAll(const SegID &segid) const
{
    return segid.selectAllFrom(*this);
}

/** Return all of the chains from this set that match the ID 'chainid'.
    The returned chains are arranged by molecule, and only one copy
    of each chain is returned, regardless of how many times it appears
    in this set.
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Chain> > MolGroupsBase::selectAll(const ChainID &chainid) const
{
    return chainid.selectAllFrom(*this);
}

/** Return all of the residues from this set that match the ID 'resid'.
    The returned residues are arranged by molecule, and only one copy
    of each residue is returned, regardless of how many times it appears
    in this set.
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Residue> > MolGroupsBase::selectAll(const ResID &resid) const
{
    return resid.selectAllFrom(*this);
}

/** Return all of the CutGroups from this set that match the ID 'cgid'.
    The returned CutGroups are arranged by molecule, and only one copy
    of each CutGroup is returned, regardless of how many times it appears
    in this set.
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<CutGroup> > MolGroupsBase::selectAll(const CGID &cgid) const
{
    return cgid.selectAllFrom(*this);
}

/** Return all of the atoms from this set that match the ID 'atomid'.
    The returned atoms are arranged by molecule, and only one copy
    of each atom is returned, regardless of how many times it appears
    in this set.
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Atom> > 
MolGroupsBase::selectAll(const AtomID &atomid) const
{
    return atomid.selectAllFrom(*this);
}

/** Return the molecule group that has number 'mgnum'

    \throw SireMol::missing_group
*/
const MolGroup& MolGroupsBase::group(MGNum mgnum) const
{
    return this->at(mgnum);
}

/** Return the molecule group that has name 'mgname'

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
*/
const MolGroup& MolGroupsBase::group(const MGName &mgname) const
{
    return this->at(mgname);
}

/** Return the molecule group at index 'mgidx'

    \throw SireError::invalid_index
*/
const MolGroup& MolGroupsBase::group(MGIdx mgidx) const
{
    return this->at(mgidx);
}

/** Return the molecule group that matches the ID 'mgid'

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
    \throw SireError::invalid_index
*/
const MolGroup& MolGroupsBase::group(const MGID &mgid) const
{
    return this->at(mgid);
}

/** Obvious shortcut for groups(const MGID&)

    \throw SireMol::missing_group
*/
QList<MoleculeGroup> MolGroupsBase::groups(MGNum mgnum) const
{
    return this->selectAll(mgnum);
}

/** Obvious shortcut for groups(const MGID&)

    \throw SireMol::invalid_index
*/
QList<MoleculeGroup> MolGroupsBase::groups(MGIdx mgidx) const
{
    return this->selectAll(mgidx);
}

/** Return all of the groups called 'mgname'

    \throw SireMol::missing_group
*/
QList<MoleculeGroup> MolGroupsBase::groups(const MGName &mgname) const
{
    return this->selectAll(mgname);
}

/** Return all of the groups that match the ID 'mgid'

    \throw SireMol::missing_group
*/
QList<MoleculeGroup> MolGroupsBase::groups(const MGID &mgid) const
{
    return this->selectAll(mgid);
}

/** Return all of the views of the molecule that has number 'molnum'

    \throw SireMol::missing_molecule
*/
ViewsOfMol MolGroupsBase::molecule(MolNum molnum) const
{
    return this->at(molnum);
}

/** Return all of the views of the molecule that matches 'molid'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
ViewsOfMol MolGroupsBase::molecule(const MolID &molid) const
{
    return this->at(molid);
}

/** Obvious shortcut for molecules(const MolID&)

    \throw SireMol::missing_molecule
*/
QList<ViewsOfMol> MolGroupsBase::molecules(MolNum molnum) const
{
    return this->selectAll(molnum);
}

/** Return all of the molecules that match the ID 'molid'

    \throw SireMol::missing_molecule
*/
QList<ViewsOfMol> MolGroupsBase::molecules(const MolID &molid) const
{
    return this->selectAll(molid);
}

/** Return the segment that matches the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
Segment MolGroupsBase::segment(const SegID &segid) const
{
    return this->at(segid);
}

/** Return the chain that matches the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
Chain MolGroupsBase::chain(const ChainID &chainid) const
{
    return this->at(chainid);
}

/** Return the residue that matches the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
Residue MolGroupsBase::residue(const ResID &resid) const
{
    return this->at(resid);
}

/** Return the CutGroup that matches the ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CutGroup MolGroupsBase::cutGroup(const CGID &cgid) const
{
    return this->at(cgid);
}

/** Return the atom that matches the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Atom MolGroupsBase::atom(const AtomID &atomid) const
{
    return this->at(atomid);
}

/** Return all of the segments from this set that match the ID 'segid'.
    The returned segments are arranged by molecule, and only one copy
    of each segment is returned, regardless of how many times it appears
    in this set.
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Segment> > MolGroupsBase::segments(const SegID &segid) const
{
    return this->selectAll(segid);
}

/** Return all of the chains from this set that match the ID 'chainid'.
    The returned chains are arranged by molecule, and only one copy
    of each chain is returned, regardless of how many times it appears
    in this set.
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Chain> > MolGroupsBase::chains(const ChainID &chainid) const
{
    return this->selectAll(chainid);
}

/** Return all of the residues from this set that match the ID 'resid'.
    The returned residues are arranged by molecule, and only one copy
    of each residue is returned, regardless of how many times it appears
    in this set.
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Residue> > MolGroupsBase::residues(const ResID &resid) const
{
    return this->selectAll(resid);
}

/** Return all of the CutGroups from this set that match the ID 'cgid'.
    The returned CutGroups are arranged by molecule, and only one copy
    of each CutGroup is returned, regardless of how many times it appears
    in this set.
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<CutGroup> > MolGroupsBase::cutGroups(const CGID &cgid) const
{
    return this->selectAll(cgid);
}

/** Return all of the atoms from this set that match the ID 'atomid'.
    The returned atoms are arranged by molecule, and only one copy
    of each atom is returned, regardless of how many times it appears
    in this set.
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Atom> > MolGroupsBase::atoms(const AtomID &atomid) const
{
    return this->selectAll(atomid);
}

/** Return whether or not this set contains the group with number 'mgnum' */
bool MolGroupsBase::contains(MGNum mgnum) const
{
    return mgidx_to_num.contains(mgnum);
}

/** Return whether any of the groups contain any view of the molecule
    with number 'molnum' */
bool MolGroupsBase::contains(MolNum molnum) const
{
    return molnum_to_mgnum.contains(molnum);
}

/** Return whether or not any of the groups contains the view 'molview' */
bool MolGroupsBase::contains(const MoleculeView &molview) const
{
    QHash< MolNum,QList<MGNum> >::const_iterator 
                                       it = molnum_to_mgnum.find(molview.number());
                                       
    if (it == molnum_to_mgnum.end())
        return false;
        
    QVarLengthArray<const MolGroup*,10> groups;
    
    this->getGroups(*it, groups);
    
    foreach (const MolGroup *group, groups)
    {
        if (group->contains(molview))
            return true;
    }
    
    return false;
}

/** Return whether or not this set contains all of the views of 
    the molecule in 'molviews'. The views can be contained in 
    multiple groups. */
bool MolGroupsBase::contains(const ViewsOfMol &molviews) const
{
    QHash< MolNum,QList<MGNum> >::const_iterator 
                                       it = molnum_to_mgnum.find(molview.number());
                                       
    if (it == molnum_to_mgnum.end())
        return false;
        
    QVarLengthArray<const MolGroup*,10> groups;
    
    this->getGroups(*it, groups);

    for (int i=0; i<molviews.nViews(); ++i)
    {
        PartialMolecule view = molviews.at(i);
        
        bool found_view = false;
        
        foreach (const MolGroup *group, groups)
        {
            if (group.contains(view))
            {
                found_view = true;
                break;
            }
        }
        
        if (not found_view)
            return false;
    }
    
    return true;
}

/** Return whether or not this set of groups contains all of the views
    of all of the molecules in 'molecules'. These views can be spread
    over lots of groups */
bool MolGroupsBase::contains(const Molecules &molecules) const
{
    for (Molecules::const_iterator it = molecules.begin();
         it != molecules.end();
         ++it)
    {
        if (not this->contains(*it))
            return false;
    }
    
    return true;
}

/** Return whether or not any of the groups in this set contain any
    of the atoms of the view of the molecule in 'molview' */
bool MolGroupsBase::intersects(const MoleculeView &molview) const
{
    QHash< MolNum,QList<MGNum> >::const_iterator 
                                       it = molnum_to_mgnum.find(molview.number());
                                       
    if (it == molnum_to_mgnum.end())
        return false;
        
    QVarLengthArray<const MolGroup*,10> groups;
    
    this->getGroups(*it, groups);
    
    foreach (const MolGroup *group, groups)
    {
        if (group->intersects(molview))
            return true;
    }
    
    return false;
}

/** Return whether any of the groups in this set contain any of the
    atoms of any of the views of any of the molecules in 'molecules' */
bool MolGroupsBase::intersects(const Molecules &molecules) const
{
    for (Molecules::const_iterator it = molecules.begin();
         it != molecules.end();
         ++it)
    {
        if (this->intersects(*it))
            return true;
    }
    
    return false;
}

/** Return the list of molecule groups numbers of groups that 
    contain at least one atom of the molecule with number 'molnum'
    
    \throw SireMol::missing_molecule
*/
const QList<MGNum>& MolGroupsBase::groupsContaining(MolNum molnum) const
{
    QHash< MolNum,QList<MGNum> >::const_iterator it = molnum_to_mgnum.find(molnum);
    
    if (it == molnum_to_mgnum.end())
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule with number %1 is the groups in this set.")
                .arg(molnum), CODELOC );
                
    return *it;
}

/** Return the total number of molecules in the groups in this set */
int MolGroupsBase::nMolecules() const
{
    return molnum_to_mgnum.count();
}

/** Return the total number of views of molecules in the groups in this set.
    Note that if a view appears multiple times, then it will be counted
    multiple times */
int MolGroupsBase::nViews() const
{
    const QHash<MGNum,const MolGroup*> groups = this->getGroups();
    
    int nviews = 0;
    
    for (QHash<MGNum,const MolGroup*>::const_iterator it = groups.constBegin();
         it != groups.constEnd();
         ++it)
    {
        nviews += it->nViews();
    }
    
    return nviews;
}

/** Return the total number of views of the molecule with number  
    'molnum' in the groups in this set. If a view appears multiple
    times then it will be counted multiple times.
    
    \throw SireMol::missing_molecule
*/
int MolGroupsBase::nViews(MolNum molnum) const
{
    const QHash<MGNum,const MolGroup*> groups = 
                                this->getGroups(groupsContaining(molnum));
    
    int nviews = 0;
    
    for (QHash<MGNum,const MolGroup*>::const_iterator it = groups.constBegin();
         it != groups.constEnd();
         ++it)
    {
        nviews += it->nViews(molnum);
    }
    
    return nviews;
}

/** Return whether or not this set is empty (contains no groups) */
bool MolGroupsBase::isEmpty() const
{
    return mgidx_to_num.isEmpty();
}

/** Return the complete set of all molecules in this group. If a view of a 
    molecule appears multiple times in this set then multiple copies of 
    that view will be placed into the returned Molecules object. 
    Note that this is a potentially very slow operation! */
Molecules MolGroupsBase::molecules() const
{
    Molecules all_mols;
    
    const QHash<MGNum,const MolGroup*> groups = this->getGroups();
    
    for (QHash<MGNum,const MolGroup*>::const_iterator = groups.constBegin();
         it != groups.constEnd();
         ++it)
    {
        all_mols += it->molecules();
    }
    
    return all_mols;
}

/** Return the numbers of all molecules that contain at least
    one atom in any of the groups of this set */
QSet<MolNum> MolGroupsBase::molNums() const
{
    return molnum_to_mgnum.keys().toSet();
}

/** Return the numbers of all molecule groups in this set */
QSet<MGNum> MolGroupsBase::mgNums() const
{
    return mgidx_to_num.toSet();
}

/** Assert that this set contains at least one atom of the 
    molecule with number 'molnum'
    
    \throw SireMol::missing_molecule
*/
void MolGroupsBase::assertContains(MolNum molnum) const
{
    if (not molnum_to_mgnum.contains(molnum))
        throw SireMol::missing_molecule( QObject::tr(
            "None of the groups in this set contain the molecule with "
            "number %1.")
                .arg(molnum), CODELOC );
}

/** Assert that this set contains at least one atom of any
    molecule that is identified by the ID 'molid'
    
    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
void MolGroupsBase::assertContains(const MolID &molid) const
{
    this->map(molid);
}

/** Assert that this contains the molecule group with number 'mgnum'

    \throw SireMol::missing_group
*/
void MolGroupsBase::assertContains(MGNum mgnum) const
{
    if (not mgidx_to_num.contains(mgnum))
        throw SireMol::missing_group( QObject::tr(
            "This set does not contain the molecule group with "
            "number %1. Contained groups have numbers %2.")
                .arg(mgnum).arg(Sire::toString(mgidx_to_num)), CODELOC );
}

/** Assert that this contains at least one molecule group that
    is identified by the ID 'mgid'
    
    \throw SireMol:missing_group
    \throw SireError::invalid_index
*/
void MolGroupsBase::assertContains(const MGID &mgid) const
{
    this->map(mgid);
}

/** Synonym for MolGroupsBase::addIfUnique

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroupsBase::unite(const MoleculeView &molview, const MGID &mgid)
{
    this->addIfUnique(molview, mgid);
}

/** Synonym for MolGroupsBase::addIfUnique

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroupsBase::unite(const ViewsOfMol &molviews, const MGID &mgid)
{
    this->addIfUnique(molviews, mgid);
}

/** Synonym for MolGroupsBase::addIfUnique

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroupsBase::unite(const Molecules &molecules, const MGID &mgid)
{
    this->addIfUnique(molecules, mgid);
}

/** Synonym for MolGroupsBase::addIfUnique

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroupsBase::unite(const MolGroup &molgroup, const MGID &mgid)
{
    this->addIfUnique(molgroup, mgid);
}

/** Remove the view of the molecule in 'molview' from all of 
    the groups in this set. This does nothing if this exact
    view is not contained by any of the groups. If a group has
    multiple copies of this view, then this removes only
    the first copy. */
void MolGroupsBase::remove(const MoleculeView &molview)
{
    //get the groups containing this molecule
    QList<MGNum> mgnums = molnum_to_mgnum.value( molview.data().number() );

    if (not mgnums.isEmpty())
        this->remove( molview, MGNumList(mgnums) );
}

/** Remove the views of the molecule in 'molviews' from all of
    the groups in this set. This does nothing if none of these
    views are contained by any of the groups. If a group contains
    multiple copies of a view, then only the first copy is removed. */
void MolGroupsBase::remove(const ViewsOfMol &molviews)
{
    //get the groups containing this molecules
    QList<MGNum> mgnums = molnum_to_mgnum.value( molviews.number() );
    
    if (not mgnums.isEmpty())
        this->remove( molviews, MGNumList(mgnums) );
}

/** Remove all of the views of all of the molecules in 'molecules'
    from all of the groups in this set. If a group contains multiple
    copies of a view then only the first copy is removed */
void MolGroupsBase::remove(const Molecules &molecules)
{
    if (molecules.isEmpty() or this->isEmpty())
        return;

    this->remove( molecules, MGNumList(mgidx_to_num) );
}

/** Remove all of the views of all of the molecules in 'molgroup'
    from this set. Note that this does not remove the actual molecule
    group. If you want to remove the group, then use the 
    MolGroupsBase::remove(MGNum) function. Note that this also
    only removes the first copy of any duplicated views. */
void MolGroupsBase::remove(const MolGroup &molgroup)
{
    this->remove(molgroup.molecules());
}

/** Remove all copies of the view of the molecule in 'molview'
    from all of the groups in this set. This removes all copies
    of a view (even duplicate copies) */
void MolGroupsBase::removeAll(const MoleculeView &molview)
{
    QList<MGNum> mgnums = molnum_to_mgnum.value(molview.data().number());
    
    if (not mgnums.isEmpty())
        this->removeAll( molview, MGNumList(mgnums) );
}

/** Remove all copies of the views of the molecule in 'molviews'
    from all of the groups in this set. This removes all
    copies of the views (even duplicate copies) */
void MolGroupsBase::removeAll(const ViewsOfMol &molviews)
{
    QList<MGNum> mgnus = molnum_to_mgnum.value(molviews.number());
    
    if (not mgnums.isEmpty())
        this->removeAll( molviews, MGNumList(mgnums) );
}

/** Remove all copies of all views of all molecules in 'molecules'
    from all of the groups in this set. This removes all copies
    of the views (even duplicate copies) */
void MolGroupsBase::removeAll(const Molecules &molecules)
{
    if (molecules.isEmpty() or this->isEmpty())
        return;
        
    this->removeAll(molecules, MGNumList(mgidx_to_num));
}

/** Remove all copies of all views of all molecules in the 
    group 'molgroup' from this set. Note that this removes
    the molecules, not the group. Note also that all copies
    of the views are removed (even duplicate copies) */
void MolGroupsBase::removeAll(const MolGroup &molgroup)
{
    this->removeAll(molgroup.molecules());
}

/** Completely remove all views of the molecule with number
    'molnum' from all of the groups from this set. This
    does nothing if there are no views of this molecule
    in any of the groups  */
void MolGroupsBase::remove(MolNum molnum)
{
    QList<MGNum> mgnums = molnum_to_mgnum.value(molnum);
    
    if (not mgnums.isEmpty())
        this->remove(molnum, MGNumList(mgnums));
}

/** Completely remove all views of the molecules whose numbers
    are in 'molnums' from all of the groups in this set. This
    does nothing if there are no views of these molecules in
    any of the groups */
void MolGroupsBase::remove(const QSet<MolNum> &molnums)
{
    if (molnums.isEmpty())
        return;
    else if (molnums.count() == 1)
    {
        this->remove( *(molnums.begin()) );
        return;
    }

    //get the list of groups that contain these molecules
    QList<MGNum> mgnums;
    
    foreach (MolNum molnum, molnums)
    {
        mgnums.append( molnum_to_mgnum.value(molnum) );
    }
    
    mgnums = mgnums.toSet().toList();
    
    if (not mgnums.isEmpty())
        this->remove(molnums, MGNumList(mgnums));
}

/** Completely clear all of the groups in this set */
void MolGroupsBase::removeAll()
{
    if (this->nMolecules() > 1)
    {
        this->removeAll( MGNumList(mgidx_to_num) );
    }
}

/** Update the copies in this set of the molecule viewed in 'molview' 
    to use the same version as 'molview' */
void MolGroupsBase::update(const MoleculeView &molview)
{
    this->update(molview.data());
}

/** Return a reference to the molecule data for the molecule whose data
    is in 'moldata' that is at the same version as the data that is
    already present in the groups of this set. This only returns
    'moldata' if this molecule is not in this set */
const MoleculeData& 
MolGroupsBase::matchToExistingVersion(const MoleculeData &moldata) const
{
    MolNum molnum = moldata.number();

    QHash< MolNum,QList<MGNum> >::const_iterator 
                                    it = molnum_to_mgnum.find(molnum);
                                    
    if (it == molnum_to_mgnum.end())
        //the molecule is not in this set - just return moldata
        return moldata;
        
    //get this molecule...
    const MoleculeData &existing_data = this->getGroup(it->first())
                                                 .molecule(molnum).data();

    return existing_data;
}

/** Return whether or not we need to update the molecule whose
    data is in 'moldata' */
bool MolGroupsBase::needToUpdate(const MoleculeData &moldata) const
{
    return moldata != this->matchToExistingVersion(moldata);
}

/** Return the set of molecules that is a copy of 'molecules', but where
    each molecule has been updated to match the version that exists in 
    the groups in this set */
Molecules MolGroupsBase::matchToExistingVersion(const Molecules &molecules) const
{
    Molecules updated_mols = molecules;
    
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        updated_mols.update( this->matchToExistingVersion(it->data()) );
    }
    
    return updated_mols;
}

/** Add the molecule group and all of its contained molecules to the index */
void MolGroupsBase::addToIndex(const MolGroup &molgroup)
{
    if (mgidx_to_num.contains(molgroup.number()))
        //this group is already in the index!
        return;
    
    //record the number of the group - this allows the groups
    //to be indexed    
    mgidx_to_num.append(molgroup.number());
    
    //record the name of the group - this allows the groups
    //to be accessed by name
    mgname_to_mgnum[molgroup.name()].append(molgroup.number());
    
    //now add the contents of this group to the index
    this->addToIndex(molgroup.number(), molgroup.molNums());
}

/** Add the molecule with number 'molnum' to the index of the group
    with number 'mgnum' */
void MolGroupsBase::addToIndex(MGNum mgnum, MolNum molnum)
{
    BOOST_ASSERT( mgidx_to_num.contains(mgnum) );
    
    QHash< MolNum,QList<MGNum> >::const_iterator 
                                it = molnum_to_mgnum.constFind(molnum);
                                
    if (it != molnum_to_mgnum.constEnd() and 
        it->contains(mgnum))
    {
        //this molecule is already recorded to be in the group
        return;
    }
                                
    molnum_to_mgnum[molnum].append(mgnum);
}

/** Add the set of molecules whose numbers are in 'molnums' to the 
    index of the group whose number is in 'mgnum' */
void MolGroupsBase::addToIndex(MGNum mgnum, const QSet<MolNum> &molnums)
{
    BOOST_ASSERT( mgidx_to_num.contains(mgnum) );
    molnum_to_mgnum[mgnum].insert
}

/** Completely remove the group with number 'mgnum' from the index */
void MolGroupsBase::removeFromIndex(MGNum mgnum)
{
    if (mgidx_to_num.removeAll(mgnum) == 0)
        //no groups were removed
        return;
        
    //ok, the group was removed - remove it from the other indicies
    this->clearIndex(mgnum);

    QMutableHashIterator< MGName,QList<MGNum> > it( mgname_to_mgnum );
    
    while (it.hasNext())
    {
        it.next();
        
        it.value().removeAll(mgnum);
        
        if (it.value().isEmpty())
            it.remove();
    }
}

/** Completely remove the molecule with number 'molnum' from the index */
void MolGroupsBase::removeFromIndex(MolNum molnum)
{
    molnum_to_mgnum.remove(molnum);
}

/** Remove the molecule with number 'molnum' from the index
    of the group with number 'mgnum' */
void MolGroupsBase::removeFromIndex(MGNum mgnum, MolNum molnum)
{
    QHash< MolNum,QList<MGNum> >::iterator it = molnum_to_mgnum.find(molnum);
    
    if (it != molnum_to_mgnum.end())
    {
        it->removeAll(mgnum);
        
        if (it->isEmpty())
            //there are no more groups that contain this molecule
            molnum_to_mgnum.remove(molnum);
    }
}

/** Remove all of the molecules whose numbers are in 'molnums' from the
    group with number 'mgnum' */
void MolGroupsBase::removeFromIndex(MGNum mgnum, const QSet<MolNum> &molnums)
{
    foreach (MolNum molnum, molnums)
    {
        this->removeFromIndex(mgnum, molnum);
    }
}

/** Completely clear the index of the group with number 'mgnum' - this
    removes the link between all molecules in this index and the group */
void MolGroupsBase::clearIndex(MGNum mgnum)
{
    QMutableHashIterator< MolNum,QList<MGNum> > it( molnum_to_mgnum );
    
    while (it.hasNext())
    {
        it.next();
        
        it.value().removeAll(mgnum);
        
        if (it.value().isEmpty())
            it.remove();
    }
}

/////////////
///////////// Implementation of MolGroups
/////////////

RegisterMetaType<MolGroups> r_molgroups;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, 
                                       const MolGroups &molgroups)
{
    writeHeader(ds, r_molgroups, 1);
    
    SharedDataStream sds(ds);
    sds << molgroups.mgroups 
        << static_cast<const MolGroupsBase&>(molgroups);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, MolGroups &molgroups)
{
    VersionID v = readHeader(ds, r_molgroups);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> molgroups.mgroups
            >> static_cast<MolGroupsBase&>(molgroups);
    }
    else
        throw version_error(r_molgroups, "1", CODELOC);
        
    return ds;
}

/** Construct an empty set of groups */
MolGroups::MolGroups()
          : ConcreteProperty<MolGroups,MolGroupsBase>()
{}

/** Construct a set of groups that contains only the single group
    'molgroup' */
MolGroups::MolGroups(const MolGroup &molgroup)
          : ConcreteProperty<MolGroups,MolGroupsBase>()
{
    this->add(molgroup);
}

/** Construct a set of groups that contains the groups in 'molgroups'.
    The groups will be added in the order that they appear in this list. */
MolGroups::MolGroups(const QList<MoleculeGroup> &molgroups)
          : ConcreteProperty<MolGroups,MolGroupsBase>()
{
    foreach (const MoleculeGroup &molgroup, molgroups)
    {
        this->add(molgroup);
    }
}

/** Copy constructor */
MolGroups::MolGroups(const MolGroups &other)
          : ConcreteProperty<MolGroups,MolGroupsBase>(other),
            mgroups(other.mgroups)
{}

/** Destructor */
MolGroups::~MolGroups()
{}

/** Copy assignment operator */
MolGroups& MolGroups::operator=(const MolGroups &other)
{
    if (this != &other)
    {
        mgroups = other.mgroups;
        MolGroupsBase::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool MolGroups::operator==(const MolGroups &other) const
{
    return this == &other or mgroups == other.mgroups;
}

/** Comparison operator */
bool MolGroups::operator!=(const MolGroups &other) const
{
    return this != &other and mgroups != other.mgroups;
}

/** Return the molecule group with number 'mgnum'

    \throw SireMol::missing_group
*/
const MolGroup& MolGroups::at(MGNum mgnum) const
{
    QHash<MGNum,MoleculeGroup>::const_iterator it = mgroups.find(mgnum);
    
    if (it == mgroups.end())
        throw SireMol::missing_group( QObject::tr(
            "There is no molecule group with number %1 in this set. "
            "The available groups are [ %2 ].")
                .arg(mgnum)
                .arg( Sire::toString(this->mgNums()) ), CODELOC );
                
    return it->base();
}

/** Add the molecule group 'molgroup' to this set. This does
    nothing if this group is already in this set. This updates
    the molecules in 'molgroup' so that they are at the 
    same version as any existing copies of the molecules
    in this set. */
void MolGroups::add(const MolGroup &molgroup)
{
    if (this->contains(molgroup.number()))
        return;

    //copy this group
    MoleculeGroup new_group(molgroup);
    
    //update the group so that it has the same version of the 
    //molecules as the other groups in this set
    for (QHash<MGNum,MoleculeGroup>::const_iterator it = mgroups.constBegin();
         it != mgroups.constEnd();
         ++it)
    {
        new_group.update(it->molecules());
    }
    
    //add the group to the index
    this->addToIndex(new_group.base());
    
    mgroups.insert(new_group.number(), new_group);
}

/** Update the group 'molgroup'. If this group is in this set,
    then it updates the group to the same version. Then, regardless
    of whether the group is in this set, it then updates all
    molecules in all of the groups so that they have the same
    version number as 'molgroup'. This does nothing if
    molgroup and none of its molecules are in this set */
void MolGroups::update(const MolGroup &molgroup)
{
    QHash<MGNum,MoleculeGroup>::iterator it = mgroups.find(molgroup.number());
    
    if (it != mgroups.end() and 
        molgroup.majorVersion() != it->majorVersion())
    {
        //this group exists in this set and it is at a different
        //major version - this means that its complement of
        //molecules has changed, so there needs to be a change 
        //in the index...
        
        //get the list of current molecules in this group
        QSet<MolNum> old_molnums = it->molNums();
        
        //now the new set of numbers...
        QSet<MolNum> new_molnums = molgroup.molNums();
        
        this->addToIndex(molgroup.number(), new_molnums - old_molnums);
        this->removeFromIndex(molgroup.number(), old_molnums - new_molnums);
    }
    
    //update all of the groups;
    for (it = mgroups.begin() ; it != mgroups.end(); ++it)
    {
        it->update(molgroup);
    }
}

/** Completely remove the group with number 'mgnum'. This 
    does nothing if there is no such group in this set */
void MolGroups::remove(MGNum mgnum)
{
    if (mgroups.contains(mgnum))
    {
        mgroups.remove(mgnum);
        this->removeFromIndex(mgnum);
    }
}

/** Remove the group 'molgroup' from this set. This does nothing
    if this group doesn't exist in this set */
void MolGroups::remove(const MolGroup &molgroup)
{
    this->remove(molgroup.number());
}

/** Remove the groups that match the ID 'mgid' from this set. This
    does nothing if there are no such groups. */
void MolGroups::remove(const MGID &mgid)
{
    try
    {
        QSet<MGNum> mgnums = this->map(mgid);
    
        foreach (MGNum mgnum, mgnums)
        {
            this->remove(mgnum);
        }
    }
    catch(...)
    {}
}

/** Add the view of the molecule in 'molview' to the groups
    identified by 'mgid'. This adds the view as a duplicate
    if it already exists in the group. The version
    of the molecule added is the version already present
    in this set, if it exists.
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::add(const MoleculeView &molview, const MGID &mgid)
{
    //get the number of the groups identified by 'mgid'
    QList<MGNum> mgnums = this->map(mgid);

    PartialMolecule view(molview);
    view.update( this->matchToExistingData(view.data()) );
    
    foreach (MGNum mgnum, mgnums)
    {
        mgroups.find(mgnum)->add(view);
        this->addToIndex(mgnum, view.number());
    }
}

/** Add the views of the molecule in 'molviews' to the groups
    identified by 'mgid'. This adds the view as a duplicate if
    it already exists in a group.
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::add(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    ViewsOfMol views(molviews);
    views.update( this->matchToExistingVersion(views.data()) );
    
    foreach (MGNum mgnum, mgnums)
    {
        mgroups.find(mgnum)->add(views);
        this->addToIndex(mgnum, views.number());
    }
}

/** Add each of the molecules in 'molecules' to the groups
    identified by 'mgid'. This adds the views as duplicates
    if they exist already in a group. Any molecules that
    already exist in any of the groups in this set are
    updated to the versions that are already present
    in this set.
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::add(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    Molecules mols = this->matchToExistingVersion(molecules);
    QSet<MolNum> molnums = mols.molNums();
    
    foreach (MGNum mgnum, mgnums)
    {
        mgroups.find(mgnum)->add(mols);
        this->addToIndex(mgnum, molnums);
    }
}

/** Add the molecules in the group 'molgroup' to the groups
    identified by 'mgid'. This adds the views as duplicates 
    if they already exist, and adds the views in the same
    order as they appear in 'molgroup'. This is slightly less
    efficient than MolGroups::add(const Molecules&), so use
    that function if you don't care about the order.
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::add(const MolGroup &molgroup, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    //update the group to match the molecule versions
    //already present in this group...
    MoleculeGroup group(molgroup);
    group.update( this->matchToExistingVersions(group.molecules()) );
    
    QSet<MolNum> molnums = group.molNums();
    
    foreach (MGNum mgnum, mgnums)
    {
        mgroups.find(mgnum)->add(group);
        this->addToIndex(mgnum,molnums);
    }
}

/** Add the view of the molecule in 'molview' to the groups
    identified by 'mgid'. This only adds the view to a group
    if it doesn't already exist in the group. The version
    of the molecule already present in this set is used if 
    such a molecule already exists.
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::addIfUnique(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    PartialMolecule view(molview);
    view.update( this->matchToExistingVersion(view.data()) );
    
    foreach (MGNum mgnum, mgnums)
    {
        if ( mgroups.find(mgnum)->addIfUnique(view) )
            this->addToIndex(mgnum,view.number());
    }
}

/** Add the views of the molecule in 'molviews' ot the groups
    identified by 'mgid'. This only adds views to groups that
    don't already exist in that group, and uses the existing
    version of the molecule is it already exists in one 
    of the groups of this set.
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::addIfUnique(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    ViewsOfMol views(molviews);
    views.update( this->matchToExistingVersion(views.data()) );
    
    foreach (MGNum mgnum, mgnums)
    {
        ViewsOfMol added_views = groups.find(mgnum)->addIfUnique(views);
        
        if (not added_views.isEmpty())
            this->addToIndex(mgnum,views.number());
    }
}

/** Add all of the views of the molecules in 'molecules' to the groups
    identified by 'mgid'. This only adds views that don't already 
    exist in the group, and uses the version of the molecules that already
    exists in one of the groups of this set (if one exists)
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::addIfUnique(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    Molecules mols = this->matchExistingVersion(molecules);
    QSet<MolNum> molnums = mols.molNums();
    
    foreach (MGNum mgnum, mgnums)
    {
        QList<ViewsOfMol> added_mols = mgroups.find(mgnum)->addIfUnique(mols);
        
        if (not added_mols.isEmpty())
        {
            QSet<MolNum> added_molnums;
            
            foreach (const ViewsOfMol &mol, added_mols)
            {
                added_molnums.insert(mol.number());
            }
        
            this->addToIndex(mgnum,added_molnums);
        }
    }
}

/** This adds all of the views of the molecules in the group
    'molgroup', in the same order as they exist in this group,
    to all of the groups identified by 'mgid'. This only
    adds views to a group that don't already exist in that
    group and uses the existing version of the molecule if 
    it exists anywhere in this set. 
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::addIfUnique(const MolGroup &molgroup, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    //update the group...
    MolGroup group(molgroup);
    group.update( this->matchExistingVersion(group.molecules()) );
    
    //now add the molecules...
    foreach (MGNum mgnum, mgnums)
    {
        QList<ViewsOfMol> added_mols = mgroups.find(mgnum)->addIfUnique(group);
        
        if (not added_mols.isEmpty())
        {
            QSet<MolNum> added_molnums;
            
            foreach (const ViewsOfMol &mol, added_mols)
            {
                added_molnums.insert(mol.number());
            }
            
            this->addToIndex(mgnum,added_molnums);
        }
    }
}

/** Remove the view of the molecule in 'molview' from the groups
    identified by 'mgid'. This only removes the first copy
    of the view from each group, if multiple copies exist.
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::remove(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    MolNum molnum = molview.data().number();
    
    foreach (MGNum mgnum, mgnums)
    {
        MoleculeGroup &mgroup = *(mgroups.find(mgnum));
        
        if (mgroup.remove(molview))
        {
            if (not mgroup.contains(molnum))
                this->removeFromIndex(mgnum,molnum);
        }
    }
}

/** Remove the views of the molecule in 'molviews' from the groups
    identified by 'mgid'. This only removes the first copy of the
    views from each group if they exist multiple times.
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::remove(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    MolNum molnum = molviews.number();
    
    foreach (MGNum mgnum, mgnums)
    {
        MoleculeGroup &mgroup = *(mgroups.find(mgnum));
        
        ViewsOfMol removed_views = mgroup.remove(molviews);
        
        if (not removed_views.isEmpty())
        {
            if (not mgroup.contains(molnum))
                this->removeFromIndex(mgnum,molnum);
        }
    }
}

/** Remove all of the views of the molecules in 'molecules' from 
    the groups identified by 'mgid'. This removes only the first
    copies of the views in each group.
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::remove(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    foreach (MGNum mgnum, mgnums)
    {
        MoleculeGroup &mgroup = *(mgroups.find(mgnum));
    
        QList<ViewsOfMol> removed_mols = mgroup.remove(molecules);

        QSet<MolNum> removed_molnums;
        
        foreach (const ViewsOfMol &removed_mol, removed_mols)
        {
            if (not mgroup.contains(removed_mol.number()))
                removed_molnums.insert(removed_mol.number());
        }
        
        if (not removed_molnums.isEmpty())
            this->removeFromIndex(mgnum, removed_molnums);
    }
}

/** Remove all of the views of the molecules in the group 'molgroup' from 
    the groups identified by 'mgid'. This removes only the first
    copies of the views in each group.
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::remove(const MolGroup &molgroup, const MGID &mgid)
{
    this->remove(molgroup.molecules(), mgid);
}

/** Remove all copies of the view of the molecule in 'molview' from
    the groups identified by 'mgid'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::removeAll(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    MolNum molnum = molview.data().number();
    
    foreach (MGNum mgnum, mgnums)
    {
        MoleculeGroup &mgroup = *(mgroups.find(mgnum));
        
        if (mgroup.removeAll(molview))
        {
            if (not mgroup.contains(molnum))
                this->removeFromIndex(mgnum, molnum);
        }
    }
}

/** Remove all copies of the views of the molecule in 'molviews' from
    the groups identified by 'mgid'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::removeAll(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    MolNum molnum = molviews.number();
    
    foreach (MGNum mgnum, mgnums)
    {
        MoleculeGroup &mgroup = *(mgroups.find(mgnum));
    
        ViewsOfMol removed_views = mgroup.removeAll(molviews);
        
        if (not removed_views.isEmpty())
        {
            if (not mgroup.contains(molnum))
                this->removeFromIndex(mgnum, molnum);
        }
    }
}

/** Remove all copies of the views of the molecules in 'molecules'
    from the groups identified by 'mgid'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::removeAll(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    foreach (MGNum mgnum, mgnums)
    {
        MoleculeGroup &mgroup = *(mgroups.find(mgnum));
    
        QList<ViewsOfMol> removed_mols = mgroup.removeAll(molecules);
        
        QSet<MolNum> removed_molnums;
        
        foreach (const ViewsOfMol &removed_mol, removed_mols)
        {
            if (not mgroup.contains(removed_mol.number()))
                removed_molnums.insert(removed_mol.number());
        }
        
        this->removeFromIndex(mgnum,removed_molnums);
    }
}

/** Remove all of the molecules from all of the groups identified by
    the ID 'mgid'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::removeAll(const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    foreach (MGNum mgnum, mgnums)
    {
        mgroups.find(mgnum)->removeAll();
        this->removeFromIndex(mgnum);
    }
}

/** Remove all of the views of the molecules in the group 'molgroup'
    from the groups identified by 'mgid'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::removeAll(const MolGroup &molgroup, const MGID &mgid)
{
    this->removeAll(molgroup.molecules(), mgid);
}

/** Remove all views of the molecule with number 'molnum' from the 
    groups identified by 'mgid'. This does nothing to any groups
    that don't contain this molecule
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::remove(MolNum molnum, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    foreach (MGNum mgnum, mgnums)
    {
        MoleculeGroup &mgroup = *(mgroups.find(mgnum));
    
        ViewsOfMol removed_views = mgroup.remove(molnum);
        
        if (not removed_views.isEmpty())
            this->removeFromIndex(molnum);
    }
}

/** Remove the molecules whose numbers are in 'molnums' from the 
    groups identified by 'mgid'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::remove(const QSet<MolNum> &molnums, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    foreach (MGNum mgnum, mgnums)
    {
        MoleculeGroup &mgroup = *(mgroups.find(mgnum));
    
        QList<ViewsOfMol> removed_mols = mgroup.remove(molnums);
        
        QSet<MolNum> removed_nums;
        
        foreach (const ViewsOfMol &removed_mol, removed_mols)
        {
            removed_nums.insert(removed_mol.number());
        }
        
        this->removeFromIndex(mgnum, removed_nums);
    }
}

/** Update all of the groups to use the version of the molecule
    present in 'moldata' */
void MolGroups::update(const MoleculeData &moldata)
{
    //get the current copy of the molecule...
    if (this->needToUpdate(moldata))
    {
        const QList<MGNum> &mgnums = this->groupsContaining(moldata.number());
        
        foreach (MGNum mgnum, mgnums)
        {
            mgroups.find(mgnum)->update(moldata);
        }
    }
}

/** Update all of the groups to use the versions of the molecules
    held in 'molecules' */
void MolGroups::update(const Molecules &molecules)
{
    for (QHash<MGNum,MoleculeGroup>::iterator it = mgroups.begin();
         it != mgroups.end();
         ++it)
    {
        it->update(molecules);
    }
}

/** Update all of the groups to use the versions of the molecules
    help in the group 'molgroup' */
void MolGroups::update(const MolGroup &molgroup)
{
    QHash<MGNum,MoleculeGroup>::const_iterator 
                                  it = mgroups.constFind(molgroup.number());
                                  
    if (it != mgroups.constEnd() and it->majorVersion() == molgroup.majorVersion()
                                 and it->minorVersion() == molgroup.minorVersion())
    {
        //the group is the same as the existing one, and at the same
        //version, so nothing to do :-)
        return;
    }
    
    this->update(molgroup.molecules());
}

/** Set the contents of the groups identified by 'mgid' so that
    they only contain the view in 'molview'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::setContents(const MGID &mgid, const MoleculeView &molview)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    PartialMolecule view(molview);
    view.update( this->matchToExistingVersion(view.data()) );

    foreach (MGNum mgnum, mgnums)
    {
        if (molgroup.find(mgnum)->setContents(view))
        {
            this->clearIndex(mgnum);
            this->addToIndex(mgnum, view.number());
        }
    }
}

/** Set the contents of the groups identified by 'mgid' so that
    they only contain the views of the molecule in 'molviews'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::setContents(const MGID &mgid, const ViewsOfMol &molviews)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    ViewsOfMol views(molviews);
    views.update( this->matchToExistingVersion(views.data()) );
    
    foreach (MGNum mgnum, mgnums)
    {
        if (molgroup.find(mgnum)->setContents(views))
        {
            this->clearIndex(mgnum);
            this->addToIndex(mgnum, views.number());
        }
    }
}

/** Set the contents of the groups identified by 'mgid' so that
    they contain only the views of the molecules contained in 'molecules'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::setContents(const MGID &mgid, const Molecules &molecules)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    Molecules mols = this->matchToExistingVersion(molecules);
    
    QSet<MolNum> molnums = mols.molNums();
    
    foreach (MGNum mgnum, mgnums)
    {
        if (molgroup.find(mgnum)->setContents(mols))
        {
            this->clearIndex(mgnum);
            this->addToIndex(mgnum,molnums);
        }
    }
}

/** Set the contents of the groups identified by 'mgid' so that
    they contain the same views of the same molecules in the 
    same order as in the group 'molgroup'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroups::setContents(const MGID &mgid, const MolGroup &molgroup)
{
    QList<MGNum> mgnums = this->map(mgid);
    
    QHash<MGNum,MoleculeGroup>::const_iterator 
                                  it = mgroups.constFind(molgroup.number());
                                  
    if (it != mgroups.constEnd() and it->majorVersion() == molgroup.majorVersion()
                                 and it->minorVersion() == molgroup.minorVersion())
    {
        //the group is the same as the existing one, and at the same
        //version, so nothing to do :-)
        return;
    }

    //ok we need to set the contents, but not necessarily at the same
    //versions... (we need to update this group to use the versions
    //of any molecules already present in this group)
    Molecules mols = this->matchToExistingVersion(molgroup.molecules());
    
    QSet<MolNum> molnums = mols.molNums();
    
    foreach (MGNum mgnum, mgnums)
    {
        //try to set the contents using the wrong molecule version...
        MolGroup &group = *(molgroups.find(mgnum));
        
        if (group.setContents(molgroup))
        {
            //ok, we need to update the group to ensure that it uses
            //the same version of molecules as are present in the 
            //rest of the set
            group.update(mols);
            
            //now update the index
            this->clearIndex(mgnum);
            this->addToIndex(mgnum, molnums);
        }
    }
}

/** Protected function used to return a modifiable reference to the 
    group with number 'mgnum'
    
    \throw SireMol::missing_group
*/
MolGroup& MolGroups::getGroup(MGNum mgnum)
{
    QHash<MGNum,MoleculeGroup>::iterator it = mgroups.find(mgnum);
    
    if (it == mgroups.end())
        throw SireMol::missing_group( QObject::tr(
            "Cannot find the MoleculeGroup with number %1. Available "
            "groups are %2")
                .arg(mgnum).arg(Sire::toString(mgroups.keys())), CODELOC );

    return it->base();
}

/** Protected function used to return modifiable pointers to the
    groups whose numbers are in 'mgnums'
    
    \throw SireMol::missing_group
*/
void MolGroups::getGroups(const QList<MGNum> &mgnums,
                          QVarLengthArray<MolGroup*,10> &groups)
{
    groups.clear();
    
    QHash<MGNum,MoleculeGroup>::iterator it;
    
    foreach (MGNum mgnum, mgnums)
    {
        it = mgroups.find(mgnum);
        
        if (it == mgroups.end())
            throw SireMol::missing_group( QObject::tr(
                "Cannot find the MoleculeGroup with number %1. Available "
                "groups are %2")
                    .arg(mgnum).arg(Sire::toString(mgroups.keys())), CODELOC );
        
    
        groups.append( &(it->base()) );
    }
}

/** Protected function used to return a hash of modifiable 
    pointers to all of the groups in this set. */
QHash<MGNum,MolGroup*> MolGroups::getGroups()
{
    QHash<MGNum,MolGroup*> groups;
    
    groups.reserve(mgroups.count());
    
    for (QHash<MGNum,MoleculeGroup>::iterator it = mgroups.begin();
         it != mgroups.end();
         ++it)
    {
        groups.insert( it.key(), &(it->base()) );
    }
    
    return groups;
}

/** Protected function used to return a hash of modifiable 
    pointers to all of the groups in this set. */
QHash<MGNum,const MolGroup*> MolGroups::getGroups() const
{
    QHash<MGNum,const MolGroup*> groups;
    
    groups.reserve(mgroups.count());
    
    for (QHash<MGNum,MoleculeGroup>::const_iterator it = mgroups.begin();
         it != mgroups.end();
         ++it)
    {
        groups.insert( it.key(), &(it.base()) );
    }
    
    return groups;
}
