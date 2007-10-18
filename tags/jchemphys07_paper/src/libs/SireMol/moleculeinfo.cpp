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

#include <QSharedData>

#include "qhash_siremol.h"

#include "moleculeinfo.h"
#include "residueinfo.h"

#include "editmoldata.h"

#include "atominfo.h"
#include "atominfogroup.h"

#include "resnum.h"
#include "resid.h"
#include "atomid.h"
#include "cutgroupid.h"

#include "cgatomid.h"
#include "cgnumatomid.h"
#include "resnumatomid.h"
#include "residatomid.h"
#include "idmolatom.h"

#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireStream;

QDataStream& operator<<(QDataStream&, const SireMol::MoleculeInfoPvt&);
QDataStream& operator>>(QDataStream&, SireMol::MoleculeInfoPvt&);

namespace SireMol
{

/** Private class used to implement MoleculeInfo.

    @author Christopher Woods
*/
class MoleculeInfoPvt : public QSharedData
{

friend QDataStream& ::operator<<(QDataStream&, const MoleculeInfoPvt&);
friend QDataStream& ::operator>>(QDataStream&, MoleculeInfoPvt&);

public:
    MoleculeInfoPvt();
    MoleculeInfoPvt(const EditMolData &moldata);

    MoleculeInfoPvt(const MoleculeInfoPvt &other);

    ~MoleculeInfoPvt();

    MoleculeInfoPvt& operator=(const MoleculeInfoPvt &other);

    bool operator==(const MoleculeInfoPvt &other) const;
    bool operator!=(const MoleculeInfoPvt &other) const;

    void regenerateIndex();

    const AtomInfoGroup& _unsafe_atomInfoGroup(CutGroupID cgid) const;
    const ResidueInfo& _unsafe_residue(ResNum resnum) const;
    const ResidueInfo& _unsafe_residue(ResID resid) const;
    const AtomInfo& _unsafe_atom(const CGAtomID &cgatomid) const;

    const CGAtomID& _unsafe_index(AtomID atomid) const;

    void assertResidueExists(ResNum resnum) const;
    void assertResidueExists(ResID resid) const;
    void assertCutGroupExists(CutGroupID cgid) const;

    void assertAtomExists(AtomID atomid) const;
    void assertAtomExists(const AtomIndex &atom) const;
    void assertAtomExists(const CGAtomID &cgatomid) const;
    void assertAtomExists(const ResNumAtomID &resatomid) const;
    void assertAtomExists(const ResIDAtomID &resatomid) const;

    void assertNAtoms(int n) const;

    /** The name of this molecule */
    QString molname;

    /** Vector of residue numbers of the residues in this molecule, in the
        order that they were added to this molecule. */
    QVector<ResNum> resnums;

    /** Metainfo for all of the atoms in this molecule, indexed by
        CGAtomID (use a hash so that we save memory by copying this
        hash directly to the residues for single-cutgroup molecules
        - e.g. most solvents!) */
    QHash<CutGroupID, AtomInfoGroup> atominfos;

    /** Metainfo for all of the residues, indexed by residue number */
    QHash<ResNum, ResidueInfo> resinfos;

    /** QMap mapping the highest index of the atom in each residue
        to the residue number. This is used to quickly find the 'ith'
        atom in the molecule. */
    QMap<AtomID, ResNum> idx2resnum;

    /** The total number of atoms in this molecule */
    quint32 nats;
};

}

Q_DECLARE_METATYPE(SireMol::MoleculeInfoPvt);

using namespace SireMol;

static const RegisterMetaType<MoleculeInfoPvt> r_molpvt;

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const MoleculeInfoPvt &molinfo)
{
    writeHeader(ds, r_molpvt, 1)
          << molinfo.molname << molinfo.resnums << molinfo.atominfos
          << molinfo.resinfos;

    //don't need to save idx2resnum and nats as these can be regenerated

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, MoleculeInfoPvt &molinfo)
{
    VersionID v = readHeader(ds, r_molpvt);

    if (v == 1)
    {
        ds >> molinfo.molname >> molinfo.resnums >> molinfo.atominfos
           >> molinfo.resinfos;

        //regenerate idx2resnums and nats
        molinfo.regenerateIndex();
    }
    else
        throw version_error(v, "1", r_molpvt, CODELOC);

    return ds;
}

////////////
//////////// Implementation of MoleculeInfoPvt
////////////

/** Null constructor */
MoleculeInfoPvt::MoleculeInfoPvt()
                : QSharedData(),
                  molname(QObject::tr("null")), nats(0)
{}

/** Construct from an EditMolData object */
MoleculeInfoPvt::MoleculeInfoPvt(const EditMolData &moldata)
                : QSharedData(),
                  molname(moldata.name()),
                  resnums(moldata.residueNumbers()),
                  nats(moldata.nAtoms())
{
    //first, get the AtomInfos of all of the atoms in the CutGroups
    int ncg = moldata.nCutGroups();

    if (ncg == 0)
        return;

    atominfos.reserve(ncg);

    foreach (CutGroupNum cgnum, moldata.cutGroupNums())
    {
        atominfos.insert( moldata.cutGroupID(cgnum), moldata.atomInfoGroup(cgnum) );
    }

    //now that we have the 'atominfos' hash, we can pass it to the
    //all of the child ResidueInfo objects that we will now create,
    //one for each residue in the molecule
    resinfos.reserve(resnums.count());

    //also, get the highest index of the atom in each residue
    //so that we can locate it quickly
    int maxidx = -1;

    foreach (ResNum resnum, resnums)
    {
        ResidueInfo resinfo(resnum, moldata, atominfos);

        resinfos.insert( resnum, resinfo );

        int nats = resinfo.nAtoms();

        if (nats > 0)
        {
            maxidx += nats;
            idx2resnum.insert( AtomID(maxidx), resnum );
        }
    }

    #warning Need to update MoleculeInfo to use CutGroupNums!
}

/** Copy constructor */
MoleculeInfoPvt::MoleculeInfoPvt(const MoleculeInfoPvt &other)
                : QSharedData(),
                  molname(other.molname), resnums(other.resnums),
                  atominfos(other.atominfos), resinfos(other.resinfos),
                  idx2resnum(other.idx2resnum), nats(other.nats)
{}

/** Destructor */
MoleculeInfoPvt::~MoleculeInfoPvt()
{}

/** Assignment operator */
MoleculeInfoPvt& MoleculeInfoPvt::operator=(const MoleculeInfoPvt &other)
{
    if (&other != this)
    {
        molname = other.molname;
        resnums = other.resnums;
        atominfos = other.atominfos;
        resinfos = other.resinfos;
        idx2resnum = other.idx2resnum;
        nats = other.nats;
    }

    return *this;
}

/** Comparison operator */
bool MoleculeInfoPvt::operator==(const MoleculeInfoPvt &other) const
{
    return (this == &other) or
           (molname == other.molname and resinfos == other.resinfos);
}

/** Comparison operator */
bool MoleculeInfoPvt::operator!=(const MoleculeInfoPvt &other) const
{
    return (this != &other) and
           (molname != other.molname or resinfos != other.resinfos);
}

/** Return the AtomInfoGroup for the CutGroup with ID == cgid - this does
    not check that cgid is valid! */
const AtomInfoGroup& MoleculeInfoPvt::_unsafe_atomInfoGroup(CutGroupID cgid) const
{
    return atominfos.find(cgid).value();
}

/** Return the ResidueInfo for the residue with number 'resnum' - this
    does not check that 'resnum' is valid */
const ResidueInfo& MoleculeInfoPvt::_unsafe_residue(ResNum resnum) const
{
    return resinfos.find(resnum).value();
}

/** Return the ResidueInfo for the residue with index 'resid' - this
    does not check that 'resid' is valid */
const ResidueInfo& MoleculeInfoPvt::_unsafe_residue(ResID resid) const
{
    return resinfos.find( resnums.constData()[resid] ).value();
}

/** Return the AtomInfo of the atom with index 'cgatomid' - this
    does not check that this is a valid index */
const AtomInfo& MoleculeInfoPvt::_unsafe_atom(const CGAtomID &cgatomid) const
{
    return _unsafe_atomInfoGroup(cgatomid.cutGroupID()).constData()[cgatomid.atomID()];
}

/** Convert the index 'atomid' into a CGAtomID index - this does
    not check that the index is valid */
const CGAtomID& MoleculeInfoPvt::_unsafe_index(AtomID atomid) const
{
    //find the nearest index to 'atomid'
    QMap<AtomID,ResNum>::const_iterator it = idx2resnum.lowerBound(atomid);

    BOOST_ASSERT( it != idx2resnum.end() ); //this should be impossible

    //get the ResidueInfo that contains this index
    const ResidueInfo &resinfo = this->_unsafe_residue(it.value());

    //ok - it.key() contains the index of the last atom in this residue
    // The real index == atomid + nAtoms(res) - it.key() - 1;

    atomid += resinfo.nAtoms();
    atomid -= it.key();
    --atomid;

    return resinfo[atomid];
}

/** Internal function used to recreate idx2resnum and nats */
void MoleculeInfoPvt::regenerateIndex()
{
    nats = 0;
    idx2resnum.clear();

    int nres = resnums.count();
    const ResNum *resarray = resnums.constData();

    for (int i=0; i<nres; ++i)
    {
        ResNum resnum = resarray[i];

        nats += _unsafe_residue(resnum).nAtoms();
        idx2resnum.insert( AtomID(nats-1), resnum );
    }
}

/** Assert that the residue 'resnum' exists - else throw an exception

    \throw SireMol::missing_residue
*/
void MoleculeInfoPvt::assertResidueExists(ResNum resnum) const
{
    if (not resinfos.contains(resnum))
        throw SireMol::missing_residue(QObject::tr(
            "Molecule \"%1\" does not contain a residue with number \"%2\"")
                .arg(molname).arg(resnum), CODELOC);
}

/** Assert that the residue with index 'resid' exists - else throw
    an exception

    \throw SireError::invalid_index
*/
void MoleculeInfoPvt::assertResidueExists(ResID resid) const
{
    if (resid >= resnums.count())
        throw SireError::invalid_index(QObject::tr(
            "Molecule \"%1\" has no residue with index \"%2\" (nResidues() == %3)")
                .arg(molname).arg(resid).arg(resnums.count()), CODELOC);
}

/** Assert the CutGroupID is valid

    \throw SireMol::missing_cutgroup
*/
void MoleculeInfoPvt::assertCutGroupExists(CutGroupID cgid) const
{
    if ( not atominfos.contains(cgid) )
        throw SireMol::missing_cutgroup(QObject::tr(
            "There is no CutGroup with ID == %1 in Molecule \"%2\" (nCutGroups() == %3)")
                .arg(cgid).arg(molname).arg(atominfos.count()), CODELOC);
}

/** Assert the Atom with AtomIndex 'atom' exists

    \throw SireMol::missing_residue
    \throw SireMol::missing_atom
*/
void MoleculeInfoPvt::assertAtomExists(const AtomIndex &atom) const
{
    assertResidueExists(atom.resNum());
    _unsafe_residue(atom.resNum()).assertAtomExists(atom.name());
}

/** Assert that the atom with index 'cgatomid' exists

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
void MoleculeInfoPvt::assertAtomExists(const CGAtomID &cgatomid) const
{
    assertCutGroupExists(cgatomid.cutGroupID());
    _unsafe_atomInfoGroup(cgatomid.cutGroupID()).assertAtomExists(cgatomid.atomID());
}

/** Assert that the atom with index 'resatomid' exists

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
void MoleculeInfoPvt::assertAtomExists(const ResNumAtomID &resatomid) const
{
    assertResidueExists(resatomid.resNum());
    _unsafe_residue(resatomid.resNum()).assertAtomExists(resatomid.atomID());
}

/** Assert that the atom with index 'resatomid' exists

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
void MoleculeInfoPvt::assertAtomExists(const ResIDAtomID &resatomid) const
{
    assertResidueExists(resatomid.resID());
    _unsafe_residue(resnums[resatomid.resID()]).assertAtomExists(resatomid.atomID());
}

/** Assert that the atom with index 'atomid' exists

    \throw SireError::invalid_index
*/
void MoleculeInfoPvt::assertAtomExists(AtomID atomid) const
{
    if ( atomid >= nats )
        throw SireError::invalid_index( QObject::tr(
            "Molecule \"%1\" has no atom with index \"%2\" (nAtoms() == %3)")
                .arg(molname).arg(atomid).arg(nats), CODELOC);

}

/** Assert that the number of atoms in the molecule is equal to 'n'

    \throw SireError::incompatible_error
*/
void MoleculeInfoPvt::assertNAtoms(int n) const
{
    if (nats != n)
        throw SireError::incompatible_error( QObject::tr(
                "The number of atoms in the molecule \"%1\", %2, is not "
                "equal to %3.")
                    .arg(molname).arg(nats).arg(n), CODELOC );
}

////////////
//////////// Implementation of MoleculeInfo
////////////

static const RegisterMetaType<SireMol::MoleculeInfo> r_molinfo;

/** Serialise to a binary data stream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const MoleculeInfo &molinfo)
{
    writeHeader(ds, r_molinfo, 1);

    SharedDataStream(ds) << molinfo.d;

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, MoleculeInfo &molinfo)
{
    VersionID v = readHeader(ds, r_molinfo);

    if (v == 1)
    {
        SharedDataStream(ds) >> molinfo.d;
    }
    else
        throw version_error(v, "1", r_molinfo, CODELOC);

    return ds;
}

/** Shared null object */
static QSharedDataPointer<MoleculeInfoPvt> shared_null;

static const QSharedDataPointer<MoleculeInfoPvt>& getSharedNull()
{
    if ( not shared_null )
        shared_null = new MoleculeInfoPvt();
    
    return shared_null;
}

/** Null constructor */
MoleculeInfo::MoleculeInfo() : d( getSharedNull() )
{}

/** Construct from the passed EditMolData */
MoleculeInfo::MoleculeInfo(const EditMolData &moldata)
             : d( new MoleculeInfoPvt(moldata) )
{}

/** Copy constructor - fast as this class is implicitly shared */
MoleculeInfo::MoleculeInfo(const MoleculeInfo &other) : d(other.d)
{}

/** Destructor */
MoleculeInfo::~MoleculeInfo()
{}

/** Assignment operator - fast as this class is implicitly shared */
MoleculeInfo& MoleculeInfo::operator=(const MoleculeInfo &other)
{
    d = other.d;
    return *this;
}

/** Comparison operator */
bool MoleculeInfo::operator==(const MoleculeInfo &other) const
{
    return *d == *(other.d);
}

/** Comparison operator */
bool MoleculeInfo::operator!=(const MoleculeInfo &other) const
{
    return *d != *(other.d);
}

/** Return a copy of the AtomInfoGroup for the CutGroup with ID == cgid

    \throw SireMol::missing_cutgroup
*/
const AtomInfoGroup& MoleculeInfo::operator[](CutGroupID cgid) const
{
    d->assertCutGroupExists(cgid);
    return d->_unsafe_atomInfoGroup(cgid);
}

/** Return a copy of the AtomInfoGroup for the CutGroup with number 'cgnum'

    \throw SireMol::missing_cutgroup
*/
const AtomInfoGroup& MoleculeInfo::operator[](CutGroupNum cgnum) const
{
    #warning CutGroupNum broken
    throw SireError::incomplete_code( QObject::tr(
                          "Need to update MoleculeInfo to work with CutGroupNum"),
                              CODELOC );

    return AtomInfoGroup();
}

/** Return the ResidueInfo for the Residue at index 'resid'

    \throw SireError::invalid_index
*/
const ResidueInfo& MoleculeInfo::operator[](ResID resid) const
{
    return d->_unsafe_residue( residueNumber(resid) );
}

/** Return the ResidueInfo for the Residue with number 'resnum'

    \throw SireMol::missing_residue
*/
const ResidueInfo& MoleculeInfo::operator[](ResNum resnum) const
{
    d->assertResidueExists(resnum);
    return d->_unsafe_residue(resnum);
}

/** Return the CGAtomID index of the atom with index 'atmid'

    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::operator[](AtomID atmid) const
{
    d->assertAtomExists(atmid);
    return d->_unsafe_index(atmid);
}

/** Return the CGAtomID index of the atom 'atom'

    \throw SireMol::missing_residue
    \throw SireMol::missing_atom
*/
const CGAtomID& MoleculeInfo::operator[](const AtomIndex &atom) const
{
    return this->operator[](atom.resNum()).at(atom.name());
}

/** Return the CGAtomID index of the atom with index 'resatomid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::operator[](const ResNumAtomID &resatomid) const
{
    return this->operator[](resatomid.resNum()).at(resatomid.atomID());
}

/** Return the CGAtomID index of the atom with index 'resatomid'

    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::operator[](const ResIDAtomID &resatomid) const
{
    return this->operator[](resatomid.resID()).at(resatomid.atomID());
}

/** Return the CGAtomID index of the atom with index 'cgatomid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::operator[](const CGNumAtomID &cgatomid) const
{
    #warning CutGroupNum broken
    throw SireError::incomplete_code( QObject::tr(
                          "Need to update ResidueInfo to work with CutGroupNum"),
                              CODELOC );

    return CGAtomID();
}

/** Return the CGAtomID index of the atom with index 'cgatomid'
    (this is a convienience function for use with templates -
    it just returns 'cgatomid' after checking that it is valid)

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::operator[](const CGAtomID &cgatomid) const
{
    this->assertAtomExists(cgatomid);
    return cgatomid;
}

/** Return the CGAtomID index of the atom identified by 'atomid'

    \throw SireMol::missing_residue
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::operator[](const IDMolAtom &atomid) const
{
    return atomid.index(*this);
}

/** Return the AtomInfoGroup for the CutGroup with ID == cgid */
const AtomInfoGroup& MoleculeInfo::at(CutGroupID cgid) const
{
    return this->operator[](cgid);
}

/** Return the AtomInfoGroup for the CutGroup with number 'cgnum' */
const AtomInfoGroup& MoleculeInfo::at(CutGroupNum cgnum) const
{
    return this->operator[](cgnum);
}

/** Return the ResidueInfo for the Residue at index 'resid'

    \throw SireError::invalid_index
*/
const ResidueInfo& MoleculeInfo::at(ResID resid) const
{
    return this->operator[](resid);
}

/** Return the ResidueInfo for the Residue with number 'resnum'

    \throw SireMol::missing_residue
*/
const ResidueInfo& MoleculeInfo::at(ResNum resnum) const
{
    return this->operator[](resnum);
}

/** Return the CGAtomID index of the atom with index 'atmid'

    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::at(AtomID atmid) const
{
    return this->operator[](atmid);
}

/** Return the CGAtomID index of the atom 'atom'

    \throw SireMol::missing_residue
    \throw SireMol::missing_atom
*/
const CGAtomID& MoleculeInfo::at(const AtomIndex &atom) const
{
    return this->operator[](atom);
}

/** Return the CGAtomID index of the atom with index 'resatomid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::at(const ResNumAtomID &resatomid) const
{
    return this->operator[](resatomid);
}

/** Return the CGAtomID index of the atom with index 'resatomid'

    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::at(const ResIDAtomID &resatomid) const
{
    return this->operator[](resatomid);
}

/** Return the CGAtomID index of the atom with index 'cgatomid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::at(const CGNumAtomID &cgatomid) const
{
    return this->operator[](cgatomid);
}

/** Return the CGAtomID index of the atom with index 'cgatomid'
    (this is a convienience function for use with templates -
    it just returns 'cgatomid' after checking that it is valid)

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::at(const CGAtomID &cgatomid) const
{
    return this->operator[](cgatomid);
}

/** Return the CGAtomID index of the atom identified by 'atomid'

    \throw SireMol::missing_residue
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
const CGAtomID& MoleculeInfo::at(const IDMolAtom &atomid) const
{
    return atomid.index(*this);
}

/** Return the AtomInfo for the atom with index 'cgatomid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
const AtomInfo& MoleculeInfo::atom(const CGAtomID &cgatomid) const
{
    this->assertAtomExists(cgatomid);
    return d->_unsafe_atom(cgatomid);
}

/** Return the AtomInfo for the atom with index 'atomid'

    \throw SireMol::missing_residue
    \throw SireMol::missing_cutgroup
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
const AtomInfo& MoleculeInfo::atom(const IDMolAtom &atomid) const
{
    return d->_unsafe_atom( atomid.index(*this) );
}

/** Return the AtomInfoGroup for the CutGroup with ID == cgid

    \throw SireMol::missing_cutgroup
*/
const AtomInfoGroup& MoleculeInfo::atomGroup(CutGroupID cgid) const
{
    return this->at(cgid);
}

/** Return the AtomInfoGroup for the CutGroup with number 'cgnum'

    \throw SireMol::missing_cutgroup
*/
const AtomInfoGroup& MoleculeInfo::atomGroup(CutGroupNum cgnum) const
{
    return this->at(cgnum);
}

/** Return the AtomInfoGroups for all CutGroups, in a hash indexed by
    CutGroupID */
QHash<CutGroupID,AtomInfoGroup> MoleculeInfo::atomGroups() const
{
    return d->atominfos;
}

/** Return the AtomInfoGroups for all of the CutGroups that contain atoms
    in the residue with number 'resnum', in a hash indexed by CutGroupID

    \throw SireMol::missing_residue
*/
QHash<CutGroupID,AtomInfoGroup> MoleculeInfo::atomGroups(ResNum resnum) const
{
    return this->at(resnum).atomGroups();
}

/** Return the AtomInfoGroups for all of the CutGroups that contain atoms
    in the residue with index 'resid', in a hash indexed by CutGroupID

    \throw SireError::invalid_index
*/
QHash<CutGroupID,AtomInfoGroup> MoleculeInfo::atomGroups(ResID resid) const
{
    return this->at(resid).atomGroups();
}

/** Return the AtomInfoGroups for the CutGroups whose IDs are in 'cgids', returning
    a hash indexed by CutGroupID

    \throw SireMol::missing_cutgroup
*/
QHash<CutGroupID,AtomInfoGroup> MoleculeInfo::atomGroups(
                                          const QSet<CutGroupID> &cgids) const
{
    QHash<CutGroupID,AtomInfoGroup> groups;
    groups.reserve(cgids.count());

    for (QSet<CutGroupID>::const_iterator it = cgids.begin();
         it != cgids.end();
         ++it)
    {
        groups.insert( *it, this->at(*it) );
    }

    return groups;
}

/** Return the AtomInfoGroups for all of the CutGroups that contain atoms
    in the residues with numbers in 'resnums', in a hash indexed by CutGroupID

    \throw SireMol::missing_residue
*/
QHash<CutGroupID,AtomInfoGroup> MoleculeInfo::atomGroups(
                                      const QSet<ResNum> &resnums) const
{
    QHash<CutGroupID,AtomInfoGroup> groups;

    for (QSet<ResNum>::const_iterator it = resnums.begin();
         it != resnums.end();
         ++it)
    {
        groups.unite( this->atomGroups(*it) );
    }

    return groups;
}

/** Return the AtomInfoGroups for all of the CutGroups that contain atoms
    in the residues with indexes in 'resids', in a hash indexed by CutGroupID

    \throw SireMol::missing_residue
*/
QHash<CutGroupID,AtomInfoGroup> MoleculeInfo::atomGroups(
                                      const QSet<ResID> &resids) const
{
    QHash<CutGroupID,AtomInfoGroup> groups;

    for (QSet<ResID>::const_iterator it = resids.begin();
         it != resids.end();
         ++it)
    {
        groups.unite( this->atomGroups(*it) );
    }

    return groups;
}

/** Return an array that contains a copy of all of the AtomInfos of the atoms
    in the molecule. The atoms are returned in the order as they are in the
    molecule. */
QVector<AtomInfo> MoleculeInfo::atoms() const
{
    QVector<AtomInfo> atms;
    atms.reserve( this->nAtoms() );

    int nres = this->nResidues();
    const ResNum *resnums = d->resnums.constData();

    for (int i=0; i<nres; ++i)
    {
        atms += this->at(resnums[i]).atoms();
    }

    return atms;
}

/** Return an array that contains a copy of all of the AtomInfos of the
    atoms in the residue with number 'resnum', in the same order as they
    were added to the residue.

    \throw SireMol::missing_residue
*/
QVector<AtomInfo> MoleculeInfo::atoms(ResNum resnum) const
{
    return this->at(resnum).atoms();
}

/** Return an array that contains a copy of all of the AtomInfos of the
    atoms in the residue with index 'resid', in the same order as they
    were added to the residue.

    \throw SireError::invalid_index
*/
QVector<AtomInfo> MoleculeInfo::atoms(ResID resid) const
{
    return this->at(resid).atoms();
}

/** Return an array that contains a copy of all of the AtomInfos of the
    atoms in the CutGroup with ID == cgid, in the same order as they
    exist in the CutGroup

    \throw SireMol::missing_cutgroup
*/
QVector<AtomInfo> MoleculeInfo::atoms(CutGroupID cgid) const
{
    return this->at(cgid).atoms();
}

/** Return an array that contains a copy of all of the AtomInfos of the
    atoms in the CutGroup with number 'cgnum', in the same order as they
    exist in the CutGroup

    \throw SireMol::missing_cutgroup
*/
QVector<AtomInfo> MoleculeInfo::atoms(CutGroupNum cgnum) const
{
    return this->at(cgnum).atoms();
}

template<class T>
QHash<T,AtomInfo> getAtoms(const MoleculeInfo &mol, const QSet<T> &ids)
{
    QHash<T,AtomInfo> atms;
    atms.reserve(ids.count());

    for (typename QSet<T>::const_iterator it = ids.begin();
         it != ids.end();
         ++it)
    {
        atms.insert( *it, mol.atom(*it) );
    }

    return atms;
}

/** Return the AtomInfos for atoms whose IDs are in 'idxs', in a hash indexed
    by AtomID

    \throw SireError::invalid_index
*/
QHash<AtomID,AtomInfo> MoleculeInfo::atoms( const QSet<AtomID> &idxs ) const
{
    return getAtoms<AtomID>(*this, idxs);
}

/** Return the AtomInfos for atoms whose IDs are in 'cgids', in a hash indexed
    by CGAtomID

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QHash<CGAtomID,AtomInfo> MoleculeInfo::atoms( const QSet<CGAtomID> &cgids ) const
{
    return getAtoms<CGAtomID>(*this, cgids);
}

/** Return the AtomInfos for atoms whose numbers are in 'cgnums', in a hash indexed
    by CGNumAtomID

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QHash<CGNumAtomID,AtomInfo> MoleculeInfo::atoms(
                                          const QSet<CGNumAtomID> &cgnums ) const
{
    return getAtoms<CGNumAtomID>(*this, cgnums);
}

/** Return the AtomInfos for atoms whose IDs are in 'rsids', in a hash indexed
    by ResNumAtomID

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QHash<ResNumAtomID,AtomInfo> MoleculeInfo::atoms( const QSet<ResNumAtomID> &rsids) const
{
    return getAtoms<ResNumAtomID>(*this, rsids);
}

/** Return the AtomInfos for atoms whose IDs are in 'rsids', in a hash indexed
    by ResIDAtomID

    \throw SireError::invalid_index
*/
QHash<ResIDAtomID,AtomInfo> MoleculeInfo::atoms( const QSet<ResIDAtomID> &rsids) const
{
    return getAtoms<ResIDAtomID>(*this, rsids);
}

/** Return the AtomInfos for atoms whose AtomIndexes are in 'atms', in a hash indexed
    by AtomIndex

    \throw SireMol::missing_residue
    \throw SireMol::missing_atom
*/
QHash<AtomIndex,AtomInfo> MoleculeInfo::atoms( const QSet<AtomIndex> &atms ) const
{
    return getAtoms<AtomIndex>(*this, atms);
}

template<class T>
QHash< T, QVector<AtomInfo> > getAtomVectors(const MoleculeInfo &molinfo,
                                             const QSet<T> &idxs)
{
    QHash< T, QVector<AtomInfo> > allatms;
    allatms.reserve(idxs.count());

    for (typename QSet<T>::const_iterator it = idxs.begin();
         it != idxs.end();
         ++it)
    {
        allatms.insert( *it, molinfo.atoms(*it) );
    }

    return allatms;
}

/** Return the arrays of AtomInfos of the atoms in the CutGroups whose IDs are in
    'cgids', in a hash indexed by CutGroupID

    \throw SireMol::missing_cutgroup
*/
QHash< CutGroupID, QVector<AtomInfo> > MoleculeInfo::atoms(const QSet<CutGroupID> &cgids) const
{
    return getAtomVectors<CutGroupID>(*this, cgids);
}

/** Return the arrays of AtomInfos of the atoms in the CutGroups whose numbers are in
    'cgnums', in a hash indexed by CutGroupNum

    \throw SireMol::missing_cutgroup
*/
QHash< CutGroupNum, QVector<AtomInfo> > MoleculeInfo::atoms(
                                                const QSet<CutGroupNum> &cgnums) const
{
    return getAtomVectors<CutGroupNum>(*this, cgnums);
}

/** Return the arrays of AtomInfos of the atoms in the residues whose numbers
    are in 'resnums', in a hash indexed by ResNum

    \throw SireMol::missing_residue
*/
QHash< ResNum, QVector<AtomInfo> > MoleculeInfo::atoms(const QSet<ResNum> &resnums) const
{
    return getAtomVectors<ResNum>(*this, resnums);
}

/** Return the arrays of AtomInfos of the atoms in the residues whose indicies
    are in 'resids', in a hash indexed by ResID

    \throw SireError::invalid_index
*/
QHash< ResID, QVector<AtomInfo> > MoleculeInfo::atoms(const QSet<ResID> &resids) const
{
    return getAtomVectors<ResID>(*this, resids);
}

/** Return the ResidueInfo for the residue that contains the atom
    with index 'atomid'

    \throw SireError::invalid_index
*/
const ResidueInfo& MoleculeInfo::residue(AtomID atmid) const
{
    d->assertAtomExists(atmid);

    //find the nearest index to 'atomid'
    QMap<AtomID,ResNum>::const_iterator it = d->idx2resnum.lowerBound(atmid);

    BOOST_ASSERT( it != d->idx2resnum.end() ); //this should be impossible

    //return the ResidueInfo that contains this index
    return d->_unsafe_residue(it.value());
}

/** Return the ResidueInfo for the residue with index 'resid'

    \throw SireError::invalid_index
*/
const ResidueInfo& MoleculeInfo::residue(ResID resid) const
{
    return this->at(resid);
}

/** Return the ResidueInfo for the residue with number 'resnum'

    \throw SireMol::missing_residue
*/
const ResidueInfo& MoleculeInfo::residue(ResNum resnum) const
{
    return this->at(resnum);
}

/** Return a string representation of this object */
QString MoleculeInfo::toString() const
{
    return QObject::tr("Molecule(\"%1\"): nAtoms() == %2, nResidues() == %3, "
                       "nCutGroups() == %4")
                    .arg(name()).arg(nAtoms()).arg(nResidues()).arg(nCutGroups());
}

/** Return the name of the molecule */
QString MoleculeInfo::name() const
{
    return d->molname;
}

/** Return the name of the residue with number 'resnum'

    \throw SireMol::missing_residue
*/
QString MoleculeInfo::residueName(ResNum resnum) const
{
    return this->at(resnum).name();
}

/** Return the name of the residue with index 'resid'

    \throw SireError::invalid_index
*/
QString MoleculeInfo::residueName(ResID resid) const
{
    return this->at(resid).name();
}

/** Return the number of the residue with index 'resid'

    \throw SireError::invalid_index
*/
ResNum MoleculeInfo::residueNumber(ResID resid) const
{
    d->assertResidueExists(resid);
    return d->resnums.constData()[resid];
}

/** Return the number of the first residue in the molecule that is
    called 'resname'

    \throw SireMol::missing_residue
*/
ResNum MoleculeInfo::residueNumber(const QString &resname) const
{
    for (QVector<ResNum>::const_iterator it = d->resnums.begin();
         it != d->resnums.end();
         ++it)
    {
        if (d->_unsafe_residue(*it).name() == resname)
            return *it;
    }

    throw SireMol::missing_residue( QObject::tr(
                  "There is no residue called \"%1\" in the molecule called \"%2\"")
                        .arg(resname).arg(this->name()), CODELOC );

    return ResNum(0);
}

/** Return the CutGroupID of the CutGroup with number 'cgnum'

    \throw SireMol::missing_cutgroup
*/
CutGroupID MoleculeInfo::cutGroupID(CutGroupNum cgnum) const
{
    #warning CutGroupNum broken
    throw SireError::incomplete_code( QObject::tr(
                          "Need to update ResidueInfo to work with CutGroupNum"),
                              CODELOC );

    return CutGroupID(0);
}

/** Return the CutGroupNum of the CutGroup with ID == cgid

    \throw SireMol::missing_cutgroup
*/
CutGroupNum MoleculeInfo::cutGroupNum(CutGroupID cgid) const
{
    #warning CutGroupNum broken
    throw SireError::incomplete_code( QObject::tr(
                          "Need to update ResidueInfo to work with CutGroupNum"),
                              CODELOC );

    return CutGroupNum(0);
}

/** Return the number of residues in this molecule */
int MoleculeInfo::nResidues() const
{
    return d->resnums.count();
}

/** Return the number of CutGroups in the molecule */
int MoleculeInfo::nCutGroups() const
{
    return d->atominfos.count();
}

/** Return the number of atoms in the molecule */
int MoleculeInfo::nAtoms() const
{
    return d->nats;
}

/** Return the number of atoms in the residue with number 'resnum'

    \throw SireMol::missing_residue
*/
int MoleculeInfo::nAtoms(ResNum resnum) const
{
    return this->at(resnum).nAtoms();
}

/** Return the number of atoms in the residue at index 'resid'

    \throw SireError::invalid_index
*/
int MoleculeInfo::nAtoms(ResID resid) const
{
    return this->at(resid).nAtoms();
}

/** Return the number of atoms in the CutGroup with ID == 'cgid'

    \throw SireMol::missing_cutgroup
*/
int MoleculeInfo::nAtoms(CutGroupID cgid) const
{
    return this->at(cgid).nAtoms();
}

/** Return the number of atoms in the CutGroup with number 'cgnum'

    \throw SireMol::missing_cutgroup
*/
int MoleculeInfo::nAtoms(CutGroupNum cgnum) const
{
    return this->nAtoms( this->cutGroupID(cgnum) );
}

/** Return the numbers of the CutGroups in this molecule, in
    CutGroupID order */
QVector<CutGroupNum> MoleculeInfo::cutGroupNums() const
{
    #warning CutGroupNum broken
    throw SireError::incomplete_code( QObject::tr(
                          "Need to update ResidueInfo to work with CutGroupNum"),
                              CODELOC );

    return QVector<CutGroupNum>();
}

/** Return the array of the numbers of the residues, in the order that
    they appear in the molecule */
QVector<ResNum> MoleculeInfo::residueNumbers() const
{
    return d->resnums;
}

/** Return the numbers of all residues that are called 'resname'. This
    returns an empty array if there are no residues with this name */
QVector<ResNum> MoleculeInfo::residueNumbers(const QString &resname) const
{
    QVector<ResNum> resnums;

    for (QVector<ResNum>::const_iterator it = d->resnums.begin();
         it != d->resnums.end();
         ++it)
    {
        if (d->_unsafe_residue(*it).name() == resname)
            resnums.append(*it);
    }

    return resnums;
}

/** Return the residue numbers of all residues that have atoms in the
    CutGroup with ID == 'cgid'

    \throw SireMol::missing_cutgroup
*/
QVector<ResNum> MoleculeInfo::residueNumbers(CutGroupID cgid) const
{
    return this->at(cgid).residueNumbers();
}

/** Return the residue numbers of all residues that have atoms in the
    CutGroup with number 'cgnum'

    \throw SireMol::missing_cutgroup
*/
QVector<ResNum> MoleculeInfo::residueNumbers(CutGroupNum cgnum) const
{
    return this->residueNumbers( this->cutGroupID(cgnum) );
}

/** Return the list of the names of all of the residues in this molecule,
    in the order that they exist in this molecule */
QStringList MoleculeInfo::residueNames() const
{
    QStringList resnames;

    for (QVector<ResNum>::const_iterator it = d->resnums.begin();
         it != d->resnums.end();
         ++it)
    {
        resnames.append( d->_unsafe_residue(*it).name() );
    }

    return resnames;
}

/** Return the names of all of the atoms in the residue with number 'resnum'

    \throw SireMol::missing_residue
*/
QStringList MoleculeInfo::atomNames(ResNum resnum) const
{
    return this->at(resnum).atomNames();
}

/** Return the names of all of the atoms in the residue with index 'resid'

    \throw SireError::invalid_index
*/
QStringList MoleculeInfo::atomNames(ResID resid) const
{
    return this->at(resid).atomNames();
}

/** Return whether or not this molecule contains a CutGroup with ID == 'cgid' */
bool MoleculeInfo::contains(CutGroupID cgid) const
{
    return cgid < this->nCutGroups();
}

/** Return whether or not this molecule contains a CutGroup with number 'cgnum' */
bool MoleculeInfo::contains(CutGroupNum cgnum) const
{
    #warning CutGroupNum broken
    throw SireError::incomplete_code( QObject::tr(
                          "Need to update ResidueInfo to work with CutGroupNum"),
                              CODELOC );

    return false;
}

/** Return whether or not this molecule contains a residue with number 'resnum' */
bool MoleculeInfo::contains(ResNum resnum) const
{
    return d->resinfos.contains(resnum);
}

/** Return whether or not this molecule contains a residue with index 'resid' */
bool MoleculeInfo::contains(ResID resid) const
{
    return resid < this->nResidues();
}

/** Return whether or not the molecule contains an atom called 'atomname'
    in a residue with number 'resnum' */
bool MoleculeInfo::contains(ResNum resnum, const QString &atomname) const
{
    return this->contains(resnum) and d->_unsafe_residue(resnum).contains(atomname);
}

/** Return whether or not the molecule contains an atom with AtomIndex 'atm' */
bool MoleculeInfo::contains(const AtomIndex &atm) const
{
    return this->contains(atm.resNum(), atm.name());
}

/** Return whether or not the molecule contains an atom with index 'atomid' in the
    CutGroup with ID == cgid */
bool MoleculeInfo::contains(CutGroupID cgid, AtomID atomid) const
{
    return this->contains(cgid) and d->_unsafe_atomInfoGroup(cgid).contains(atomid);
}

/** Return whether or not the molecule contains an atom with index 'cgatomid' */
bool MoleculeInfo::contains(const CGAtomID &cgatomid) const
{
    return this->contains(cgatomid.cutGroupID(), cgatomid.atomID());
}

/** Return whether or not the molecule contains an atom with index 'atomid' in the
    CutGroup with number 'cgnum' */
bool MoleculeInfo::contains(CutGroupNum cgnum, AtomID atomid) const
{
    #warning CutGroupNum broken
    throw SireError::incomplete_code( QObject::tr(
                          "Need to update ResidueInfo to work with CutGroupNum"),
                              CODELOC );

    return false;
}

/** Return whether or not the molecule contains an atom with index 'cgatomid' */
bool MoleculeInfo::contains(const CGNumAtomID &cgatomid) const
{
    return this->contains(cgatomid.cutGroupNum(), cgatomid.atomID());
}

/** Return whether or not the molecule contains an atom with index 'atomid'
    in a residue with number 'resnum' */
bool MoleculeInfo::contains(ResNum resnum, AtomID atomid) const
{
    return this->contains(resnum) and d->_unsafe_residue(resnum).contains(atomid);
}

/** Return whether or not the molecule contains an atom with index 'resatomid' */
bool MoleculeInfo::contains(const ResNumAtomID &resatomid) const
{
    return this->contains( resatomid.resNum(), resatomid.atomID() );
}

/** Return whether or not the molecule contains an atom with index 'atomid'
    in a residue with index 'resid' */
bool MoleculeInfo::contains(ResID resid, AtomID atomid) const
{
    return this->contains(resid) and d->_unsafe_residue(resid).contains(atomid);
}

/** Return whether or not the molecule contains an atom with index 'resatomid' */
bool MoleculeInfo::contains(const ResIDAtomID &resatomid) const
{
    return this->contains(resatomid.resID(), resatomid.atomID());
}

/** Return whether or not this molecule contains the atom at 'atomid' */
bool MoleculeInfo::contains(const IDMolAtom &atomid) const
{
    return atomid.contains(*this);
}

/** Return whether or not this molecule is empty (contains no atoms) */
bool MoleculeInfo::isEmpty() const
{
    return this->nAtoms() == 0;
}

/** Return whether the residue with number 'resnum' is empty

    \throw SireMol::missing_residue
*/
bool MoleculeInfo::isEmpty(ResNum resnum) const
{
    return this->at(resnum).isEmpty();
}

/** Return whether the residue at index 'resid' is empty

    \throw SireError::invalid_index
*/
bool MoleculeInfo::isEmpty(ResID resid) const
{
    return this->at(resid).isEmpty();
}

/** Return whether the CutGroup with ID == cgid is empty

    \throw SireMol::missing_cutgroup
*/
bool MoleculeInfo::isEmpty(CutGroupID cgid) const
{
    return this->at(cgid).isEmpty();
}

/** Return whether the CutGroup with number 'cgnum' is empty

    \throw SireMol::missing_cutgroup
*/
bool MoleculeInfo::isEmpty(CutGroupNum cgnum) const
{
    return this->isEmpty( this->cutGroupID(cgnum) );
}

/** Assert that the residue 'resnum' exists - else throw an exception

    \throw SireMol::missing_residue
*/
void MoleculeInfo::assertResidueExists(ResNum resnum) const
{
    d->assertResidueExists(resnum);
}

/** Assert that the residue with index 'resid' exists - else throw
    an exception

    \throw SireError::invalid_index
*/
void MoleculeInfo::assertResidueExists(ResID resid) const
{
    d->assertResidueExists(resid);
}

/** Assert the CutGroupID is valid

    \throw SireMol::missing_cutgroup
*/
void MoleculeInfo::assertCutGroupExists(CutGroupID cgid) const
{
    d->assertCutGroupExists(cgid);
}

/** Assert the CutGroupNum is valid

    \throw SireMol::missing_cutgroup
*/
void MoleculeInfo::assertCutGroupExists(CutGroupNum cgnum) const
{
    #warning CutGroupNum broken
    throw SireError::incomplete_code( QObject::tr(
                          "Need to update ResidueInfo to work with CutGroupNum"),
                              CODELOC );
}

/** Assert that the AtomIndex is valid

    \throw SireMol::missing_residue
    \throw SireMol::missing_atom
*/
void MoleculeInfo::assertAtomExists(const AtomIndex &atom) const
{
    d->assertAtomExists(atom);
}

/** Assert that the index is valid

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
void MoleculeInfo::assertAtomExists(const CGAtomID &cgatomid) const
{
    d->assertAtomExists(cgatomid);
}

/** Assert that the index is valid

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
void MoleculeInfo::assertAtomExists(const CGNumAtomID &cgatomid) const
{
    #warning CutGroupNum broken
    throw SireError::incomplete_code( QObject::tr(
                          "Need to update ResidueInfo to work with CutGroupNum"),
                              CODELOC );
}

/** Assert that the index is valid

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
void MoleculeInfo::assertAtomExists(const ResNumAtomID &resatomid) const
{
    d->assertAtomExists(resatomid);
}

/** Assert that the index is valid

    \throw SireError::invalid_index
*/
void MoleculeInfo::assertAtomExists(const ResIDAtomID &resatomid) const
{
    d->assertAtomExists(resatomid);
}

/** Assert that the index is valid

    \throw SireError::invalid_index
*/
void MoleculeInfo::assertAtomExists(AtomID atomid) const
{
    d->assertAtomExists(atomid);
}

/** Assert that the atom at index 'atomid' exists

    \throw SireMol::missing_residue
    \throw SireMol::missing_cutgroup
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void MoleculeInfo::assertAtomExists(const IDMolAtom &atomid) const
{
    atomid.index(*this);
}

/** Assert that the number of atoms in the molecule is equal
    to 'nats'

    \throw SireError::incompatible_error
*/
void MoleculeInfo::assertNAtoms(int nats) const
{
    d->assertNAtoms( nats );
}
