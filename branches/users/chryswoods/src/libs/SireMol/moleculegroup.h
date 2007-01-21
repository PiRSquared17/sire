#ifndef SIREMOL_MOLECULEGROUP_H
#define SIREMOL_MOLECULEGROUP_H

#include <QVector>
#include <QHash>
#include <QSharedDataPointer>

#include "moleculegroupid.h"

#include "molecule.h"
#include "moleculeid.h"

#include "SireBase/idmajminversion.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

namespace detail
{
class MoleculeGroupPvt;
}

class MoleculeGroup;
}

QDataStream& operator<<(QDataStream&, const SireMol::detail::MoleculeGroupPvt&);
QDataStream& operator>>(QDataStream&, SireMol::detail::MoleculeGroupPvt&);

QDataStream& operator<<(QDataStream&, const SireMol::MoleculeGroup&);
QDataStream& operator>>(QDataStream&, SireMol::MoleculeGroup&);

namespace SireMol
{

using SireBase::Version;

namespace detail
{

/** This is the private implementation of MoleculeGroup - this
    allows the MoleculeGroup class to be implicitly shared.

    @author Christopher Woods
*/
class SIREMOL_EXPORT MoleculeGroupPvt : public QSharedData
{

friend QDataStream& ::operator<<(QDataStream&, const MoleculeGroupPvt&);
friend QDataStream& ::operator>>(QDataStream&, MoleculeGroupPvt&);

public:
    MoleculeGroupPvt();
    MoleculeGroupPvt(const QString &name);

    MoleculeGroupPvt(const QString &name,
                     const QVector<Molecule> &molecules);

    MoleculeGroupPvt(const MoleculeGroupPvt &other);

    ~MoleculeGroupPvt();

    MoleculeGroupPvt& operator=(const MoleculeGroupPvt &other);

    bool operator==(const MoleculeGroupPvt &other) const;
    bool operator!=(const MoleculeGroupPvt &other) const;

    const QVector<Molecule>& molecules() const;
    const QHash<MoleculeID,int>& index() const;

    QSet<MoleculeID> moleculeIDs() const;

    const Molecule& molecule(MoleculeID molid) const;

    const QString& name() const;

    MoleculeGroupID ID() const;
    const Version& version() const;

    bool add(const Molecule &molecule);
    bool change(const Molecule &molecule);
    bool remove(const Molecule &molecule);

    bool add(const QVector<Molecule> &molecules);
    bool change(const QVector<Molecule> &molecules);
    bool remove(const QVector<Molecule> &molecules);

    bool change(const QHash<MoleculeID,Molecule> &molecules);

    void rename(const QString &newname);

    bool contains(MoleculeID molid) const;

    void assertContains(MoleculeID molid) const;

private:
    /** Incremint used to increment that ID number */
    static SireBase::Incremint molgroup_incremint;

    void reindex();

    /** The array of all of the molecules in this group */
    QVector<Molecule> mols;

    /** Hash indexing the molecules in the array, allowing
        rapid lookup by MoleculeID */
    QHash<MoleculeID, int> idx;

    /** The name of this group */
    QString nme;

    /** The ID and version number of this group */
    SireBase::IDMajMinVersion id_and_version;
};

/** Return the name of the group */
inline const QString& MoleculeGroupPvt::name() const
{
    return nme;
}

/** Return the ID of the group */
inline MoleculeGroupID MoleculeGroupPvt::ID() const
{
    return MoleculeGroupID(id_and_version.ID());
}

/** Return the version of the group */
inline const Version& MoleculeGroupPvt::version() const
{
    return id_and_version.version();
}

/** Return the array of all molecules in the group */
inline const QVector<Molecule>& MoleculeGroupPvt::molecules() const
{
    return mols;
}

/** Return whether or not we contain the molecule 'molecule' */
inline bool MoleculeGroupPvt::contains(MoleculeID molid) const
{
    return idx.contains(molid);
}

} //end of namespace detail

/** This class holds a group of molecules. The group is versioned,
    and holds the molecules in a manner that allows efficient linear
    access and indexing (by MoleculeID)

    @author Christopher Woods
*/
class SIREMOL_EXPORT MoleculeGroup
{

friend QDataStream& ::operator<<(QDataStream&, const MoleculeGroup&);
friend QDataStream& ::operator>>(QDataStream&, MoleculeGroup&);

public:
    MoleculeGroup();

    MoleculeGroup(const QString &name);
    MoleculeGroup(const QString &name, const QVector<Molecule> &molecules);

    MoleculeGroup(const MoleculeGroup &other);

    ~MoleculeGroup();

    MoleculeGroup& operator=(const MoleculeGroup &other);

    bool operator==(const MoleculeGroup &other) const;
    bool operator!=(const MoleculeGroup &other) const;

    const Molecule& operator[](MoleculeID molid) const;

    const Molecule& at(MoleculeID molid) const;

    const Molecule& molecule(MoleculeID molid) const;

    const Molecule& operator[](const Molecule &molecule) const;
    const Molecule& at(const Molecule &molecule) const;
    const Molecule& molecule(const Molecule &molecule) const;

    const QString& name() const;

    MoleculeGroupID ID() const;
    const Version& version() const;

    const QVector<Molecule>& molecules() const;

    QSet<MoleculeID> moleculeIDs() const;

    int count() const;

    bool add(const Molecule &molecule);
    bool change(const Molecule &molecule);
    bool remove(const Molecule &molecule);

    bool add(const QVector<Molecule> &molecules);
    bool change(const QVector<Molecule> &molecules);
    bool remove(const QVector<Molecule> &molecules);

    bool change(const QHash<MoleculeID,Molecule> &molecules);

    void rename(const QString &newname);

    bool contains(const Molecule &molecule) const;
    bool contains(MoleculeID molid) const;

private:
    /** Implicit pointer to the data of this object */
    QSharedDataPointer<detail::MoleculeGroupPvt> d;
};

/** Return the name of the group */
inline const QString& MoleculeGroup::name() const
{
    return d->name();
}

/** Return the molecule with ID == molid

    \throw SireMol::missing_molecule
*/
inline const Molecule& MoleculeGroup::operator[](MoleculeID molid) const
{
    return d->molecule(molid);
}

/** Return the molecule with ID == molid

    \throw SireMol::missing_molecule
*/
inline const Molecule& MoleculeGroup::at(MoleculeID molid) const
{
    return d->molecule(molid);
}

/** Return the molecule with ID == molid

    \throw SireMol::missing_molecule
*/
inline const Molecule& MoleculeGroup::molecule(MoleculeID molid) const
{
    return d->molecule(molid);
}

/** Return the copy of 'molecule' that is held in this group

    \throw SireMol::missing_molecule
*/
inline const Molecule& MoleculeGroup::operator[](const Molecule &molecule) const
{
    return d->molecule(molecule.ID());
}

/** Return the copy of 'molecule' that is held in this group

    \throw SireMol::missing_molecule
*/
inline const Molecule& MoleculeGroup::at(const Molecule &molecule) const
{
    return d->molecule(molecule.ID());
}

/** Return the set of IDs of all molecules in this group */
inline QSet<MoleculeID> MoleculeGroup::moleculeIDs() const
{
    return d->moleculeIDs();
}

/** Return the copy of 'molecule' that is held in this group

    \throw SireMol::missing_molecule
*/
inline const Molecule& MoleculeGroup::molecule(const Molecule &molecule) const
{
    return d->molecule(molecule.ID());
}

/** Return the ID number of this group */
inline MoleculeGroupID MoleculeGroup::ID() const
{
    return d->ID();
}

/** Return the version of this group */
inline const Version& MoleculeGroup::version() const
{
    return d->version();
}

/** Return the array of all molecules in this group */
inline const QVector<Molecule>& MoleculeGroup::molecules() const
{
    return d->molecules();
}

/** Return whether or not this group contains the molecule
    with ID == molid */
inline bool MoleculeGroup::contains(MoleculeID molid) const
{
    return d->contains(molid);
}

/** Return whether or not this group contains the molecule
    'molecule' */
inline bool MoleculeGroup::contains(const Molecule &molecule) const
{
    return d->contains(molecule.ID());
}

/** Return the number of molecules in the group */
inline int MoleculeGroup::count() const
{
    return molecules().count();
}

}

Q_DECLARE_METATYPE(SireMol::detail::MoleculeGroupPvt);
Q_DECLARE_METATYPE(SireMol::MoleculeGroup);

SIRE_END_HEADER

#endif
