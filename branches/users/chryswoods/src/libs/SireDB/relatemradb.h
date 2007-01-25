/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREDB_RELATEMRADB_H
#define SIREDB_RELATEMRADB_H

#include "SireMol/atomindex.h"
#include "SireMol/molecule.h"

#include "relationshipdb.h"
#include "relateidmap.h"
#include "relatemradatalist.h"
#include "nameidset.h"

#include <QCache>

SIRE_BEGIN_HEADER

namespace SireDB
{

using SireMol::AtomIndex;
using SireMol::Molecule;
using SireMol::CutGroupID;
using SireMol::ResNum;

class MatchMRAData;
class MatchMRA;

class NameDB;
class RelateMRDB;
class RelateMRADB;

/**
This class provides a database that is used to hold atom matching criteria (relationships), which use matching criteria based on the molecule, residue and atom (MRA), e.g. match atoms called "CA1" in residue called "ALA". This class also contains the code necessary to find all of the relationships that match the atoms in a molecule. Each relationship is given a unique ID number (a RelateID) which is can be matched to parameters by other parts of the database, e.g. atoms called "CA1" in residues called "ALA" could have RelateID == 3, which could be used by the CLJDB to find a matching CLJ parameter.

@author Christopher Woods
*/
class SIREDB_EXPORT RelateMRADB : public RelationshipDB
{

public:
    RelateMRADB();
    RelateMRADB(const RelateMRADB &other);

    ~RelateMRADB();

    /** Return the full name of this database component */
    static const char* typeName()
    {
        return "SireDB::RelateMRADB";
    }

    const char* what() const
    {
        return RelateMRADB::typeName();
    }

    RelateID add(const MatchMRAData &match, uint n);
    RelateID get(const MatchMRAData &match, uint n);

    RelateID add(const MatchMRA &match);
    RelateID get(const MatchMRA &match);

    MatchMRAData get(RelateID relateid, uint n);

protected:
    void initialise();
    void prepareToDump();
    void postLoad();

    RelateIDMap findMatches(const MatchMRAData &match, uint n);

    RelateIDMap findMatches(const QList<AtomIndex> &atoms,
                            const Molecule &molecule,
                            const ParameterTable &param_table,
                            const MatchMRData &matchmr);

private:

    static QString relateTable(uint n);
    static QString queryString(const RelateIDMap &molresmap, uint i);

    void createRelateTable(uint n);

    RelateMRAData convert(const MatchMRA &data);
    RelateMRADataList convert(const MatchMRAData &data, uint n);

    MatchMRA convert(const RelateMRAData &data);
    MatchMRAData convert(const RelateMRADataList &data);

    RelateID add(const RelateMRADataList &relatedata);
    RelateID get(const RelateMRADataList &relatedata);

    /** Set of values of 'n' for the relationship tables that exist
        in this database */
    QSet<quint32> relatetables;

    /** Reference to the RelateMRDB that is used to relate molecule and residue
        information to a RelateID */
    RelateMRDB *relatemrdb;

    /** Reference to the NameDB that is used to relate atom names
        to NameID numbers */
    NameDB *namedb;
};

}

Q_DECLARE_METATYPE(SireDB::RelateMRADB)

SIRE_END_HEADER

#endif
