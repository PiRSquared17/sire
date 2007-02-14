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

#include "SireCAS/qhash_sirecas.h"
#include "SireMol/qhash_siremol.h"
#include "SireMM/qhash_siremm.h"

#include "angledb.h"

#include "SireCAS/expression.h"
#include "SireCAS/values.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireMM;
using namespace SireDB;
using namespace SireCAS;

static const RegisterMetaType<AngleDB> r_angledb(NO_STREAM);

/** Constructor */
AngleDB::AngleDB() : Term13DB(), InternalDB<MolAngleInfo>(),
                     _theta("theta"), _r("r")
{}

/** Copy constructor */
AngleDB::AngleDB(const AngleDB &other)
        : Term13DB(other), InternalDB<MolAngleInfo>(other),
          _theta(other._theta), _r(other._r)
{}

/** Destructor */
AngleDB::~AngleDB()
{}

/** Dump the version number of this database */
void AngleDB::prepareToDump()
{
    Term13DB::prepareToDump();

    this->saveParameter<AngleDB>( "version" , 0 );
}

/** Check the version number */
void AngleDB::postLoad()
{
    Term13DB::postLoad();

    int v = this->loadParameter<AngleDB>( "version" ).toInt();

    if (v != 0)
        throw version_error( v, "0", "SireMM::AngleDB", CODELOC );
}

/** Add the angle function 'anglefunc' and associate it with the user identification
    string 'userid'. Note that 'anglefunc' must be a function of the variable
    representing the angle size ("theta", obtainable via Symbol("theta") or
    AngleDB::theta()), or the distance between the 1-3 atoms, ("r", via
    Symbol("r") or AngleDB::r()) or else it will be treated as a constant. */
void AngleDB::addAngle(const QString &userid, const Expression &anglefunc)
{
    if ( not (anglefunc.isFunction(_theta) or anglefunc.isFunction(_r)) )
    {
        //evaluate this constant...
        Expression const_val = anglefunc.evaluate(Values());

        ParamID paramid = addFunc(const_val);
        relateParameter(userid, paramid);
    }
    else
    {
        ParamID paramid = addFunc(anglefunc);
        relateParameter(userid, paramid);
    }
}

/** Return the angle function associated with the user identification ID 'userid',
    or a zero function if there is no such function. foundangle is set to whether
    a function was found. */
Expression AngleDB::getAngle(const QString &userid, bool *foundangle)
{
    ParamID paramid = getParameter(userid);

    if (foundangle)
        *foundangle = (paramid != 0);

    return retrieveFunc(paramid);
}

/** Relate the angle function associated with the user identification ID 'userid'
    with the relationship matching a triple of atoms 'matchangle' */
void AngleDB::relateAngle(const AssertMatch<3> &matchangle, const QString &userid)
{
    relateFunc(matchangle, userid);
}

/** Relate the angle function 'anglefunc' with the relationship matching a triple of atoms.
    Note that the anglefunc must be a function of "theta" or "r" or else the function
    will be treated as a constant */
void AngleDB::relateAngle(const AssertMatch<3> &matchangle, const Expression &anglefunc)
{
    if ( not (anglefunc.isFunction(_theta) or anglefunc.isFunction(_r)) )
    {
        Expression const_val = anglefunc.evaluate(Values());
        relateFunc(matchangle, const_val);
    }
    else
        relateFunc(matchangle, anglefunc);
}

/** Return the angle that matches the atom-triple matching relationship 'relateid',
    or a zero function if there is no match. foundangle is set to whether or
    not a function was found */
Expression AngleDB::getAngle(RelateID relateid, bool *foundangle)
{
    return getFunc(relateid, foundangle);
}

/** Return the angle that best matches the atom-triple matching relationships
    in 'relateids', or a zero function if there is no match. foundangle is set
    to whether or not a function was found. */
Expression AngleDB::getAngle(const RelateIDMap &relateids, bool *foundangle)
{
    return getFunc(relateids, foundangle);
}

/** Relate the angle parameter identified by the user ID string 'userid' to the matching criteria
    identified by the relationship ID 'relateid' */
void AngleDB::relateAngle(RelateID relateid, const QString &userid)
{
    relateParameter(relateid, userid);
}

/** Relate the angle function 'anglefunc' to the matching criteria identified by the relationship
    ID 'relateid' */
void AngleDB::relateAngle(RelateID relateid, const Expression &anglefunc)
{
    //add the function
    ParamID paramid = addFunc(anglefunc);
    //associate this relationship with the parameter
    relateParameter(relateid, paramid);
}

/** Assign the parameter for the angle 'angle' using the relationship IDs in
    'relateids', and place the parameter into the table 'param_table'. Return
    whether or not a parameter was found and assigned. */
bool AngleDB::assignParameter(const Angle &angle, const RelateIDMap &relateids,
                             ParameterTable &param_table)
{
/*    BOOST_ASSERT(param_table.isA<AngleTable>());

    bool found;
    Expression anglefunc = this->getFunc(relateids, &found);

    if (found)
        //add the expression to the parameter table
        param_table.asA<AngleTable>().setParameter(angle, anglefunc);

    return found;*/

    return false;
}
