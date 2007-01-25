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

/**
  * This file contains the boost::python wrapping of MatchMR
  *
  * @author Christopher Woods
  *
  */

#include <Python.h>
#include <boost/python.hpp>

#include "SireDB/matchdata.h"
#include "SireDB/matchatom.h"
#include "SireDB/matchres.h"
#include "SireDB/matchmol.h"
#include "SireDB/matchatomtype.h"
#include "SireDB/matchatomtypedata.h"
#include "SireDB/matchmr.h"
#include "SireDB/matchmrdata.h"
#include "SireDB/matchmra.h"
#include "SireDB/matchmradata.h"

using namespace boost::python;

namespace SireDB
{

MatchMR mr_and_res(const MatchMR &mr, const MatchRes &res)
{
    return mr and res;
}

MatchMR mr_and_mol(const MatchMR &mr, const MatchMol &mol)
{
    return mr and mol;
}

MatchMR res_and_mr(const MatchMR &mr, const MatchRes &res)
{
    return res and mr;
}

MatchMR mol_and_mr(const MatchMR &mr, const MatchMol &mol)
{
    return mol and mr;
}

MatchMRData mrdata_and_res(const MatchMRData &mrdata, const MatchRes &res)
{
    return mrdata and res;
}

MatchMRData mrdata_and_mol(const MatchMRData &mrdata, const MatchMol &mol)
{
    return mrdata and mol;
}

MatchMRData mrdata_and_mr(const MatchMRData &mrdata, const MatchMR &mr)
{
    return mrdata and mr;
}

MatchMRData res_and_mrdata(const MatchMRData &mrdata, const MatchRes &res)
{
    return res and mrdata;
}

MatchMRData mol_and_mrdata(const MatchMRData &mrdata, const MatchMol &mol)
{
    return mol and mrdata;
}

MatchMRData mr_and_mrdata(const MatchMRData &mrdata, const MatchMR &mr)
{
    return mr and mrdata;
}

void  
SIREDB_EXPORT
export_MatchMR()
{
    class_<MatchMR>("MatchMR", init<>())
        .def( init<const MatchMol&>() )
        .def( init<const MatchRes&>() )
        .def( init<const MatchMol&, const MatchRes&>() )
        .def( "__str__", &MatchMR::toString )
        .def( "toString", &MatchMR::toString )
        .def( "__and__", &MatchMR::operator&& )
        .def( "__and__", &mr_and_mol )
        .def( "__and__", &mr_and_res )
        .def( "__rand__", &mol_and_mr )
        .def( "__rand__", &res_and_mr )
        .def( "isEmpty", &MatchMR::isEmpty )
        .def( "molecule", &MatchMR::molecule, return_value_policy<copy_const_reference>() )
        .def( "residue", &MatchMR::residue, return_value_policy<copy_const_reference>() )
        .def( "hasMoleculeCriteria", &MatchMR::hasMoleculeCriteria )
        .def( "hasResidueCriteria", &MatchMR::hasResidueCriteria )
    ;
    
    implicitly_convertible<MatchMR, MatchMRData>();
    
    class_<MatchMRData>("MatchMRData", init<>())
        .def( init<const MatchMol&>() )
        .def( init<const MatchRes&>() )
        .def( init<const MatchMol&, const MatchRes&>() )
        .def( init<uint, const MatchRes&>() )
        .def( "__str__", &MatchMRData::toString )
        .def( "toString", &MatchMRData::toString )
        .def( "__and__", &MatchMRData::operator&& )
        .def( "__and__", &mrdata_and_res )
        .def( "__and__", &mrdata_and_mol )
        .def( "__and__", &mrdata_and_mr )
        .def( "__rand__", &res_and_mrdata )
        .def( "__rand__", &mol_and_mrdata )
        .def( "__rand__", &mr_and_mrdata )
        .def( "isEmpty", &MatchMRData::isEmpty )
        .def( "getMatch", &MatchMRData::getMatch )
        .def( "nMatches", &MatchMRData::nMatches )
        .def( "hasMoleculeCriteria", &MatchMRData::hasMoleculeCriteria )
        .def( "hasResidueCriteria", &MatchMRData::hasResidueCriteria )
        .def( "set", &MatchMRData::set )
        .def( "simplify", &MatchMRData::simplify )
    ;
    
    implicitly_convertible<MatchMRData, MatchMRAData>();

}

}
