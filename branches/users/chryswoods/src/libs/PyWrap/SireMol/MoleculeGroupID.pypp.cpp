// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MoleculeGroupID.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

SireMol::MoleculeGroupID __copy__(const SireMol::MoleculeGroupID &other){ return SireMol::MoleculeGroupID(other); }

#include "SirePy/str.hpp"

void register_MoleculeGroupID_class(){

    bp::class_< SireMol::MoleculeGroupID, bp::bases< SireMol::IDBase > >( "MoleculeGroupID" )    
        .def( bp::init< >() )    
        .def( bp::init< quint32 >(( bp::arg("id") )) )    
        .def( "__copy__", &__copy__)    
        .def( "__str__", &SirePy::__str__< ::SireMol::MoleculeGroupID > );

}
