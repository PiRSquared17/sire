// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MoleculeCutting.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMol/moleculedata.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"
#include "SireVol/space.h"

namespace bp = boost::python;

SireMol::MoleculeCutting __copy__(const SireMol::MoleculeCutting &other){ return SireMol::MoleculeCutting(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMol::MoleculeCutting&){ return "SireMol::MoleculeCutting";}

void register_MoleculeCutting_class(){

    bp::class_< SireMol::MoleculeCutting, bp::bases< SireMol::CuttingFunctionBase > >( "MoleculeCutting" )    
        .def( bp::init< >() )    
        .def( 
            "__call__"
            , &::SireMol::MoleculeCutting::operator()
            , ( bp::arg("atom"), bp::arg("moldata") ) )    
        .def( 
            "typeName"
            , &::SireMol::MoleculeCutting::typeName )    
        .def( 
            "what"
            , &::SireMol::MoleculeCutting::what )    
        .staticmethod( "typeName" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::MoleculeCutting >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::MoleculeCutting >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
