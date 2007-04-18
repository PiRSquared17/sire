// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ResNumIndexID.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMol/moleculedata.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

SireMol::ResNumIndexID __copy__(const SireMol::ResNumIndexID &other){ return SireMol::ResNumIndexID(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_ResNumIndexID_class(){

    bp::class_< SireMol::ResNumIndexID >( "ResNumIndexID", bp::init< bp::optional< SireMol::ResNum, SireMol::Index > >(( bp::arg("resnum")=::SireMol::ResNum( 0 ), bp::arg("idx")=::SireMol::Index( 0 ) )) )    
        .def( bp::init< boost::tuples::tuple<SireMol::ResNum, SireMol::Index, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type> const & >(( bp::arg("t") )) )    
        .def( 
            "index"
            , &::SireMol::ResNumIndexID::index )    
        .def( bp::self != bp::self )    
        .def( bp::self < bp::self )    
        .def( bp::self <= bp::self )    
        .def( bp::self == bp::self )    
        .def( bp::self > bp::self )    
        .def( bp::self >= bp::self )    
        .def( 
            "resNum"
            , &::SireMol::ResNumIndexID::resNum )    
        .def( 
            "toString"
            , &::SireMol::ResNumIndexID::toString )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::ResNumIndexID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::ResNumIndexID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireMol::ResNumIndexID > );

}
