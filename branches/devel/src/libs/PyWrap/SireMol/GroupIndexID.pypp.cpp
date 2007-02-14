// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "GroupIndexID.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_GroupIndexID_class(){

    bp::class_< SireMol::GroupIndexID >( "GroupIndexID", bp::init< bp::optional< SireMol::GroupID, SireMol::Index > >(( bp::arg("gid")=::SireMol::GroupID( 0 ), bp::arg("idx")=::SireMol::Index( 0 ) )) )    
        .def( bp::init< boost::tuples::tuple<SireMol::GroupID, SireMol::Index, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type> const & >(( bp::arg("t") )) )    
        .def( 
            "groupID"
            , &::SireMol::GroupIndexID::groupID )    
        .def( 
            "index"
            , &::SireMol::GroupIndexID::index )    
        .def( bp::self != bp::self )    
        .def( bp::self < bp::self )    
        .def( bp::self <= bp::self )    
        .def( bp::self == bp::self )    
        .def( bp::self > bp::self )    
        .def( bp::self >= bp::self )    
        .def( 
            "toString"
            , &::SireMol::GroupIndexID::toString );

}
