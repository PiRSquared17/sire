// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ResIDAtomID.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

SireMol::ResIDAtomID __copy__(const SireMol::ResIDAtomID &other){ return SireMol::ResIDAtomID(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_ResIDAtomID_class(){

    bp::class_< SireMol::ResIDAtomID >( "ResIDAtomID", bp::init< bp::optional< SireMol::ResID, SireMol::AtomID > >(( bp::arg("resid")=::SireMol::ResID( 0 ), bp::arg("atmid")=::SireMol::AtomID( 0 ) )) )    
        .def( bp::init< boost::tuples::tuple<SireMol::ResID, SireMol::AtomID, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type> const & >(( bp::arg("t") )) )    
        .def( 
            "atomID"
            , &::SireMol::ResIDAtomID::atomID )    
        .def( bp::self != bp::self )    
        .def( bp::self < bp::self )    
        .def( bp::self <= bp::self )    
        .def( bp::self == bp::self )    
        .def( bp::self > bp::self )    
        .def( bp::self >= bp::self )    
        .def( 
            "resID"
            , &::SireMol::ResIDAtomID::resID )    
        .def( 
            "toString"
            , &::SireMol::ResIDAtomID::toString )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::ResIDAtomID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::ResIDAtomID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireMol::ResIDAtomID > );

}
