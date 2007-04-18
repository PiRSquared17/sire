// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ResNumAtomID.pypp.hpp"
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

SireMol::ResNumAtomID __copy__(const SireMol::ResNumAtomID &other){ return SireMol::ResNumAtomID(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_ResNumAtomID_class(){

    bp::class_< SireMol::ResNumAtomID >( "ResNumAtomID", bp::init< bp::optional< SireMol::ResNum, SireMol::AtomID > >(( bp::arg("resnum")=::SireMol::ResNum( 0 ), bp::arg("atmid")=::SireMol::AtomID( 0 ) )) )    
        .def( bp::init< boost::tuples::tuple<SireMol::ResNum, SireMol::AtomID, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type> const & >(( bp::arg("t") )) )    
        .def( 
            "atomID"
            , &::SireMol::ResNumAtomID::atomID )    
        .def( bp::self != bp::self )    
        .def( bp::self < bp::self )    
        .def( bp::self <= bp::self )    
        .def( bp::self == bp::self )    
        .def( bp::self > bp::self )    
        .def( bp::self >= bp::self )    
        .def( 
            "resNum"
            , &::SireMol::ResNumAtomID::resNum )    
        .def( 
            "toString"
            , &::SireMol::ResNumAtomID::toString )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::ResNumAtomID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::ResNumAtomID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireMol::ResNumAtomID > );

}
