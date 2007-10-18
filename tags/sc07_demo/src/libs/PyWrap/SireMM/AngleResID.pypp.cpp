// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "AngleResID.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

SireMM::AngleResID __copy__(const SireMM::AngleResID &other){ return SireMM::AngleResID(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::AngleResID&){ return "SireMM::AngleResID";}

void register_AngleResID_class(){

    bp::class_< SireMM::AngleResID >( "AngleResID" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::ResNum >(( bp::arg("resnum") )) )    
        .def( bp::init< SireMol::ResNum, SireMol::ResNum, SireMol::ResNum >(( bp::arg("resnum0"), bp::arg("resnum1"), bp::arg("resnum2") )) )    
        .def( bp::init< SireMol::Angle const & >(( bp::arg("angle") )) )    
        .def( 
            "isInterResidue"
            , (bool ( ::SireMM::AngleResID::* )(  ) const)( &::SireMM::AngleResID::isInterResidue ) )    
        .def( 
            "isIntraResidue"
            , (bool ( ::SireMM::AngleResID::* )(  ) const)( &::SireMM::AngleResID::isIntraResidue ) )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "__getitem__"
            , (::SireMol::ResNum ( ::SireMM::AngleResID::* )( int ) const)( &::SireMM::AngleResID::operator[] )
            , ( bp::arg("i") ) )    
        .def( 
            "resNum0"
            , (::SireMol::ResNum ( ::SireMM::AngleResID::* )(  ) const)( &::SireMM::AngleResID::resNum0 ) )    
        .def( 
            "resNum1"
            , (::SireMol::ResNum ( ::SireMM::AngleResID::* )(  ) const)( &::SireMM::AngleResID::resNum1 ) )    
        .def( 
            "resNum2"
            , (::SireMol::ResNum ( ::SireMM::AngleResID::* )(  ) const)( &::SireMM::AngleResID::resNum2 ) )    
        .def( 
            "residueNumbers"
            , (::QSet<SireMol::ResNum> ( ::SireMM::AngleResID::* )(  ) const)( &::SireMM::AngleResID::residueNumbers ) )    
        .def( 
            "size"
            , (int ( ::SireMM::AngleResID::* )(  ) const)( &::SireMM::AngleResID::size ) )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMM::AngleResID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMM::AngleResID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
