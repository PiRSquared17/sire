// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "BondResID.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

void register_BondResID_class(){

    bp::class_< SireMM::BondResID >( "BondResID" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::ResNum >(( bp::arg("resnum") )) )    
        .def( bp::init< SireMol::ResNum, SireMol::ResNum >(( bp::arg("resnum0"), bp::arg("resnum1") )) )    
        .def( bp::init< SireMol::Bond const & >(( bp::arg("bond") )) )    
        .def( 
            "isInterResidue"
            , (bool ( ::SireMM::BondResID::* )(  ) const)( &::SireMM::BondResID::isInterResidue ) )    
        .def( 
            "isIntraResidue"
            , (bool ( ::SireMM::BondResID::* )(  ) const)( &::SireMM::BondResID::isIntraResidue ) )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "__getitem__"
            , (::SireMol::ResNum ( ::SireMM::BondResID::* )( int ) const)( &::SireMM::BondResID::operator[] )
            , ( bp::arg("i") ) )    
        .def( 
            "resNum0"
            , (::SireMol::ResNum ( ::SireMM::BondResID::* )(  ) const)( &::SireMM::BondResID::resNum0 ) )    
        .def( 
            "resNum1"
            , (::SireMol::ResNum ( ::SireMM::BondResID::* )(  ) const)( &::SireMM::BondResID::resNum1 ) )    
        .def( 
            "residueNumbers"
            , (::QSet<SireMol::ResNum> ( ::SireMM::BondResID::* )(  ) const)( &::SireMM::BondResID::residueNumbers ) )    
        .def( 
            "size"
            , (int ( ::SireMM::BondResID::* )(  ) const)( &::SireMM::BondResID::size ) );

}
