// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "RigidBodyMC.pypp.hpp"
#include "boost/python.hpp"
#include "siremove_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/moleculegroup.h"
#include "SireSystem/simsystem.h"

namespace bp = boost::python;

SireMove::RigidBodyMC __copy__(const SireMove::RigidBodyMC &other){ return SireMove::RigidBodyMC(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMove::RigidBodyMC&){ return "SireMove::RigidBodyMC";}

void register_RigidBodyMC_class(){

    bp::class_< SireMove::RigidBodyMC, bp::bases< SireMove::MonteCarlo > >( "RigidBodyMC", bp::init< bp::optional< SireMaths::RanGenerator const & > >(( bp::arg("generator")=::SireMaths::RanGenerator( ) )) )    
        .def( bp::init< SireMove::Sampler const &, bp::optional< SireMaths::RanGenerator const & > >(( bp::arg("sampler"), bp::arg("generator")=::SireMaths::RanGenerator( ) )) )    
        .def( bp::init< SireMove::Sampler const &, SireUnits::Dimension::Length, SireUnits::Dimension::Angle, bp::optional< SireMaths::RanGenerator const & > >(( bp::arg("sampler"), bp::arg("max_translation"), bp::arg("max_rotation"), bp::arg("generator")=::SireMaths::RanGenerator( ) )) )    
        .def( 
            "assertCompatibleWith"
            , (void ( ::SireMove::RigidBodyMC::* )( ::SireSystem::QuerySystem & ) const)( &::SireMove::RigidBodyMC::assertCompatibleWith )
            , ( bp::arg("system") ) )    
        .def( 
            "maximumRotation"
            , (::SireUnits::Dimension::Angle ( ::SireMove::RigidBodyMC::* )(  ) const)( &::SireMove::RigidBodyMC::maximumRotation ) )    
        .def( 
            "maximumTranslation"
            , (::SireUnits::Dimension::Length ( ::SireMove::RigidBodyMC::* )(  ) const)( &::SireMove::RigidBodyMC::maximumTranslation ) )    
        .def( 
            "move"
            , (void ( ::SireMove::RigidBodyMC::* )( ::SireSystem::SimSystem & ) )( &::SireMove::RigidBodyMC::move )
            , ( bp::arg("system") ) )    
        .def( 
            "setMaximumRotation"
            , (void ( ::SireMove::RigidBodyMC::* )( ::SireUnits::Dimension::Angle ) )( &::SireMove::RigidBodyMC::setMaximumRotation )
            , ( bp::arg("max_rotation") ) )    
        .def( 
            "setMaximumTranslation"
            , (void ( ::SireMove::RigidBodyMC::* )( ::SireUnits::Dimension::Length ) )( &::SireMove::RigidBodyMC::setMaximumTranslation )
            , ( bp::arg("max_translation") ) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireMove::RigidBodyMC::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireMove::RigidBodyMC::* )(  ) const)( &::SireMove::RigidBodyMC::what ) )    
        .staticmethod( "typeName" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMove::RigidBodyMC >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMove::RigidBodyMC >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
