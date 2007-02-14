// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "IntegrationConstant.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

void register_IntegrationConstant_class(){

    bp::class_< SireCAS::IntegrationConstant, bp::bases< SireCAS::Symbol > >( "IntegrationConstant" )    
        .def( bp::init< >() )    
        .def( 
            "hash"
            , (::uint ( ::SireCAS::IntegrationConstant::* )(  ) const)( &::SireCAS::IntegrationConstant::hash ) )    
        .def( 
            "integrate"
            , (::SireCAS::Expression ( ::SireCAS::IntegrationConstant::* )( ::SireCAS::Symbol const & ) const)( &::SireCAS::IntegrationConstant::integrate )
            , ( bp::arg("symbol") ) )    
        .def( bp::self == bp::other< SireCAS::ExBase >() )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireCAS::IntegrationConstant::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireCAS::IntegrationConstant::* )(  ) const)( &::SireCAS::IntegrationConstant::what ) )    
        .staticmethod( "typeName" );

}
