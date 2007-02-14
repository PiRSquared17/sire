// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Cosh.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

void register_Cosh_class(){

    bp::class_< SireCAS::Cosh, bp::bases< SireCAS::SingleFunc > >( "Cosh" )    
        .def( bp::init< >() )    
        .def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) )    
        .def( 
            "evaluate"
            , (double ( ::SireCAS::Cosh::* )( ::SireCAS::Values const & ) const)( &::SireCAS::Cosh::evaluate )
            , ( bp::arg("values") ) )    
        .def( 
            "evaluate"
            , (::SireMaths::Complex ( ::SireCAS::Cosh::* )( ::SireCAS::ComplexValues const & ) const)( &::SireCAS::Cosh::evaluate )
            , ( bp::arg("values") ) )    
        .def( bp::self == bp::other< SireCAS::ExBase >() )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireCAS::Cosh::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireCAS::Cosh::* )(  ) const)( &::SireCAS::Cosh::what ) )    
        .staticmethod( "typeName" );

}
