// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Ln.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

void register_Ln_class(){

    bp::class_< SireCAS::Ln, bp::bases< SireCAS::SingleFunc > >( "Ln" )    
        .def( bp::init< >() )    
        .def( bp::init< SireCAS::Expression const & >(( bp::arg("expression") )) )    
        .def( 
            "evaluate"
            , (double ( ::SireCAS::Ln::* )( ::SireCAS::Values const & ) const)( &::SireCAS::Ln::evaluate )
            , ( bp::arg("values") ) )    
        .def( 
            "evaluate"
            , (::SireMaths::Complex ( ::SireCAS::Ln::* )( ::SireCAS::ComplexValues const & ) const)( &::SireCAS::Ln::evaluate )
            , ( bp::arg("values") ) )    
        .def( bp::self == bp::other< SireCAS::ExBase >() )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireCAS::Ln::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireCAS::Ln::* )(  ) const)( &::SireCAS::Ln::what ) )    
        .staticmethod( "typeName" );

}
