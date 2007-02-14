// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ArcCsch.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

void register_ArcCsch_class(){

    bp::class_< SireCAS::ArcCsch, bp::bases< SireCAS::SingleFunc > >( "ArcCsch" )    
        .def( bp::init< >() )    
        .def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) )    
        .def( 
            "evaluate"
            , (double ( ::SireCAS::ArcCsch::* )( ::SireCAS::Values const & ) const)( &::SireCAS::ArcCsch::evaluate )
            , ( bp::arg("values") ) )    
        .def( 
            "evaluate"
            , (::SireMaths::Complex ( ::SireCAS::ArcCsch::* )( ::SireCAS::ComplexValues const & ) const)( &::SireCAS::ArcCsch::evaluate )
            , ( bp::arg("values") ) )    
        .def( bp::self == bp::other< SireCAS::ExBase >() )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireCAS::ArcCsch::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireCAS::ArcCsch::* )(  ) const)( &::SireCAS::ArcCsch::what ) )    
        .staticmethod( "typeName" );

}
