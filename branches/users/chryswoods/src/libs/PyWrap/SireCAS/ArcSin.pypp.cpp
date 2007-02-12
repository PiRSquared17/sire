// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ArcSin.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

void register_ArcSin_class(){

    bp::class_< SireCAS::ArcSin, bp::bases< SireCAS::SingleFunc > >( "ArcSin" )    
        .def( bp::init< >() )    
        .def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) )    
        .def( 
            "evaluate"
            , (double ( ::SireCAS::ArcSin::* )( ::SireCAS::Values const & ) const)( &::SireCAS::ArcSin::evaluate )
            , ( bp::arg("values") ) )    
        .def( 
            "evaluate"
            , (::SireMaths::Complex ( ::SireCAS::ArcSin::* )( ::SireCAS::ComplexValues const & ) const)( &::SireCAS::ArcSin::evaluate )
            , ( bp::arg("values") ) )    
        .def( bp::self == bp::other< SireCAS::ExBase >() )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireCAS::ArcSin::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireCAS::ArcSin::* )(  ) const)( &::SireCAS::ArcSin::what ) )    
        .staticmethod( "typeName" );

}
