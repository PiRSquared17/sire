// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "IntegerPower.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

void register_IntegerPower_class(){

    bp::class_< SireCAS::IntegerPower, bp::bases< SireCAS::ConstantPower > >( "IntegerPower" )    
        .def( bp::init< >() )    
        .def( bp::init< SireCAS::Expression const &, int >(( bp::arg("expression"), bp::arg("power") )) )    
        .def( 
            "evaluate"
            , (double ( ::SireCAS::IntegerPower::* )( ::SireCAS::Values const & ) const)( &::SireCAS::IntegerPower::evaluate )
            , ( bp::arg("values") ) )    
        .def( 
            "evaluate"
            , (::SireMaths::Complex ( ::SireCAS::IntegerPower::* )( ::SireCAS::ComplexValues const & ) const)( &::SireCAS::IntegerPower::evaluate )
            , ( bp::arg("values") ) )    
        .def( 
            "hash"
            , (::uint ( ::SireCAS::IntegerPower::* )(  ) const)( &::SireCAS::IntegerPower::hash ) )    
        .def( bp::self == bp::other< SireCAS::ExBase >() )    
        .def( 
            "power"
            , (::SireCAS::Expression ( ::SireCAS::IntegerPower::* )(  ) const)( &::SireCAS::IntegerPower::power ) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireCAS::IntegerPower::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireCAS::IntegerPower::* )(  ) const)( &::SireCAS::IntegerPower::what ) )    
        .staticmethod( "typeName" );

}
