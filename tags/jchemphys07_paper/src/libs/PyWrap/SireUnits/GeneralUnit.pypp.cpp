// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "GeneralUnit.pypp.hpp"
#include "boost/python.hpp"
#include "sireunits_headers.h"

namespace bp = boost::python;

SireUnits::Dimension::GeneralUnit __copy__(const SireUnits::Dimension::GeneralUnit &other){ return SireUnits::Dimension::GeneralUnit(other); }

#include "SirePy/str.hpp"

void register_GeneralUnit_class(){

    bp::class_< SireUnits::Dimension::GeneralUnit, bp::bases< SireUnits::Dimension::Unit > >( "GeneralUnit" )    
        .def( bp::init< >() )    
        .def( 
            "ANGLE"
            , &::SireUnits::Dimension::GeneralUnit::ANGLE )    
        .def( 
            "CHARGE"
            , &::SireUnits::Dimension::GeneralUnit::CHARGE )    
        .def( 
            "LENGTH"
            , &::SireUnits::Dimension::GeneralUnit::LENGTH )    
        .def( 
            "MASS"
            , &::SireUnits::Dimension::GeneralUnit::MASS )    
        .def( 
            "QUANTITY"
            , &::SireUnits::Dimension::GeneralUnit::QUANTITY )    
        .def( 
            "TEMPERATURE"
            , &::SireUnits::Dimension::GeneralUnit::TEMPERATURE )    
        .def( 
            "TIME"
            , &::SireUnits::Dimension::GeneralUnit::TIME )    
        .def( 
            "invert"
            , &::SireUnits::Dimension::GeneralUnit::invert )    
        .def( bp::self != bp::self )    
        .def( bp::self * bp::self )    
        .def( bp::self * bp::other< double >() )    
        .def( bp::self * bp::other< int >() )    
        .def( bp::self *= bp::self )    
        .def( bp::self *= bp::other< double >() )    
        .def( bp::self *= bp::other< int >() )    
        .def( bp::self + bp::self )    
        .def( bp::self += bp::self )    
        .def( -bp::self )    
        .def( bp::self - bp::self )    
        .def( bp::self -= bp::self )    
        .def( bp::self / bp::self )    
        .def( bp::self / bp::other< double >() )    
        .def( bp::self / bp::other< int >() )    
        .def( bp::self /= bp::self )    
        .def( bp::self /= bp::other< double >() )    
        .def( bp::self /= bp::other< int >() )    
        .def( bp::self == bp::self )    
        .def( 
            "to"
            , (double ( ::SireUnits::Dimension::GeneralUnit::* )( ::SireUnits::Dimension::TempBase const & ) const)( &::SireUnits::Dimension::GeneralUnit::to )
            , ( bp::arg("other") ) )    
        .def( 
            "to"
            , (double ( ::SireUnits::Dimension::GeneralUnit::* )( ::SireUnits::Dimension::GeneralUnit const & ) const)( &::SireUnits::Dimension::GeneralUnit::to )
            , ( bp::arg("other") ) )    
        .def( 
            "toString"
            , &::SireUnits::Dimension::GeneralUnit::toString )    
        .def( bp::other< int >() * bp::self )    
        .def( bp::other< double >() * bp::self )    
        .def( bp::other< int >() / bp::self )    
        .def( bp::other< double >() / bp::self )    
        .def( "__copy__", &__copy__)    
        .def( "__str__", &SirePy::__str__< ::SireUnits::Dimension::GeneralUnit > );

}
