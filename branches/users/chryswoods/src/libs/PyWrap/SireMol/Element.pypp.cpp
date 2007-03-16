// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Element.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

SireMol::Element __copy__(const SireMol::Element &other){ return SireMol::Element(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_Element_class(){

    bp::class_< SireMol::Element >( "Element" )    
        .def( bp::init< >() )    
        .def( bp::init< QString >(( bp::arg("element") )) )    
        .def( bp::init< char const * >(( bp::arg("element") )) )    
        .def( bp::init< unsigned int >(( bp::arg("nprotons") )) )    
        .def( bp::init< int >(( bp::arg("nprotons") )) )    
        .def( 
            "actinide"
            , &::SireMol::Element::actinide )    
        .def( 
            "alkaliEarthMetal"
            , &::SireMol::Element::alkaliEarthMetal )    
        .def( 
            "alkaliMetal"
            , &::SireMol::Element::alkaliMetal )    
        .def( 
            "biological"
            , &::SireMol::Element::biological )    
        .def( 
            "biologicalElement"
            , &::SireMol::Element::biologicalElement
            , ( bp::arg("name") ) )    
        .def( 
            "blue"
            , &::SireMol::Element::blue )    
        .def( 
            "bondOrderRadius"
            , &::SireMol::Element::bondOrderRadius )    
        .def( 
            "covalentRadius"
            , &::SireMol::Element::covalentRadius )    
        .def( 
            "electroNegativity"
            , &::SireMol::Element::electroNegativity )    
        .def( 
            "elementWithMass"
            , &::SireMol::Element::elementWithMass
            , ( bp::arg("mass") ) )    
        .def( 
            "green"
            , &::SireMol::Element::green )    
        .def( 
            "group"
            , &::SireMol::Element::group )    
        .def( 
            "halogen"
            , &::SireMol::Element::halogen )    
        .def( 
            "lanthanide"
            , &::SireMol::Element::lanthanide )    
        .def( 
            "mass"
            , &::SireMol::Element::mass )    
        .def( 
            "maxBonds"
            , &::SireMol::Element::maxBonds )    
        .def( 
            "nProtons"
            , &::SireMol::Element::nProtons )    
        .def( 
            "name"
            , &::SireMol::Element::name )    
        .def( 
            "nobleGas"
            , &::SireMol::Element::nobleGas )    
        .def( bp::self != bp::self )    
        .def( bp::self < bp::self )    
        .def( bp::self <= bp::self )    
        .def( bp::self == bp::self )    
        .def( bp::self > bp::self )    
        .def( bp::self >= bp::self )    
        .def( 
            "period"
            , &::SireMol::Element::period )    
        .def( 
            "rareEarth"
            , &::SireMol::Element::rareEarth )    
        .def( 
            "red"
            , &::SireMol::Element::red )    
        .def( 
            "symbol"
            , &::SireMol::Element::symbol )    
        .def( 
            "toString"
            , &::SireMol::Element::toString )    
        .def( 
            "transitionMetal"
            , &::SireMol::Element::transitionMetal )    
        .def( 
            "vdwRadius"
            , &::SireMol::Element::vdwRadius )    
        .staticmethod( "biologicalElement" )    
        .staticmethod( "elementWithMass" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::Element >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::Element >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireMol::Element > );

}
