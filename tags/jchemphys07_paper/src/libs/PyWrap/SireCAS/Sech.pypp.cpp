// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Sech.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

SireCAS::Sech __copy__(const SireCAS::Sech &other){ return SireCAS::Sech(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_Sech_class(){

    bp::class_< SireCAS::Sech, bp::bases< SireCAS::SingleFunc > >( "Sech" )    
        .def( bp::init< >() )    
        .def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) )    
        .def( 
            "evaluate"
            , (double ( ::SireCAS::Sech::* )( ::SireCAS::Values const & ) const)( &::SireCAS::Sech::evaluate )
            , ( bp::arg("values") ) )    
        .def( 
            "evaluate"
            , (::SireMaths::Complex ( ::SireCAS::Sech::* )( ::SireCAS::ComplexValues const & ) const)( &::SireCAS::Sech::evaluate )
            , ( bp::arg("values") ) )    
        .def( bp::self == bp::other< SireCAS::ExBase >() )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireCAS::Sech::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireCAS::Sech::* )(  ) const)( &::SireCAS::Sech::what ) )    
        .staticmethod( "typeName" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireCAS::Sech >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireCAS::Sech >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireCAS::Sech > );

}
