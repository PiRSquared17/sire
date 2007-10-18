// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Csch.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

SireCAS::Csch __copy__(const SireCAS::Csch &other){ return SireCAS::Csch(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_Csch_class(){

    bp::class_< SireCAS::Csch, bp::bases< SireCAS::SingleFunc > >( "Csch" )    
        .def( bp::init< >() )    
        .def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) )    
        .def( 
            "evaluate"
            , (double ( ::SireCAS::Csch::* )( ::SireCAS::Values const & ) const)( &::SireCAS::Csch::evaluate )
            , ( bp::arg("values") ) )    
        .def( 
            "evaluate"
            , (::SireMaths::Complex ( ::SireCAS::Csch::* )( ::SireCAS::ComplexValues const & ) const)( &::SireCAS::Csch::evaluate )
            , ( bp::arg("values") ) )    
        .def( bp::self == bp::other< SireCAS::ExBase >() )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireCAS::Csch::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireCAS::Csch::* )(  ) const)( &::SireCAS::Csch::what ) )    
        .staticmethod( "typeName" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireCAS::Csch >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireCAS::Csch >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireCAS::Csch > );

}
