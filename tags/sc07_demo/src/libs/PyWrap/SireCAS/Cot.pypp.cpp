// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Cot.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

SireCAS::Cot __copy__(const SireCAS::Cot &other){ return SireCAS::Cot(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_Cot_class(){

    bp::class_< SireCAS::Cot, bp::bases< SireCAS::SingleFunc > >( "Cot" )    
        .def( bp::init< >() )    
        .def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) )    
        .def( 
            "evaluate"
            , (double ( ::SireCAS::Cot::* )( ::SireCAS::Values const & ) const)( &::SireCAS::Cot::evaluate )
            , ( bp::arg("values") ) )    
        .def( 
            "evaluate"
            , (::SireMaths::Complex ( ::SireCAS::Cot::* )( ::SireCAS::ComplexValues const & ) const)( &::SireCAS::Cot::evaluate )
            , ( bp::arg("values") ) )    
        .def( bp::self == bp::other< SireCAS::ExBase >() )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireCAS::Cot::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireCAS::Cot::* )(  ) const)( &::SireCAS::Cot::what ) )    
        .staticmethod( "typeName" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireCAS::Cot >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireCAS::Cot >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireCAS::Cot > );

}
