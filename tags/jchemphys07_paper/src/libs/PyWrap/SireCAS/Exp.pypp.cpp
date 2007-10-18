// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Exp.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

SireCAS::Exp __copy__(const SireCAS::Exp &other){ return SireCAS::Exp(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_Exp_class(){

    bp::class_< SireCAS::Exp, bp::bases< SireCAS::PowerFunction > >( "Exp" )    
        .def( bp::init< >() )    
        .def( bp::init< SireCAS::Expression const & >(( bp::arg("power") )) )    
        .def( 
            "core"
            , (::SireCAS::Expression ( ::SireCAS::Exp::* )(  ) const)( &::SireCAS::Exp::core ) )    
        .def( 
            "differentiate"
            , (::SireCAS::Expression ( ::SireCAS::Exp::* )( ::SireCAS::Symbol const & ) const)( &::SireCAS::Exp::differentiate )
            , ( bp::arg("symbol") ) )    
        .def( 
            "evaluate"
            , (double ( ::SireCAS::Exp::* )( ::SireCAS::Values const & ) const)( &::SireCAS::Exp::evaluate )
            , ( bp::arg("values") ) )    
        .def( 
            "evaluate"
            , (::SireMaths::Complex ( ::SireCAS::Exp::* )( ::SireCAS::ComplexValues const & ) const)( &::SireCAS::Exp::evaluate )
            , ( bp::arg("values") ) )    
        .def( 
            "hash"
            , (::uint ( ::SireCAS::Exp::* )(  ) const)( &::SireCAS::Exp::hash ) )    
        .def( 
            "integrate"
            , (::SireCAS::Expression ( ::SireCAS::Exp::* )( ::SireCAS::Symbol const & ) const)( &::SireCAS::Exp::integrate )
            , ( bp::arg("symbol") ) )    
        .def( bp::self == bp::other< SireCAS::ExBase >() )    
        .def( 
            "power"
            , (::SireCAS::Expression ( ::SireCAS::Exp::* )(  ) const)( &::SireCAS::Exp::power ) )    
        .def( 
            "toString"
            , (::QString ( ::SireCAS::Exp::* )(  ) const)( &::SireCAS::Exp::toString ) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireCAS::Exp::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireCAS::Exp::* )(  ) const)( &::SireCAS::Exp::what ) )    
        .staticmethod( "typeName" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireCAS::Exp >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireCAS::Exp >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireCAS::Exp > );

}
