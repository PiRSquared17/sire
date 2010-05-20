// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Tanh.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "exp.h"

#include "expression.h"

#include "hyperbolicfuncs.h"

#include "identities.h"

#include "trigfuncs.h"

#include "hyperbolicfuncs.h"

SireCAS::Tanh __copy__(const SireCAS::Tanh &other){ return SireCAS::Tanh(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Tanh_class(){

    { //::SireCAS::Tanh
        typedef bp::class_< SireCAS::Tanh, bp::bases< SireCAS::SingleFunc, SireCAS::ExBase > > Tanh_exposer_t;
        Tanh_exposer_t Tanh_exposer = Tanh_exposer_t( "Tanh" );
        bp::scope Tanh_scope( Tanh_exposer );
        Tanh_exposer.def( bp::init< >() );
        Tanh_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) );
        Tanh_exposer.def( bp::init< SireCAS::Tanh const & >(( bp::arg("other") )) );
        { //::SireCAS::Tanh::evaluate
        
            typedef double ( ::SireCAS::Tanh::*evaluate_function_type )( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Tanh::evaluate );
            
            Tanh_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::Tanh::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::Tanh::*evaluate_function_type )( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Tanh::evaluate );
            
            Tanh_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        Tanh_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::Tanh::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::Tanh::typeName );
            
            Tanh_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireCAS::Tanh::what
        
            typedef char const * ( ::SireCAS::Tanh::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireCAS::Tanh::what );
            
            Tanh_exposer.def( 
                "what"
                , what_function_value );
        
        }
        Tanh_exposer.staticmethod( "typeName" );
        Tanh_exposer.def( "__copy__", &__copy__);
        Tanh_exposer.def( "__deepcopy__", &__copy__);
        Tanh_exposer.def( "clone", &__copy__);
        Tanh_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::Tanh >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Tanh_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::Tanh >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Tanh_exposer.def( "__str__", &__str__< ::SireCAS::Tanh > );
        Tanh_exposer.def( "__repr__", &__str__< ::SireCAS::Tanh > );
    }

}
