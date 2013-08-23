// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Coth.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "exp.h"

#include "expression.h"

#include "hyperbolicfuncs.h"

#include "identities.h"

#include "trigfuncs.h"

#include "hyperbolicfuncs.h"

SireCAS::Coth __copy__(const SireCAS::Coth &other){ return SireCAS::Coth(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Coth_class(){

    { //::SireCAS::Coth
        typedef bp::class_< SireCAS::Coth, bp::bases< SireCAS::SingleFunc, SireCAS::ExBase > > Coth_exposer_t;
        Coth_exposer_t Coth_exposer = Coth_exposer_t( "Coth", bp::init< >() );
        bp::scope Coth_scope( Coth_exposer );
        Coth_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) );
        Coth_exposer.def( bp::init< SireCAS::Coth const & >(( bp::arg("other") )) );
        { //::SireCAS::Coth::evaluate
        
            typedef double ( ::SireCAS::Coth::*evaluate_function_type )( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Coth::evaluate );
            
            Coth_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::Coth::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::Coth::*evaluate_function_type )( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Coth::evaluate );
            
            Coth_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        Coth_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::Coth::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::Coth::typeName );
            
            Coth_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireCAS::Coth::what
        
            typedef char const * ( ::SireCAS::Coth::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireCAS::Coth::what );
            
            Coth_exposer.def( 
                "what"
                , what_function_value );
        
        }
        Coth_exposer.staticmethod( "typeName" );
        Coth_exposer.def( "__copy__", &__copy__);
        Coth_exposer.def( "__deepcopy__", &__copy__);
        Coth_exposer.def( "clone", &__copy__);
        Coth_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::Coth >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Coth_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::Coth >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Coth_exposer.def( "__str__", &__str__< ::SireCAS::Coth > );
        Coth_exposer.def( "__repr__", &__str__< ::SireCAS::Coth > );
    }

}