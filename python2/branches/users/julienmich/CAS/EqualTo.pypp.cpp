// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "EqualTo.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/errors.h"

#include "SireError/errors.h"

#include "SireMaths/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "complexvalues.h"

#include "conditional.h"

#include "expressions.h"

#include "functions.h"

#include "identities.h"

#include "symbols.h"

#include "values.h"

#include "conditional.h"

SireCAS::EqualTo __copy__(const SireCAS::EqualTo &other){ return SireCAS::EqualTo(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_EqualTo_class(){

    { //::SireCAS::EqualTo
        typedef bp::class_< SireCAS::EqualTo, bp::bases< SireCAS::Condition, SireCAS::ExBase > > EqualTo_exposer_t;
        EqualTo_exposer_t EqualTo_exposer = EqualTo_exposer_t( "EqualTo" );
        bp::scope EqualTo_scope( EqualTo_exposer );
        EqualTo_exposer.def( bp::init< >() );
        EqualTo_exposer.def( bp::init< SireCAS::Expression const &, SireCAS::Expression const & >(( bp::arg("left_hand_side"), bp::arg("right_hand_side") )) );
        EqualTo_exposer.def( bp::init< SireCAS::EqualTo const & >(( bp::arg("other") )) );
        { //::SireCAS::EqualTo::alwaysFalse
        
            typedef bool ( ::SireCAS::EqualTo::*alwaysFalse_function_type )(  ) const;
            alwaysFalse_function_type alwaysFalse_function_value( &::SireCAS::EqualTo::alwaysFalse );
            
            EqualTo_exposer.def( 
                "alwaysFalse"
                , alwaysFalse_function_value );
        
        }
        { //::SireCAS::EqualTo::alwaysTrue
        
            typedef bool ( ::SireCAS::EqualTo::*alwaysTrue_function_type )(  ) const;
            alwaysTrue_function_type alwaysTrue_function_value( &::SireCAS::EqualTo::alwaysTrue );
            
            EqualTo_exposer.def( 
                "alwaysTrue"
                , alwaysTrue_function_value );
        
        }
        { //::SireCAS::EqualTo::operator=
        
            typedef ::SireCAS::EqualTo & ( ::SireCAS::EqualTo::*assign_function_type )( ::SireCAS::EqualTo const & ) ;
            assign_function_type assign_function_value( &::SireCAS::EqualTo::operator= );
            
            EqualTo_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        EqualTo_exposer.def( bp::self == bp::self );
        EqualTo_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::EqualTo::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::EqualTo::typeName );
            
            EqualTo_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireCAS::EqualTo::what
        
            typedef char const * ( ::SireCAS::EqualTo::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireCAS::EqualTo::what );
            
            EqualTo_exposer.def( 
                "what"
                , what_function_value );
        
        }
        EqualTo_exposer.staticmethod( "typeName" );
        EqualTo_exposer.def( "__copy__", &__copy__);
        EqualTo_exposer.def( "__deepcopy__", &__copy__);
        EqualTo_exposer.def( "clone", &__copy__);
        EqualTo_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::EqualTo >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        EqualTo_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::EqualTo >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        EqualTo_exposer.def( "__str__", &__str__< ::SireCAS::EqualTo > );
        EqualTo_exposer.def( "__repr__", &__str__< ::SireCAS::EqualTo > );
    }

}
