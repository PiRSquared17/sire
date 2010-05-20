// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "GreaterThan.pypp.hpp"

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

SireCAS::GreaterThan __copy__(const SireCAS::GreaterThan &other){ return SireCAS::GreaterThan(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_GreaterThan_class(){

    { //::SireCAS::GreaterThan
        typedef bp::class_< SireCAS::GreaterThan, bp::bases< SireCAS::Condition, SireCAS::ExBase > > GreaterThan_exposer_t;
        GreaterThan_exposer_t GreaterThan_exposer = GreaterThan_exposer_t( "GreaterThan" );
        bp::scope GreaterThan_scope( GreaterThan_exposer );
        GreaterThan_exposer.def( bp::init< >() );
        GreaterThan_exposer.def( bp::init< SireCAS::Expression const &, SireCAS::Expression const & >(( bp::arg("left_hand_side"), bp::arg("right_hand_side") )) );
        GreaterThan_exposer.def( bp::init< SireCAS::GreaterThan const & >(( bp::arg("other") )) );
        { //::SireCAS::GreaterThan::alwaysFalse
        
            typedef bool ( ::SireCAS::GreaterThan::*alwaysFalse_function_type )(  ) const;
            alwaysFalse_function_type alwaysFalse_function_value( &::SireCAS::GreaterThan::alwaysFalse );
            
            GreaterThan_exposer.def( 
                "alwaysFalse"
                , alwaysFalse_function_value );
        
        }
        { //::SireCAS::GreaterThan::alwaysTrue
        
            typedef bool ( ::SireCAS::GreaterThan::*alwaysTrue_function_type )(  ) const;
            alwaysTrue_function_type alwaysTrue_function_value( &::SireCAS::GreaterThan::alwaysTrue );
            
            GreaterThan_exposer.def( 
                "alwaysTrue"
                , alwaysTrue_function_value );
        
        }
        { //::SireCAS::GreaterThan::operator=
        
            typedef ::SireCAS::GreaterThan & ( ::SireCAS::GreaterThan::*assign_function_type )( ::SireCAS::GreaterThan const & ) ;
            assign_function_type assign_function_value( &::SireCAS::GreaterThan::operator= );
            
            GreaterThan_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        GreaterThan_exposer.def( bp::self == bp::self );
        GreaterThan_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::GreaterThan::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::GreaterThan::typeName );
            
            GreaterThan_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireCAS::GreaterThan::what
        
            typedef char const * ( ::SireCAS::GreaterThan::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireCAS::GreaterThan::what );
            
            GreaterThan_exposer.def( 
                "what"
                , what_function_value );
        
        }
        GreaterThan_exposer.staticmethod( "typeName" );
        GreaterThan_exposer.def( "__copy__", &__copy__);
        GreaterThan_exposer.def( "__deepcopy__", &__copy__);
        GreaterThan_exposer.def( "clone", &__copy__);
        GreaterThan_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::GreaterThan >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GreaterThan_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::GreaterThan >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GreaterThan_exposer.def( "__str__", &__str__< ::SireCAS::GreaterThan > );
        GreaterThan_exposer.def( "__repr__", &__str__< ::SireCAS::GreaterThan > );
    }

}
