// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Identities.pypp.hpp"

namespace bp = boost::python;

#include "Siren/stream.h"

#include "Siren/streamqt.h"

#include "expression.h"

#include "identities.h"

#include "symbol.h"

#include "identities.h"

#include "Siren/logger.h"

SireCAS::Identities __copy__(const SireCAS::Identities &other){ return SireCAS::Identities(other); }

#include "Siren/str.hpp"

void register_Identities_class(){

    { //::SireCAS::Identities
        typedef bp::class_< SireCAS::Identities, bp::bases< Siren::Object > > Identities_exposer_t;
        Identities_exposer_t Identities_exposer = Identities_exposer_t( "Identities" );
        bp::scope Identities_scope( Identities_exposer );
        Identities_exposer.def( bp::init< >() );
        Identities_exposer.def( bp::init< SireCAS::Symbol const &, SireCAS::Expression const & >(( bp::arg("symbol"), bp::arg("expression") )) );
        Identities_exposer.def( bp::init< QHash< SireCAS::Symbol, SireCAS::Expression > const & >(( bp::arg("expressions") )) );
        Identities_exposer.def( bp::init< SireCAS::Identities const & >(( bp::arg("other") )) );
        { //::SireCAS::Identities::contains
        
            typedef bool ( ::SireCAS::Identities::*contains_function_type )( ::SireCAS::Symbol const & ) const;
            contains_function_type contains_function_value( &::SireCAS::Identities::contains );
            
            Identities_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("sym") ) );
        
        }
        { //::SireCAS::Identities::expression
        
            typedef ::SireCAS::Expression ( ::SireCAS::Identities::*expression_function_type )( ::SireCAS::Symbol const & ) const;
            expression_function_type expression_function_value( &::SireCAS::Identities::expression );
            
            Identities_exposer.def( 
                "expression"
                , expression_function_value
                , ( bp::arg("sym") ) );
        
        }
        { //::SireCAS::Identities::hashCode
        
            typedef ::uint ( ::SireCAS::Identities::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::SireCAS::Identities::hashCode );
            
            Identities_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::SireCAS::Identities::keys
        
            typedef ::QList< SireCAS::Symbol > ( ::SireCAS::Identities::*keys_function_type )(  ) const;
            keys_function_type keys_function_value( &::SireCAS::Identities::keys );
            
            Identities_exposer.def( 
                "keys"
                , keys_function_value );
        
        }
        Identities_exposer.def( bp::self != bp::self );
        Identities_exposer.def( bp::self + bp::self );
        { //::SireCAS::Identities::operator=
        
            typedef ::SireCAS::Identities & ( ::SireCAS::Identities::*assign_function_type )( ::SireCAS::Identities const & ) ;
            assign_function_type assign_function_value( &::SireCAS::Identities::operator= );
            
            Identities_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Identities_exposer.def( bp::self == bp::self );
        { //::SireCAS::Identities::stream
        
            typedef void ( ::SireCAS::Identities::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireCAS::Identities::stream );
            
            Identities_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::SireCAS::Identities::symbols
        
            typedef ::QList< SireCAS::Symbol > ( ::SireCAS::Identities::*symbols_function_type )(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::Identities::symbols );
            
            Identities_exposer.def( 
                "symbols"
                , symbols_function_value );
        
        }
        { //::SireCAS::Identities::toString
        
            typedef ::QString ( ::SireCAS::Identities::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireCAS::Identities::toString );
            
            Identities_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireCAS::Identities::values
        
            typedef ::QHash< unsigned int, SireCAS::Expression > const & ( ::SireCAS::Identities::*values_function_type )(  ) const;
            values_function_type values_function_value( &::SireCAS::Identities::values );
            
            Identities_exposer.def( 
                "values"
                , values_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        Identities_exposer.def( "__copy__", &__copy__);
        Identities_exposer.def( "__deepcopy__", &__copy__);
        Identities_exposer.def( "clone", &__copy__);
        Identities_exposer.def( "__str__", &__str__< ::SireCAS::Identities > );
        Identities_exposer.def( "__repr__", &__str__< ::SireCAS::Identities > );
    }

}
