// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ComplexValues.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/complex.h"

#include "Siren/stream.h"

#include "Siren/streamqt.h"

#include "Siren/tostring.h"

#include "complexvalues.h"

#include "symbol.h"

#include "complexvalues.h"

#include "Siren/logger.h"

SireCAS::ComplexValues __copy__(const SireCAS::ComplexValues &other){ return SireCAS::ComplexValues(other); }

#include "Siren/str.hpp"

void register_ComplexValues_class(){

    { //::SireCAS::ComplexValues
        typedef bp::class_< SireCAS::ComplexValues, bp::bases< Siren::Object > > ComplexValues_exposer_t;
        ComplexValues_exposer_t ComplexValues_exposer = ComplexValues_exposer_t( "ComplexValues" );
        bp::scope ComplexValues_scope( ComplexValues_exposer );
        ComplexValues_exposer.def( bp::init< >() );
        ComplexValues_exposer.def( bp::init< SireCAS::Symbol const &, SireMaths::Complex const & >(( bp::arg("symbol"), bp::arg("value") )) );
        ComplexValues_exposer.def( bp::init< QHash< SireCAS::Symbol, SireMaths::Complex > const & >(( bp::arg("values") )) );
        ComplexValues_exposer.def( bp::init< SireCAS::Values const & >(( bp::arg("other") )) );
        ComplexValues_exposer.def( bp::init< SireCAS::ComplexValues const & >(( bp::arg("other") )) );
        { //::SireCAS::ComplexValues::hashCode
        
            typedef ::uint ( ::SireCAS::ComplexValues::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::SireCAS::ComplexValues::hashCode );
            
            ComplexValues_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::SireCAS::ComplexValues::keys
        
            typedef ::QList< SireCAS::Symbol > ( ::SireCAS::ComplexValues::*keys_function_type )(  ) const;
            keys_function_type keys_function_value( &::SireCAS::ComplexValues::keys );
            
            ComplexValues_exposer.def( 
                "keys"
                , keys_function_value );
        
        }
        ComplexValues_exposer.def( bp::self != bp::self );
        ComplexValues_exposer.def( bp::self + bp::self );
        { //::SireCAS::ComplexValues::operator=
        
            typedef ::SireCAS::ComplexValues & ( ::SireCAS::ComplexValues::*assign_function_type )( ::SireCAS::ComplexValues const & ) ;
            assign_function_type assign_function_value( &::SireCAS::ComplexValues::operator= );
            
            ComplexValues_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        ComplexValues_exposer.def( bp::self == bp::self );
        { //::SireCAS::ComplexValues::stream
        
            typedef void ( ::SireCAS::ComplexValues::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireCAS::ComplexValues::stream );
            
            ComplexValues_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::SireCAS::ComplexValues::symbols
        
            typedef ::QList< SireCAS::Symbol > ( ::SireCAS::ComplexValues::*symbols_function_type )(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::ComplexValues::symbols );
            
            ComplexValues_exposer.def( 
                "symbols"
                , symbols_function_value );
        
        }
        { //::SireCAS::ComplexValues::toString
        
            typedef ::QString ( ::SireCAS::ComplexValues::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireCAS::ComplexValues::toString );
            
            ComplexValues_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireCAS::ComplexValues::value
        
            typedef ::SireMaths::Complex ( ::SireCAS::ComplexValues::*value_function_type )( ::SireCAS::Symbol const & ) const;
            value_function_type value_function_value( &::SireCAS::ComplexValues::value );
            
            ComplexValues_exposer.def( 
                "value"
                , value_function_value
                , ( bp::arg("sym") ) );
        
        }
        { //::SireCAS::ComplexValues::values
        
            typedef ::QHash< unsigned int, SireMaths::Complex > const & ( ::SireCAS::ComplexValues::*values_function_type )(  ) const;
            values_function_type values_function_value( &::SireCAS::ComplexValues::values );
            
            ComplexValues_exposer.def( 
                "values"
                , values_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        ComplexValues_exposer.def( "__copy__", &__copy__);
        ComplexValues_exposer.def( "__deepcopy__", &__copy__);
        ComplexValues_exposer.def( "clone", &__copy__);
        ComplexValues_exposer.def( "__str__", &__str__< ::SireCAS::ComplexValues > );
        ComplexValues_exposer.def( "__repr__", &__str__< ::SireCAS::ComplexValues > );
    }

}
