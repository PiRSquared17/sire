// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "EqualTo.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/errors.h"

#include "SireMaths/complex.h"

#include "SireMaths/errors.h"

#include "Siren/errors.h"

#include "Siren/objref.h"

#include "Siren/stream.h"

#include "complexvalues.h"

#include "conditional.h"

#include "factor.h"

#include "identities.h"

#include "values.h"

#include "conditional.h"

#include "Siren/logger.h"

SireCAS::EqualTo __copy__(const SireCAS::EqualTo &other){ return SireCAS::EqualTo(other); }

#include "Siren/str.hpp"

void register_EqualTo_class(){

    { //::SireCAS::EqualTo
        typedef bp::class_< SireCAS::EqualTo, bp::bases< SireCAS::Condition, SireCAS::CASNode, Siren::Object > > EqualTo_exposer_t;
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
        EqualTo_exposer.def( bp::self != bp::self );
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
        { //::SireCAS::EqualTo::stream
        
            typedef void ( ::SireCAS::EqualTo::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireCAS::EqualTo::stream );
            
            EqualTo_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        EqualTo_exposer.def( "__copy__", &__copy__);
        EqualTo_exposer.def( "__deepcopy__", &__copy__);
        EqualTo_exposer.def( "clone", &__copy__);
        EqualTo_exposer.def( "__str__", &__str__< ::SireCAS::EqualTo > );
        EqualTo_exposer.def( "__repr__", &__str__< ::SireCAS::EqualTo > );
    }

}
