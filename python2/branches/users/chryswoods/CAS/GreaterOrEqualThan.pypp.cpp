// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "GreaterOrEqualThan.pypp.hpp"

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

SireCAS::GreaterOrEqualThan __copy__(const SireCAS::GreaterOrEqualThan &other){ return SireCAS::GreaterOrEqualThan(other); }

#include "Siren/str.hpp"

void register_GreaterOrEqualThan_class(){

    { //::SireCAS::GreaterOrEqualThan
        typedef bp::class_< SireCAS::GreaterOrEqualThan, bp::bases< SireCAS::Condition, SireCAS::CASNode, Siren::Object > > GreaterOrEqualThan_exposer_t;
        GreaterOrEqualThan_exposer_t GreaterOrEqualThan_exposer = GreaterOrEqualThan_exposer_t( "GreaterOrEqualThan" );
        bp::scope GreaterOrEqualThan_scope( GreaterOrEqualThan_exposer );
        GreaterOrEqualThan_exposer.def( bp::init< >() );
        GreaterOrEqualThan_exposer.def( bp::init< SireCAS::Expression const &, SireCAS::Expression const & >(( bp::arg("left_hand_side"), bp::arg("right_hand_side") )) );
        GreaterOrEqualThan_exposer.def( bp::init< SireCAS::GreaterOrEqualThan const & >(( bp::arg("other") )) );
        { //::SireCAS::GreaterOrEqualThan::alwaysFalse
        
            typedef bool ( ::SireCAS::GreaterOrEqualThan::*alwaysFalse_function_type )(  ) const;
            alwaysFalse_function_type alwaysFalse_function_value( &::SireCAS::GreaterOrEqualThan::alwaysFalse );
            
            GreaterOrEqualThan_exposer.def( 
                "alwaysFalse"
                , alwaysFalse_function_value );
        
        }
        { //::SireCAS::GreaterOrEqualThan::alwaysTrue
        
            typedef bool ( ::SireCAS::GreaterOrEqualThan::*alwaysTrue_function_type )(  ) const;
            alwaysTrue_function_type alwaysTrue_function_value( &::SireCAS::GreaterOrEqualThan::alwaysTrue );
            
            GreaterOrEqualThan_exposer.def( 
                "alwaysTrue"
                , alwaysTrue_function_value );
        
        }
        GreaterOrEqualThan_exposer.def( bp::self != bp::self );
        { //::SireCAS::GreaterOrEqualThan::operator=
        
            typedef ::SireCAS::GreaterOrEqualThan & ( ::SireCAS::GreaterOrEqualThan::*assign_function_type )( ::SireCAS::GreaterOrEqualThan const & ) ;
            assign_function_type assign_function_value( &::SireCAS::GreaterOrEqualThan::operator= );
            
            GreaterOrEqualThan_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        GreaterOrEqualThan_exposer.def( bp::self == bp::self );
        { //::SireCAS::GreaterOrEqualThan::stream
        
            typedef void ( ::SireCAS::GreaterOrEqualThan::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireCAS::GreaterOrEqualThan::stream );
            
            GreaterOrEqualThan_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        GreaterOrEqualThan_exposer.def( "__copy__", &__copy__);
        GreaterOrEqualThan_exposer.def( "__deepcopy__", &__copy__);
        GreaterOrEqualThan_exposer.def( "clone", &__copy__);
        GreaterOrEqualThan_exposer.def( "__str__", &__str__< ::SireCAS::GreaterOrEqualThan > );
        GreaterOrEqualThan_exposer.def( "__repr__", &__str__< ::SireCAS::GreaterOrEqualThan > );
    }

}
