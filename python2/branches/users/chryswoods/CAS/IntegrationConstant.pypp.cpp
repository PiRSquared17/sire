// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "IntegrationConstant.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/errors.h"

#include "SireMaths/complex.h"

#include "Siren/stream.h"

#include "expression.h"

#include "integrationconstant.h"

#include "integrationconstant.h"

#include "Siren/logger.h"

SireCAS::IntegrationConstant __copy__(const SireCAS::IntegrationConstant &other){ return SireCAS::IntegrationConstant(other); }

#include "Siren/str.hpp"

void register_IntegrationConstant_class(){

    { //::SireCAS::IntegrationConstant
        typedef bp::class_< SireCAS::IntegrationConstant, bp::bases< SireCAS::Symbol, SireCAS::CASNode, Siren::Object > > IntegrationConstant_exposer_t;
        IntegrationConstant_exposer_t IntegrationConstant_exposer = IntegrationConstant_exposer_t( "IntegrationConstant" );
        bp::scope IntegrationConstant_scope( IntegrationConstant_exposer );
        IntegrationConstant_exposer.def( bp::init< >() );
        IntegrationConstant_exposer.def( bp::init< SireCAS::IntegrationConstant const & >(( bp::arg("other") )) );
        { //::SireCAS::IntegrationConstant::integrate
        
            typedef ::SireCAS::Expression ( ::SireCAS::IntegrationConstant::*integrate_function_type )( ::SireCAS::Symbol const & ) const;
            integrate_function_type integrate_function_value( &::SireCAS::IntegrationConstant::integrate );
            
            IntegrationConstant_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("symbol") ) );
        
        }
        IntegrationConstant_exposer.def( bp::self != bp::self );
        { //::SireCAS::IntegrationConstant::operator=
        
            typedef ::SireCAS::IntegrationConstant & ( ::SireCAS::IntegrationConstant::*assign_function_type )( ::SireCAS::IntegrationConstant const & ) ;
            assign_function_type assign_function_value( &::SireCAS::IntegrationConstant::operator= );
            
            IntegrationConstant_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        IntegrationConstant_exposer.def( bp::self == bp::self );
        { //::SireCAS::IntegrationConstant::stream
        
            typedef void ( ::SireCAS::IntegrationConstant::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireCAS::IntegrationConstant::stream );
            
            IntegrationConstant_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        IntegrationConstant_exposer.def( "__copy__", &__copy__);
        IntegrationConstant_exposer.def( "__deepcopy__", &__copy__);
        IntegrationConstant_exposer.def( "clone", &__copy__);
        IntegrationConstant_exposer.def( "__str__", &__str__< ::SireCAS::IntegrationConstant > );
        IntegrationConstant_exposer.def( "__repr__", &__str__< ::SireCAS::IntegrationConstant > );
    }

}