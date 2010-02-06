// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Sin.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/complex.h"

#include "Siren/stream.h"

#include "complexvalues.h"

#include "exp.h"

#include "expression.h"

#include "identities.h"

#include "trigfuncs.h"

#include <cmath>

#include "trigfuncs.h"

#include "Siren/logger.h"

SireCAS::Sin __copy__(const SireCAS::Sin &other){ return SireCAS::Sin(other); }

#include "Siren/str.hpp"

void register_Sin_class(){

    { //::SireCAS::Sin
        typedef bp::class_< SireCAS::Sin, bp::bases< SireCAS::SingleFunc, SireCAS::CASNode, Siren::Object > > Sin_exposer_t;
        Sin_exposer_t Sin_exposer = Sin_exposer_t( "Sin" );
        bp::scope Sin_scope( Sin_exposer );
        Sin_exposer.def( bp::init< >() );
        Sin_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) );
        Sin_exposer.def( bp::init< SireCAS::Sin const & >(( bp::arg("other") )) );
        { //::SireCAS::Sin::evaluate
        
            typedef double ( ::SireCAS::Sin::*evaluate_function_type )( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Sin::evaluate );
            
            Sin_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::Sin::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::Sin::*evaluate_function_type )( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Sin::evaluate );
            
            Sin_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        Sin_exposer.def( bp::self != bp::self );
        { //::SireCAS::Sin::operator=
        
            typedef ::SireCAS::Sin & ( ::SireCAS::Sin::*assign_function_type )( ::SireCAS::Sin const & ) ;
            assign_function_type assign_function_value( &::SireCAS::Sin::operator= );
            
            Sin_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Sin_exposer.def( bp::self == bp::self );
        { //::SireCAS::Sin::stream
        
            typedef void ( ::SireCAS::Sin::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireCAS::Sin::stream );
            
            Sin_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        Sin_exposer.def( "__copy__", &__copy__);
        Sin_exposer.def( "__deepcopy__", &__copy__);
        Sin_exposer.def( "clone", &__copy__);
        Sin_exposer.def( "__str__", &__str__< ::SireCAS::Sin > );
        Sin_exposer.def( "__repr__", &__str__< ::SireCAS::Sin > );
    }

}