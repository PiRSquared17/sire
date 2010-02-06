// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Sech.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/complex.h"

#include "Siren/stream.h"

#include "complexvalues.h"

#include "exp.h"

#include "expression.h"

#include "hyperbolicfuncs.h"

#include "identities.h"

#include "trigfuncs.h"

#include "hyperbolicfuncs.h"

#include "Siren/logger.h"

SireCAS::Sech __copy__(const SireCAS::Sech &other){ return SireCAS::Sech(other); }

#include "Siren/str.hpp"

void register_Sech_class(){

    { //::SireCAS::Sech
        typedef bp::class_< SireCAS::Sech, bp::bases< SireCAS::SingleFunc, SireCAS::CASNode, Siren::Object > > Sech_exposer_t;
        Sech_exposer_t Sech_exposer = Sech_exposer_t( "Sech" );
        bp::scope Sech_scope( Sech_exposer );
        Sech_exposer.def( bp::init< >() );
        Sech_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) );
        Sech_exposer.def( bp::init< SireCAS::Sech const & >(( bp::arg("other") )) );
        { //::SireCAS::Sech::evaluate
        
            typedef double ( ::SireCAS::Sech::*evaluate_function_type )( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Sech::evaluate );
            
            Sech_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::Sech::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::Sech::*evaluate_function_type )( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Sech::evaluate );
            
            Sech_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        Sech_exposer.def( bp::self != bp::self );
        { //::SireCAS::Sech::operator=
        
            typedef ::SireCAS::Sech & ( ::SireCAS::Sech::*assign_function_type )( ::SireCAS::Sech const & ) ;
            assign_function_type assign_function_value( &::SireCAS::Sech::operator= );
            
            Sech_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Sech_exposer.def( bp::self == bp::self );
        { //::SireCAS::Sech::stream
        
            typedef void ( ::SireCAS::Sech::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireCAS::Sech::stream );
            
            Sech_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        Sech_exposer.def( "__copy__", &__copy__);
        Sech_exposer.def( "__deepcopy__", &__copy__);
        Sech_exposer.def( "clone", &__copy__);
        Sech_exposer.def( "__str__", &__str__< ::SireCAS::Sech > );
        Sech_exposer.def( "__repr__", &__str__< ::SireCAS::Sech > );
    }

}