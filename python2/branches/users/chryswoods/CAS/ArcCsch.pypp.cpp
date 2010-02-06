// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ArcCsch.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/complex.h"

#include "SireMaths/errors.h"

#include "Siren/stream.h"

#include "complexvalues.h"

#include "exp.h"

#include "expression.h"

#include "hyperbolicfuncs.h"

#include "identities.h"

#include "invhyperbolicfuncs.h"

#include "invtrigfuncs.h"

#include "trigfuncs.h"

#include <boost/math/special_functions/asinh.hpp>

#include "invhyperbolicfuncs.h"

#include "Siren/logger.h"

SireCAS::ArcCsch __copy__(const SireCAS::ArcCsch &other){ return SireCAS::ArcCsch(other); }

#include "Siren/str.hpp"

void register_ArcCsch_class(){

    { //::SireCAS::ArcCsch
        typedef bp::class_< SireCAS::ArcCsch, bp::bases< SireCAS::SingleFunc, SireCAS::CASNode, Siren::Object > > ArcCsch_exposer_t;
        ArcCsch_exposer_t ArcCsch_exposer = ArcCsch_exposer_t( "ArcCsch" );
        bp::scope ArcCsch_scope( ArcCsch_exposer );
        ArcCsch_exposer.def( bp::init< >() );
        ArcCsch_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) );
        ArcCsch_exposer.def( bp::init< SireCAS::ArcCsch const & >(( bp::arg("other") )) );
        { //::SireCAS::ArcCsch::evaluate
        
            typedef double ( ::SireCAS::ArcCsch::*evaluate_function_type )( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ArcCsch::evaluate );
            
            ArcCsch_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::ArcCsch::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::ArcCsch::*evaluate_function_type )( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ArcCsch::evaluate );
            
            ArcCsch_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        ArcCsch_exposer.def( bp::self != bp::self );
        { //::SireCAS::ArcCsch::operator=
        
            typedef ::SireCAS::ArcCsch & ( ::SireCAS::ArcCsch::*assign_function_type )( ::SireCAS::ArcCsch const & ) ;
            assign_function_type assign_function_value( &::SireCAS::ArcCsch::operator= );
            
            ArcCsch_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        ArcCsch_exposer.def( bp::self == bp::self );
        { //::SireCAS::ArcCsch::stream
        
            typedef void ( ::SireCAS::ArcCsch::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireCAS::ArcCsch::stream );
            
            ArcCsch_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        ArcCsch_exposer.def( "__copy__", &__copy__);
        ArcCsch_exposer.def( "__deepcopy__", &__copy__);
        ArcCsch_exposer.def( "clone", &__copy__);
        ArcCsch_exposer.def( "__str__", &__str__< ::SireCAS::ArcCsch > );
        ArcCsch_exposer.def( "__repr__", &__str__< ::SireCAS::ArcCsch > );
    }

}