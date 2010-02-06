// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ArcTanh.pypp.hpp"

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

SireCAS::ArcTanh __copy__(const SireCAS::ArcTanh &other){ return SireCAS::ArcTanh(other); }

#include "Siren/str.hpp"

void register_ArcTanh_class(){

    { //::SireCAS::ArcTanh
        typedef bp::class_< SireCAS::ArcTanh, bp::bases< SireCAS::SingleFunc, SireCAS::CASNode, Siren::Object > > ArcTanh_exposer_t;
        ArcTanh_exposer_t ArcTanh_exposer = ArcTanh_exposer_t( "ArcTanh" );
        bp::scope ArcTanh_scope( ArcTanh_exposer );
        ArcTanh_exposer.def( bp::init< >() );
        ArcTanh_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) );
        ArcTanh_exposer.def( bp::init< SireCAS::ArcTanh const & >(( bp::arg("other") )) );
        { //::SireCAS::ArcTanh::evaluate
        
            typedef double ( ::SireCAS::ArcTanh::*evaluate_function_type )( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ArcTanh::evaluate );
            
            ArcTanh_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::ArcTanh::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::ArcTanh::*evaluate_function_type )( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ArcTanh::evaluate );
            
            ArcTanh_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        ArcTanh_exposer.def( bp::self != bp::self );
        { //::SireCAS::ArcTanh::operator=
        
            typedef ::SireCAS::ArcTanh & ( ::SireCAS::ArcTanh::*assign_function_type )( ::SireCAS::ArcTanh const & ) ;
            assign_function_type assign_function_value( &::SireCAS::ArcTanh::operator= );
            
            ArcTanh_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        ArcTanh_exposer.def( bp::self == bp::self );
        { //::SireCAS::ArcTanh::stream
        
            typedef void ( ::SireCAS::ArcTanh::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireCAS::ArcTanh::stream );
            
            ArcTanh_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        ArcTanh_exposer.def( "__copy__", &__copy__);
        ArcTanh_exposer.def( "__deepcopy__", &__copy__);
        ArcTanh_exposer.def( "clone", &__copy__);
        ArcTanh_exposer.def( "__str__", &__str__< ::SireCAS::ArcTanh > );
        ArcTanh_exposer.def( "__repr__", &__str__< ::SireCAS::ArcTanh > );
    }

}