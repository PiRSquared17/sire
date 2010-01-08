// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Constant.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/complex.h"

#include "Siren/stream.h"

#include "complexvalues.h"

#include "constant.h"

#include "expression.h"

#include "factor.h"

#include "identities.h"

#include "symbol.h"

#include "values.h"

#include "constant.h"

#include "Siren/logger.h"

SireCAS::Constant __copy__(const SireCAS::Constant &other){ return SireCAS::Constant(other); }

#include "Siren/str.hpp"

void register_Constant_class(){

    { //::SireCAS::Constant
        typedef bp::class_< SireCAS::Constant, bp::bases< SireCAS::CASNode, Siren::Object > > Constant_exposer_t;
        Constant_exposer_t Constant_exposer = Constant_exposer_t( "Constant" );
        bp::scope Constant_scope( Constant_exposer );
        Constant_exposer.def( bp::init< >() );
        Constant_exposer.def( bp::init< SireCAS::Constant const & >(( bp::arg("other") )) );
        { //::SireCAS::Constant::children
        
            typedef ::QList< SireCAS::Expression > ( ::SireCAS::Constant::*children_function_type )(  ) const;
            children_function_type children_function_value( &::SireCAS::Constant::children );
            
            Constant_exposer.def( 
                "children"
                , children_function_value );
        
        }
        { //::SireCAS::Constant::differentiate
        
            typedef ::SireCAS::Expression ( ::SireCAS::Constant::*differentiate_function_type )( ::SireCAS::Symbol const & ) const;
            differentiate_function_type differentiate_function_value( &::SireCAS::Constant::differentiate );
            
            Constant_exposer.def( 
                "differentiate"
                , differentiate_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireCAS::Constant::evaluate
        
            typedef double ( ::SireCAS::Constant::*evaluate_function_type )( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Constant::evaluate );
            
            Constant_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::Constant::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::Constant::*evaluate_function_type )( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Constant::evaluate );
            
            Constant_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::Constant::expand
        
            typedef ::QList< SireCAS::Factor > ( ::SireCAS::Constant::*expand_function_type )( ::SireCAS::Symbol const & ) const;
            expand_function_type expand_function_value( &::SireCAS::Constant::expand );
            
            Constant_exposer.def( 
                "expand"
                , expand_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireCAS::Constant::hashCode
        
            typedef ::uint ( ::SireCAS::Constant::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::SireCAS::Constant::hashCode );
            
            Constant_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::SireCAS::Constant::integrate
        
            typedef ::SireCAS::Expression ( ::SireCAS::Constant::*integrate_function_type )( ::SireCAS::Symbol const & ) const;
            integrate_function_type integrate_function_value( &::SireCAS::Constant::integrate );
            
            Constant_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("symbol") ) );
        
        }
        Constant_exposer.def( bp::self != bp::self );
        { //::SireCAS::Constant::operator=
        
            typedef ::SireCAS::Constant & ( ::SireCAS::Constant::*assign_function_type )( ::SireCAS::Constant const & ) ;
            assign_function_type assign_function_value( &::SireCAS::Constant::operator= );
            
            Constant_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Constant_exposer.def( bp::self == bp::self );
        { //::SireCAS::Constant::stream
        
            typedef void ( ::SireCAS::Constant::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireCAS::Constant::stream );
            
            Constant_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::SireCAS::Constant::substitute
        
            typedef ::SireCAS::Expression ( ::SireCAS::Constant::*substitute_function_type )( ::SireCAS::Identities const & ) const;
            substitute_function_type substitute_function_value( &::SireCAS::Constant::substitute );
            
            Constant_exposer.def( 
                "substitute"
                , substitute_function_value
                , ( bp::arg("identities") ) );
        
        }
        { //::SireCAS::Constant::symbols
        
            typedef ::QSet< SireCAS::Symbol > ( ::SireCAS::Constant::*symbols_function_type )(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::Constant::symbols );
            
            Constant_exposer.def( 
                "symbols"
                , symbols_function_value );
        
        }
        { //::SireCAS::Constant::toString
        
            typedef ::QString ( ::SireCAS::Constant::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireCAS::Constant::toString );
            
            Constant_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        Constant_exposer.def( "__copy__", &__copy__);
        Constant_exposer.def( "__deepcopy__", &__copy__);
        Constant_exposer.def( "clone", &__copy__);
        Constant_exposer.def( "__str__", &__str__< ::SireCAS::Constant > );
        Constant_exposer.def( "__repr__", &__str__< ::SireCAS::Constant > );
    }

}
