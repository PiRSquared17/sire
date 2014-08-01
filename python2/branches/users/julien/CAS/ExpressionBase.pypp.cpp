// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ExpressionBase.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "complexvalues.h"

#include "constant.h"

#include "exbase.h"

#include "expressionbase.h"

#include "expressions.h"

#include "functions.h"

#include "identities.h"

#include "symbols.h"

#include "expressionbase.h"

SireCAS::ExpressionBase __copy__(const SireCAS::ExpressionBase &other){ return SireCAS::ExpressionBase(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ExpressionBase_class(){

    { //::SireCAS::ExpressionBase
        typedef bp::class_< SireCAS::ExpressionBase > ExpressionBase_exposer_t;
        ExpressionBase_exposer_t ExpressionBase_exposer = ExpressionBase_exposer_t( "ExpressionBase", bp::init< >() );
        bp::scope ExpressionBase_scope( ExpressionBase_exposer );
        ExpressionBase_exposer.def( bp::init< SireCAS::ExBase const & >(( bp::arg("ex") )) );
        ExpressionBase_exposer.def( bp::init< SireCAS::ExpressionBase const & >(( bp::arg("other") )) );
        { //::SireCAS::ExpressionBase::children
        
            typedef ::SireCAS::Expressions ( ::SireCAS::ExpressionBase::*children_function_type )(  ) const;
            children_function_type children_function_value( &::SireCAS::ExpressionBase::children );
            
            ExpressionBase_exposer.def( 
                "children"
                , children_function_value );
        
        }
        { //::SireCAS::ExpressionBase::conjugate
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExpressionBase::*conjugate_function_type )(  ) const;
            conjugate_function_type conjugate_function_value( &::SireCAS::ExpressionBase::conjugate );
            
            ExpressionBase_exposer.def( 
                "conjugate"
                , conjugate_function_value );
        
        }
        { //::SireCAS::ExpressionBase::differentiate
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExpressionBase::*differentiate_function_type )( ::SireCAS::Symbol const & ) const;
            differentiate_function_type differentiate_function_value( &::SireCAS::ExpressionBase::differentiate );
            
            ExpressionBase_exposer.def( 
                "differentiate"
                , differentiate_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireCAS::ExpressionBase::evaluate
        
            typedef double ( ::SireCAS::ExpressionBase::*evaluate_function_type )( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ExpressionBase::evaluate );
            
            ExpressionBase_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::ExpressionBase::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::ExpressionBase::*evaluate_function_type )( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ExpressionBase::evaluate );
            
            ExpressionBase_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::ExpressionBase::expand
        
            typedef ::QList< SireCAS::Factor > ( ::SireCAS::ExpressionBase::*expand_function_type )( ::SireCAS::Symbol const & ) const;
            expand_function_type expand_function_value( &::SireCAS::ExpressionBase::expand );
            
            ExpressionBase_exposer.def( 
                "expand"
                , expand_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireCAS::ExpressionBase::functions
        
            typedef ::SireCAS::Functions ( ::SireCAS::ExpressionBase::*functions_function_type )(  ) const;
            functions_function_type functions_function_value( &::SireCAS::ExpressionBase::functions );
            
            ExpressionBase_exposer.def( 
                "functions"
                , functions_function_value );
        
        }
        { //::SireCAS::ExpressionBase::hash
        
            typedef ::uint ( ::SireCAS::ExpressionBase::*hash_function_type )(  ) const;
            hash_function_type hash_function_value( &::SireCAS::ExpressionBase::hash );
            
            ExpressionBase_exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireCAS::ExpressionBase::integrate
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExpressionBase::*integrate_function_type )( ::SireCAS::Symbol const & ) const;
            integrate_function_type integrate_function_value( &::SireCAS::ExpressionBase::integrate );
            
            ExpressionBase_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireCAS::ExpressionBase::isComplex
        
            typedef bool ( ::SireCAS::ExpressionBase::*isComplex_function_type )(  ) const;
            isComplex_function_type isComplex_function_value( &::SireCAS::ExpressionBase::isComplex );
            
            ExpressionBase_exposer.def( 
                "isComplex"
                , isComplex_function_value );
        
        }
        { //::SireCAS::ExpressionBase::isCompound
        
            typedef bool ( ::SireCAS::ExpressionBase::*isCompound_function_type )(  ) const;
            isCompound_function_type isCompound_function_value( &::SireCAS::ExpressionBase::isCompound );
            
            ExpressionBase_exposer.def( 
                "isCompound"
                , isCompound_function_value );
        
        }
        { //::SireCAS::ExpressionBase::isConstant
        
            typedef bool ( ::SireCAS::ExpressionBase::*isConstant_function_type )(  ) const;
            isConstant_function_type isConstant_function_value( &::SireCAS::ExpressionBase::isConstant );
            
            ExpressionBase_exposer.def( 
                "isConstant"
                , isConstant_function_value );
        
        }
        { //::SireCAS::ExpressionBase::isFunction
        
            typedef bool ( ::SireCAS::ExpressionBase::*isFunction_function_type )( ::SireCAS::Symbol const & ) const;
            isFunction_function_type isFunction_function_value( &::SireCAS::ExpressionBase::isFunction );
            
            ExpressionBase_exposer.def( 
                "isFunction"
                , isFunction_function_value
                , ( bp::arg("arg0") ) );
        
        }
        ExpressionBase_exposer.def( bp::self != bp::self );
        { //::SireCAS::ExpressionBase::operator=
        
            typedef ::SireCAS::ExpressionBase & ( ::SireCAS::ExpressionBase::*assign_function_type )( ::SireCAS::ExpressionBase const & ) ;
            assign_function_type assign_function_value( &::SireCAS::ExpressionBase::operator= );
            
            ExpressionBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireCAS::ExpressionBase::operator=
        
            typedef ::SireCAS::ExpressionBase & ( ::SireCAS::ExpressionBase::*assign_function_type )( ::SireCAS::ExBase const & ) ;
            assign_function_type assign_function_value( &::SireCAS::ExpressionBase::operator= );
            
            ExpressionBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        ExpressionBase_exposer.def( bp::self == bp::self );
        { //::SireCAS::ExpressionBase::series
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExpressionBase::*series_function_type )( ::SireCAS::Symbol const &,int ) const;
            series_function_type series_function_value( &::SireCAS::ExpressionBase::series );
            
            ExpressionBase_exposer.def( 
                "series"
                , series_function_value
                , ( bp::arg("symbol"), bp::arg("n") ) );
        
        }
        { //::SireCAS::ExpressionBase::simplify
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExpressionBase::*simplify_function_type )( int ) const;
            simplify_function_type simplify_function_value( &::SireCAS::ExpressionBase::simplify );
            
            ExpressionBase_exposer.def( 
                "simplify"
                , simplify_function_value
                , ( bp::arg("options")=(int)(0) ) );
        
        }
        { //::SireCAS::ExpressionBase::substitute
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExpressionBase::*substitute_function_type )( ::SireCAS::Identities const & ) const;
            substitute_function_type substitute_function_value( &::SireCAS::ExpressionBase::substitute );
            
            ExpressionBase_exposer.def( 
                "substitute"
                , substitute_function_value
                , ( bp::arg("identities") ) );
        
        }
        { //::SireCAS::ExpressionBase::symbols
        
            typedef ::SireCAS::Symbols ( ::SireCAS::ExpressionBase::*symbols_function_type )(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::ExpressionBase::symbols );
            
            ExpressionBase_exposer.def( 
                "symbols"
                , symbols_function_value );
        
        }
        { //::SireCAS::ExpressionBase::toOpenMMString
        
            typedef ::QString ( ::SireCAS::ExpressionBase::*toOpenMMString_function_type )(  ) const;
            toOpenMMString_function_type toOpenMMString_function_value( &::SireCAS::ExpressionBase::toOpenMMString );
            
            ExpressionBase_exposer.def( 
                "toOpenMMString"
                , toOpenMMString_function_value );
        
        }
        { //::SireCAS::ExpressionBase::toString
        
            typedef ::QString ( ::SireCAS::ExpressionBase::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireCAS::ExpressionBase::toString );
            
            ExpressionBase_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireCAS::ExpressionBase::what
        
            typedef char const * ( ::SireCAS::ExpressionBase::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireCAS::ExpressionBase::what );
            
            ExpressionBase_exposer.def( 
                "what"
                , what_function_value );
        
        }
        ExpressionBase_exposer.def( "__copy__", &__copy__);
        ExpressionBase_exposer.def( "__deepcopy__", &__copy__);
        ExpressionBase_exposer.def( "clone", &__copy__);
        ExpressionBase_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::ExpressionBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ExpressionBase_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::ExpressionBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ExpressionBase_exposer.def( "__str__", &__str__< ::SireCAS::ExpressionBase > );
        ExpressionBase_exposer.def( "__repr__", &__str__< ::SireCAS::ExpressionBase > );
    }

}
