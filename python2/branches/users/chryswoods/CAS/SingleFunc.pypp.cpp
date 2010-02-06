// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SingleFunc.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/errors.h"

#include "SireMaths/complex.h"

#include "Siren/stream.h"

#include "complexvalues.h"

#include "factor.h"

#include "identities.h"

#include "singlefunc.h"

#include "values.h"

#include "singlefunc.h"

#include "Siren/logger.h"

#include "Siren/str.hpp"

void register_SingleFunc_class(){

    { //::SireCAS::SingleFunc
        typedef bp::class_< SireCAS::SingleFunc, bp::bases< SireCAS::CASNode, Siren::Object >, boost::noncopyable > SingleFunc_exposer_t;
        SingleFunc_exposer_t SingleFunc_exposer = SingleFunc_exposer_t( "SingleFunc", bp::no_init );
        bp::scope SingleFunc_scope( SingleFunc_exposer );
        { //::SireCAS::SingleFunc::argument
        
            typedef ::SireCAS::Expression const & ( ::SireCAS::SingleFunc::*argument_function_type )(  ) const;
            argument_function_type argument_function_value( &::SireCAS::SingleFunc::argument );
            
            SingleFunc_exposer.def( 
                "argument"
                , argument_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireCAS::SingleFunc::children
        
            typedef ::QList< SireCAS::Expression > ( ::SireCAS::SingleFunc::*children_function_type )(  ) const;
            children_function_type children_function_value( &::SireCAS::SingleFunc::children );
            
            SingleFunc_exposer.def( 
                "children"
                , children_function_value );
        
        }
        { //::SireCAS::SingleFunc::conjugate
        
            typedef ::SireCAS::Expression ( ::SireCAS::SingleFunc::*conjugate_function_type )(  ) const;
            conjugate_function_type conjugate_function_value( &::SireCAS::SingleFunc::conjugate );
            
            SingleFunc_exposer.def( 
                "conjugate"
                , conjugate_function_value );
        
        }
        { //::SireCAS::SingleFunc::differentiate
        
            typedef ::SireCAS::Expression ( ::SireCAS::SingleFunc::*differentiate_function_type )( ::SireCAS::Symbol const & ) const;
            differentiate_function_type differentiate_function_value( &::SireCAS::SingleFunc::differentiate );
            
            SingleFunc_exposer.def( 
                "differentiate"
                , differentiate_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireCAS::SingleFunc::expand
        
            typedef ::QList< SireCAS::Factor > ( ::SireCAS::SingleFunc::*expand_function_type )( ::SireCAS::Symbol const & ) const;
            expand_function_type expand_function_value( &::SireCAS::SingleFunc::expand );
            
            SingleFunc_exposer.def( 
                "expand"
                , expand_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireCAS::SingleFunc::hashCode
        
            typedef ::uint ( ::SireCAS::SingleFunc::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::SireCAS::SingleFunc::hashCode );
            
            SingleFunc_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::SireCAS::SingleFunc::integrate
        
            typedef ::SireCAS::Expression ( ::SireCAS::SingleFunc::*integrate_function_type )( ::SireCAS::Symbol const & ) const;
            integrate_function_type integrate_function_value( &::SireCAS::SingleFunc::integrate );
            
            SingleFunc_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireCAS::SingleFunc::isComplex
        
            typedef bool ( ::SireCAS::SingleFunc::*isComplex_function_type )(  ) const;
            isComplex_function_type isComplex_function_value( &::SireCAS::SingleFunc::isComplex );
            
            SingleFunc_exposer.def( 
                "isComplex"
                , isComplex_function_value );
        
        }
        { //::SireCAS::SingleFunc::isCompound
        
            typedef bool ( ::SireCAS::SingleFunc::*isCompound_function_type )(  ) const;
            isCompound_function_type isCompound_function_value( &::SireCAS::SingleFunc::isCompound );
            
            SingleFunc_exposer.def( 
                "isCompound"
                , isCompound_function_value );
        
        }
        { //::SireCAS::SingleFunc::isConstant
        
            typedef bool ( ::SireCAS::SingleFunc::*isConstant_function_type )(  ) const;
            isConstant_function_type isConstant_function_value( &::SireCAS::SingleFunc::isConstant );
            
            SingleFunc_exposer.def( 
                "isConstant"
                , isConstant_function_value );
        
        }
        { //::SireCAS::SingleFunc::isFunction
        
            typedef bool ( ::SireCAS::SingleFunc::*isFunction_function_type )( ::SireCAS::Symbol const & ) const;
            isFunction_function_type isFunction_function_value( &::SireCAS::SingleFunc::isFunction );
            
            SingleFunc_exposer.def( 
                "isFunction"
                , isFunction_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireCAS::SingleFunc::stream
        
            typedef void ( ::SireCAS::SingleFunc::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireCAS::SingleFunc::stream );
            
            SingleFunc_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::SireCAS::SingleFunc::substitute
        
            typedef ::SireCAS::Expression ( ::SireCAS::SingleFunc::*substitute_function_type )( ::SireCAS::Identities const & ) const;
            substitute_function_type substitute_function_value( &::SireCAS::SingleFunc::substitute );
            
            SingleFunc_exposer.def( 
                "substitute"
                , substitute_function_value
                , ( bp::arg("identities") ) );
        
        }
        { //::SireCAS::SingleFunc::symbols
        
            typedef ::QSet< SireCAS::Symbol > ( ::SireCAS::SingleFunc::*symbols_function_type )(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::SingleFunc::symbols );
            
            SingleFunc_exposer.def( 
                "symbols"
                , symbols_function_value );
        
        }
        { //::SireCAS::SingleFunc::toString
        
            typedef ::QString ( ::SireCAS::SingleFunc::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireCAS::SingleFunc::toString );
            
            SingleFunc_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireCAS::SingleFunc::typeName
        
            typedef ::QString ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::SingleFunc::typeName );
            
            SingleFunc_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireCAS::SingleFunc::x
        
            typedef ::SireCAS::Expression const & ( ::SireCAS::SingleFunc::*x_function_type )(  ) const;
            x_function_type x_function_value( &::SireCAS::SingleFunc::x );
            
            SingleFunc_exposer.def( 
                "x"
                , x_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        SingleFunc_exposer.staticmethod( "typeName" );
        SingleFunc_exposer.def( "__str__", &__str__< ::SireCAS::SingleFunc > );
        SingleFunc_exposer.def( "__repr__", &__str__< ::SireCAS::SingleFunc > );
    }

}