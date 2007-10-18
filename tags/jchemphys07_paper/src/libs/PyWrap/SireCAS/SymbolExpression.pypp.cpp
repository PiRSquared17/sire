// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "SymbolExpression.pypp.hpp"
#include "boost/python.hpp"
#include "sirecas_headers.h"

namespace bp = boost::python;

SireCAS::SymbolExpression __copy__(const SireCAS::SymbolExpression &other){ return SireCAS::SymbolExpression(other); }

const char* pvt_get_name(const SireCAS::SymbolExpression&){ return "SireCAS::SymbolExpression";}

void register_SymbolExpression_class(){

    bp::class_< SireCAS::SymbolExpression >( "SymbolExpression", bp::init< SireCAS::ExpressionBase const &, SireCAS::Expression const & >(( bp::arg("symbol"), bp::arg("expression") )) )    
        .def( 
            "expression"
            , (::SireCAS::Expression const & ( ::SireCAS::SymbolExpression::* )(  ) const)( &::SireCAS::SymbolExpression::expression )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "function"
            , (::SireCAS::Function const & ( ::SireCAS::SymbolExpression::* )(  ) const)( &::SireCAS::SymbolExpression::function )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "isFunction"
            , (bool ( ::SireCAS::SymbolExpression::* )(  ) const)( &::SireCAS::SymbolExpression::isFunction ) )    
        .def( 
            "symbol"
            , (::SireCAS::Symbol const & ( ::SireCAS::SymbolExpression::* )(  ) const)( &::SireCAS::SymbolExpression::symbol )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( "__copy__", &__copy__)    
        .def( "__str__", &pvt_get_name);

}
