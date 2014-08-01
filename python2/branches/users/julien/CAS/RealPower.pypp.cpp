// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "RealPower.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "expression.h"

#include "expressions.h"

#include "identities.h"

#include "integrationconstant.h"

#include "powerconstant.h"

#include "symbol.h"

#include "symbols.h"

#include "values.h"

#include <QDebug>

#include "powerconstant.h"

SireCAS::RealPower __copy__(const SireCAS::RealPower &other){ return SireCAS::RealPower(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_RealPower_class(){

    { //::SireCAS::RealPower
        typedef bp::class_< SireCAS::RealPower, bp::bases< SireCAS::PowerFunction, SireCAS::ExBase > > RealPower_exposer_t;
        RealPower_exposer_t RealPower_exposer = RealPower_exposer_t( "RealPower", bp::init< >() );
        bp::scope RealPower_scope( RealPower_exposer );
        RealPower_exposer.def( bp::init< SireCAS::Expression const &, double >(( bp::arg("expression"), bp::arg("power") )) );
        RealPower_exposer.def( bp::init< SireCAS::RealPower const & >(( bp::arg("other") )) );
        { //::SireCAS::RealPower::evaluate
        
            typedef double ( ::SireCAS::RealPower::*evaluate_function_type )( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::RealPower::evaluate );
            
            RealPower_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::RealPower::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::RealPower::*evaluate_function_type )( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::RealPower::evaluate );
            
            RealPower_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::RealPower::hash
        
            typedef ::uint ( ::SireCAS::RealPower::*hash_function_type )(  ) const;
            hash_function_type hash_function_value( &::SireCAS::RealPower::hash );
            
            RealPower_exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        RealPower_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::RealPower::power
        
            typedef ::SireCAS::Expression ( ::SireCAS::RealPower::*power_function_type )(  ) const;
            power_function_type power_function_value( &::SireCAS::RealPower::power );
            
            RealPower_exposer.def( 
                "power"
                , power_function_value );
        
        }
        { //::SireCAS::RealPower::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::RealPower::typeName );
            
            RealPower_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireCAS::RealPower::what
        
            typedef char const * ( ::SireCAS::RealPower::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireCAS::RealPower::what );
            
            RealPower_exposer.def( 
                "what"
                , what_function_value );
        
        }
        RealPower_exposer.staticmethod( "typeName" );
        RealPower_exposer.def( "__copy__", &__copy__);
        RealPower_exposer.def( "__deepcopy__", &__copy__);
        RealPower_exposer.def( "clone", &__copy__);
        RealPower_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::RealPower >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RealPower_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::RealPower >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RealPower_exposer.def( "__str__", &__str__< ::SireCAS::RealPower > );
        RealPower_exposer.def( "__repr__", &__str__< ::SireCAS::RealPower > );
    }

}
