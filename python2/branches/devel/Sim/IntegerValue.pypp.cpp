// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "IntegerValue.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireSim/values.h"

#include "values.h"

#include <QDomDocument>

#include <QDomElement>

#include "values.h"

SireSim::IntegerValue __copy__(const SireSim::IntegerValue &other){ return SireSim::IntegerValue(other); }

#include "Helpers/str.hpp"

void register_IntegerValue_class(){

    { //::SireSim::IntegerValue
        typedef bp::class_< SireSim::IntegerValue, bp::bases< SireSim::Value > > IntegerValue_exposer_t;
        IntegerValue_exposer_t IntegerValue_exposer = IntegerValue_exposer_t( "IntegerValue", bp::init< >() );
        bp::scope IntegerValue_scope( IntegerValue_exposer );
        IntegerValue_exposer.def( bp::init< qint64 >(( bp::arg("val") )) );
        IntegerValue_exposer.def( bp::init< QString >(( bp::arg("val") )) );
        IntegerValue_exposer.def( bp::init< QDomElement >(( bp::arg("elem") )) );
        IntegerValue_exposer.def( bp::init< SireSim::IntegerValue const & >(( bp::arg("other") )) );
        IntegerValue_exposer.def( bp::self != bp::self );
        { //::SireSim::IntegerValue::operator=
        
            typedef ::SireSim::IntegerValue & ( ::SireSim::IntegerValue::*assign_function_type )( ::SireSim::IntegerValue const & ) ;
            assign_function_type assign_function_value( &::SireSim::IntegerValue::operator= );
            
            IntegerValue_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        IntegerValue_exposer.def( bp::self == bp::self );
        { //::SireSim::IntegerValue::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSim::IntegerValue::typeName );
            
            IntegerValue_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireSim::IntegerValue::value
        
            typedef ::qint64 ( ::SireSim::IntegerValue::*value_function_type )(  ) const;
            value_function_type value_function_value( &::SireSim::IntegerValue::value );
            
            IntegerValue_exposer.def( 
                "value"
                , value_function_value );
        
        }
        { //::SireSim::IntegerValue::valueWithMaximum
        
            typedef ::SireSim::IntegerValue ( *valueWithMaximum_function_type )( ::qint64,::qint64 );
            valueWithMaximum_function_type valueWithMaximum_function_value( &::SireSim::IntegerValue::valueWithMaximum );
            
            IntegerValue_exposer.def( 
                "valueWithMaximum"
                , valueWithMaximum_function_value
                , ( bp::arg("val"), bp::arg("maximum") ) );
        
        }
        { //::SireSim::IntegerValue::valueWithMinimum
        
            typedef ::SireSim::IntegerValue ( *valueWithMinimum_function_type )( ::qint64,::qint64 );
            valueWithMinimum_function_type valueWithMinimum_function_value( &::SireSim::IntegerValue::valueWithMinimum );
            
            IntegerValue_exposer.def( 
                "valueWithMinimum"
                , valueWithMinimum_function_value
                , ( bp::arg("val"), bp::arg("minimum") ) );
        
        }
        { //::SireSim::IntegerValue::valueWithRange
        
            typedef ::SireSim::IntegerValue ( *valueWithRange_function_type )( ::qint64,::qint64,::qint64 );
            valueWithRange_function_type valueWithRange_function_value( &::SireSim::IntegerValue::valueWithRange );
            
            IntegerValue_exposer.def( 
                "valueWithRange"
                , valueWithRange_function_value
                , ( bp::arg("val"), bp::arg("minimum"), bp::arg("maximum") ) );
        
        }
        { //::SireSim::IntegerValue::what
        
            typedef char const * ( ::SireSim::IntegerValue::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireSim::IntegerValue::what );
            
            IntegerValue_exposer.def( 
                "what"
                , what_function_value );
        
        }
        IntegerValue_exposer.staticmethod( "typeName" );
        IntegerValue_exposer.staticmethod( "valueWithMaximum" );
        IntegerValue_exposer.staticmethod( "valueWithMinimum" );
        IntegerValue_exposer.staticmethod( "valueWithRange" );
        IntegerValue_exposer.def( "__copy__", &__copy__);
        IntegerValue_exposer.def( "__deepcopy__", &__copy__);
        IntegerValue_exposer.def( "clone", &__copy__);
        IntegerValue_exposer.def( "__str__", &__str__< ::SireSim::IntegerValue > );
        IntegerValue_exposer.def( "__repr__", &__str__< ::SireSim::IntegerValue > );
    }

}
