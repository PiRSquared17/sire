// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Value.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireSim/value.h"

#include "SireSim/values.h"

#include "value.h"

#include <QDebug>

#include <QDomElement>

#include "value.h"

#include "Helpers/str.hpp"

void register_Value_class(){

    { //::SireSim::Value
        typedef bp::class_< SireSim::Value, boost::noncopyable > Value_exposer_t;
        Value_exposer_t Value_exposer = Value_exposer_t( "Value", bp::no_init );
        bp::scope Value_scope( Value_exposer );
        { //::SireSim::Value::getNestedValue
        
            typedef ::SireSim::ValuePtr ( ::SireSim::Value::*getNestedValue_function_type )( ::QString ) const;
            getNestedValue_function_type getNestedValue_function_value( &::SireSim::Value::getNestedValue );
            
            Value_exposer.def( 
                "getNestedValue"
                , getNestedValue_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireSim::Value::setNestedValue
        
            typedef ::SireSim::ValuePtr ( ::SireSim::Value::*setNestedValue_function_type )( ::QString,::QString ) const;
            setNestedValue_function_type setNestedValue_function_value( &::SireSim::Value::setNestedValue );
            
            Value_exposer.def( 
                "setNestedValue"
                , setNestedValue_function_value
                , ( bp::arg("key"), bp::arg("value") ) );
        
        }
        { //::SireSim::Value::setNestedValue
        
            typedef ::SireSim::ValuePtr ( ::SireSim::Value::*setNestedValue_function_type )( ::QString,double ) const;
            setNestedValue_function_type setNestedValue_function_value( &::SireSim::Value::setNestedValue );
            
            Value_exposer.def( 
                "setNestedValue"
                , setNestedValue_function_value
                , ( bp::arg("key"), bp::arg("value") ) );
        
        }
        { //::SireSim::Value::setNestedValue
        
            typedef ::SireSim::ValuePtr ( ::SireSim::Value::*setNestedValue_function_type )( ::QString,::qint64 ) const;
            setNestedValue_function_type setNestedValue_function_value( &::SireSim::Value::setNestedValue );
            
            Value_exposer.def( 
                "setNestedValue"
                , setNestedValue_function_value
                , ( bp::arg("key"), bp::arg("value") ) );
        
        }
        { //::SireSim::Value::toConfig
        
            typedef ::QString ( ::SireSim::Value::*toConfig_function_type )( bool ) const;
            toConfig_function_type toConfig_function_value( &::SireSim::Value::toConfig );
            
            Value_exposer.def( 
                "toConfig"
                , toConfig_function_value
                , ( bp::arg("include_help")=(bool)(false) ) );
        
        }
        { //::SireSim::Value::toString
        
            typedef ::QString ( ::SireSim::Value::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireSim::Value::toString );
            
            Value_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireSim::Value::toXML
        
            typedef ::QString ( ::SireSim::Value::*toXML_function_type )(  ) const;
            toXML_function_type toXML_function_value( &::SireSim::Value::toXML );
            
            Value_exposer.def( 
                "toXML"
                , toXML_function_value );
        
        }
        { //::SireSim::Value::what
        
            typedef char const * ( ::SireSim::Value::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireSim::Value::what );
            
            Value_exposer.def( 
                "what"
                , what_function_value );
        
        }
        bp::register_ptr_to_python< boost::shared_ptr< SireSim::Value > >();
        Value_exposer.def( "__str__", &__str__< ::SireSim::Value > );
        Value_exposer.def( "__repr__", &__str__< ::SireSim::Value > );
    }

}
