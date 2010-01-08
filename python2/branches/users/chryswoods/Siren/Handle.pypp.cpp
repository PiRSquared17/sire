// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Handle.pypp.hpp"

namespace bp = boost::python;

#include "Siren/errors.h"

#include "class.h"

#include "handle.h"

#include "hanref.h"

#include "logger.h"

#include "handle.h"

#include "Siren/logger.h"

#include "Siren/str.hpp"

void register_Handle_class(){

    { //::Siren::Handle
        typedef bp::class_< Siren::Handle, boost::noncopyable > Handle_exposer_t;
        Handle_exposer_t Handle_exposer = Handle_exposer_t( "Handle", bp::no_init );
        bp::scope Handle_scope( Handle_exposer );
        { //::Siren::Handle::assertNotNull
        
            typedef void ( ::Siren::Handle::*assertNotNull_function_type )(  ) const;
            assertNotNull_function_type assertNotNull_function_value( &::Siren::Handle::assertNotNull );
            
            Handle_exposer.def( 
                "assertNotNull"
                , assertNotNull_function_value );
        
        }
        { //::Siren::Handle::copy
        
            typedef void ( ::Siren::Handle::*copy_function_type )( ::Siren::Handle const & ) ;
            copy_function_type copy_function_value( &::Siren::Handle::copy );
            
            Handle_exposer.def( 
                "copy"
                , copy_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::Siren::Handle::equals
        
            typedef bool ( ::Siren::Handle::*equals_function_type )( ::Siren::Handle const & ) const;
            equals_function_type equals_function_value( &::Siren::Handle::equals );
            
            Handle_exposer.def( 
                "equals"
                , equals_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::Siren::Handle::getClass
        
            typedef ::Siren::Class const & ( ::Siren::Handle::*getClass_function_type )(  ) const;
            getClass_function_type getClass_function_value( &::Siren::Handle::getClass );
            
            Handle_exposer.def( 
                "getClass"
                , getClass_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Siren::Handle::hashCode
        
            typedef ::uint ( ::Siren::Handle::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::Siren::Handle::hashCode );
            
            Handle_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::Siren::Handle::isNull
        
            typedef bool ( ::Siren::Handle::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::Siren::Handle::isNull );
            
            Handle_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::Siren::Handle::lock
        
            typedef void ( ::Siren::Handle::*lock_function_type )(  ) ;
            lock_function_type lock_function_value( &::Siren::Handle::lock );
            
            Handle_exposer.def( 
                "lock"
                , lock_function_value );
        
        }
        Handle_exposer.def( bp::self != bp::self );
        { //::Siren::Handle::operator=
        
            typedef ::Siren::Handle & ( ::Siren::Handle::*assign_function_type )( ::Siren::Handle const & ) ;
            assign_function_type assign_function_value( &::Siren::Handle::operator= );
            
            Handle_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Handle_exposer.def( bp::self == bp::self );
        { //::Siren::Handle::test
        
            typedef bool ( ::Siren::Handle::*test_function_type )(  ) const;
            test_function_type test_function_value( &::Siren::Handle::test );
            
            Handle_exposer.def( 
                "test"
                , test_function_value );
        
        }
        { //::Siren::Handle::test
        
            typedef bool ( ::Siren::Handle::*test_function_type )( ::Siren::Logger & ) const;
            test_function_type test_function_value( &::Siren::Handle::test );
            
            Handle_exposer.def( 
                "test"
                , test_function_value
                , ( bp::arg("logger") ) );
        
        }
        { //::Siren::Handle::toString
        
            typedef ::QString ( ::Siren::Handle::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::Siren::Handle::toString );
            
            Handle_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::Siren::Handle::tryLock
        
            typedef bool ( ::Siren::Handle::*tryLock_function_type )(  ) ;
            tryLock_function_type tryLock_function_value( &::Siren::Handle::tryLock );
            
            Handle_exposer.def( 
                "tryLock"
                , tryLock_function_value );
        
        }
        { //::Siren::Handle::tryLock
        
            typedef bool ( ::Siren::Handle::*tryLock_function_type )( int ) ;
            tryLock_function_type tryLock_function_value( &::Siren::Handle::tryLock );
            
            Handle_exposer.def( 
                "tryLock"
                , tryLock_function_value
                , ( bp::arg("ms") ) );
        
        }
        { //::Siren::Handle::typeName
        
            typedef ::QString ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Siren::Handle::typeName );
            
            Handle_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::Siren::Handle::unlock
        
            typedef void ( ::Siren::Handle::*unlock_function_type )(  ) ;
            unlock_function_type unlock_function_value( &::Siren::Handle::unlock );
            
            Handle_exposer.def( 
                "unlock"
                , unlock_function_value );
        
        }
        { //::Siren::Handle::what
        
            typedef ::QString ( ::Siren::Handle::*what_function_type )(  ) const;
            what_function_type what_function_value( &::Siren::Handle::what );
            
            Handle_exposer.def( 
                "what"
                , what_function_value );
        
        }
        Handle_exposer.staticmethod( "typeName" );
        Handle_exposer.def( "__str__", &__str__< ::Siren::Handle > );
        Handle_exposer.def( "__repr__", &__str__< ::Siren::Handle > );
    }

}
