// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "HanRef.pypp.hpp"

namespace bp = boost::python;

#include "Siren/errors.h"

#include "class.h"

#include "handle.h"

#include "hanref.h"

#include "logger.h"

#include "stream.h"

#include "tester.h"

#include "hanref.h"

#include "Siren/logger.h"

Siren::HanRef __copy__(const Siren::HanRef &other){ return Siren::HanRef(other); }

const char* pvt_get_name(const Siren::HanRef&){ return "Siren::HanRef";}

void register_HanRef_class(){

    { //::Siren::HanRef
        typedef bp::class_< Siren::HanRef > HanRef_exposer_t;
        HanRef_exposer_t HanRef_exposer = HanRef_exposer_t( "HanRef" );
        bp::scope HanRef_scope( HanRef_exposer );
        HanRef_exposer.def( bp::init< >() );
        HanRef_exposer.def( bp::init< Siren::Handle * >(( bp::arg("handle") )) );
        HanRef_exposer.def( bp::init< Siren::Handle const & >(( bp::arg("handle") )) );
        HanRef_exposer.def( bp::init< Siren::HanRef const & >(( bp::arg("other") )) );
        { //::Siren::HanRef::copy
        
            typedef void ( ::Siren::HanRef::*copy_function_type )( ::Siren::Handle const & ) ;
            copy_function_type copy_function_value( &::Siren::HanRef::copy );
            
            HanRef_exposer.def( 
                "copy"
                , copy_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::Siren::HanRef::equals
        
            typedef bool ( ::Siren::HanRef::*equals_function_type )( ::Siren::Handle const & ) const;
            equals_function_type equals_function_value( &::Siren::HanRef::equals );
            
            HanRef_exposer.def( 
                "equals"
                , equals_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::Siren::HanRef::getClass
        
            typedef ::Siren::Class const & ( ::Siren::HanRef::*getClass_function_type )(  ) const;
            getClass_function_type getClass_function_value( &::Siren::HanRef::getClass );
            
            HanRef_exposer.def( 
                "getClass"
                , getClass_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Siren::HanRef::hashCode
        
            typedef ::uint ( ::Siren::HanRef::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::Siren::HanRef::hashCode );
            
            HanRef_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::Siren::HanRef::isNull
        
            typedef bool ( ::Siren::HanRef::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::Siren::HanRef::isNull );
            
            HanRef_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::Siren::HanRef::lock
        
            typedef void ( ::Siren::HanRef::*lock_function_type )(  ) ;
            lock_function_type lock_function_value( &::Siren::HanRef::lock );
            
            HanRef_exposer.def( 
                "lock"
                , lock_function_value );
        
        }
        HanRef_exposer.def( bp::self != bp::self );
        { //::Siren::HanRef::operator=
        
            typedef ::Siren::HanRef & ( ::Siren::HanRef::*assign_function_type )( ::Siren::HanRef const & ) ;
            assign_function_type assign_function_value( &::Siren::HanRef::operator= );
            
            HanRef_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        HanRef_exposer.def( bp::self == bp::self );
        { //::Siren::HanRef::test
        
            typedef bool ( ::Siren::HanRef::*test_function_type )(  ) const;
            test_function_type test_function_value( &::Siren::HanRef::test );
            
            HanRef_exposer.def( 
                "test"
                , test_function_value );
        
        }
        { //::Siren::HanRef::test
        
            typedef bool ( ::Siren::HanRef::*test_function_type )( ::Siren::Logger & ) const;
            test_function_type test_function_value( &::Siren::HanRef::test );
            
            HanRef_exposer.def( 
                "test"
                , test_function_value
                , ( bp::arg("logger") ) );
        
        }
        { //::Siren::HanRef::tryLock
        
            typedef bool ( ::Siren::HanRef::*tryLock_function_type )(  ) ;
            tryLock_function_type tryLock_function_value( &::Siren::HanRef::tryLock );
            
            HanRef_exposer.def( 
                "tryLock"
                , tryLock_function_value );
        
        }
        { //::Siren::HanRef::tryLock
        
            typedef bool ( ::Siren::HanRef::*tryLock_function_type )( int ) ;
            tryLock_function_type tryLock_function_value( &::Siren::HanRef::tryLock );
            
            HanRef_exposer.def( 
                "tryLock"
                , tryLock_function_value
                , ( bp::arg("ms") ) );
        
        }
        { //::Siren::HanRef::unlock
        
            typedef void ( ::Siren::HanRef::*unlock_function_type )(  ) ;
            unlock_function_type unlock_function_value( &::Siren::HanRef::unlock );
            
            HanRef_exposer.def( 
                "unlock"
                , unlock_function_value );
        
        }
        { //::Siren::HanRef::what
        
            typedef ::QString ( ::Siren::HanRef::*what_function_type )(  ) const;
            what_function_type what_function_value( &::Siren::HanRef::what );
            
            HanRef_exposer.def( 
                "what"
                , what_function_value );
        
        }
        HanRef_exposer.def( "__copy__", &__copy__);
        HanRef_exposer.def( "__deepcopy__", &__copy__);
        HanRef_exposer.def( "clone", &__copy__);
        HanRef_exposer.def( "__str__", &pvt_get_name);
        HanRef_exposer.def( "__repr__", &pvt_get_name);
    }

}
