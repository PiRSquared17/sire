// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "WeakHandle.pypp.hpp"

namespace bp = boost::python;

#include "Siren/errors.h"

#include "class.h"

#include "handle.h"

#include "hanref.h"

#include "logger.h"

#include "waitcondition.h"

#include "handle.h"

Siren::WeakHandle __copy__(const Siren::WeakHandle &other){ return Siren::WeakHandle(other); }

const char* pvt_get_name(const Siren::WeakHandle&){ return "Siren::WeakHandle";}

void register_WeakHandle_class(){

    { //::Siren::WeakHandle
        typedef bp::class_< Siren::WeakHandle > WeakHandle_exposer_t;
        WeakHandle_exposer_t WeakHandle_exposer = WeakHandle_exposer_t( "WeakHandle" );
        bp::scope WeakHandle_scope( WeakHandle_exposer );
        WeakHandle_exposer.def( bp::init< >() );
        WeakHandle_exposer.def( bp::init< Siren::Handle const & >(( bp::arg("handle") )) );
        WeakHandle_exposer.def( bp::init< Siren::WeakHandle const & >(( bp::arg("other") )) );
        { //::Siren::WeakHandle::expired
        
            typedef bool ( ::Siren::WeakHandle::*expired_function_type )(  ) const;
            expired_function_type expired_function_value( &::Siren::WeakHandle::expired );
            
            WeakHandle_exposer.def( 
                "expired"
                , expired_function_value );
        
        }
        { //::Siren::WeakHandle::lock
        
            typedef ::Siren::HanRef ( ::Siren::WeakHandle::*lock_function_type )(  ) const;
            lock_function_type lock_function_value( &::Siren::WeakHandle::lock );
            
            WeakHandle_exposer.def( 
                "lock"
                , lock_function_value );
        
        }
        WeakHandle_exposer.def( bp::self != bp::self );
        { //::Siren::WeakHandle::operator=
        
            typedef ::Siren::WeakHandle & ( ::Siren::WeakHandle::*assign_function_type )( ::Siren::WeakHandle const & ) ;
            assign_function_type assign_function_value( &::Siren::WeakHandle::operator= );
            
            WeakHandle_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        WeakHandle_exposer.def( bp::self == bp::self );
        WeakHandle_exposer.def( "__copy__", &__copy__);
        WeakHandle_exposer.def( "__deepcopy__", &__copy__);
        WeakHandle_exposer.def( "clone", &__copy__);
        WeakHandle_exposer.def( "__str__", &pvt_get_name);
        WeakHandle_exposer.def( "__repr__", &pvt_get_name);
    }

}