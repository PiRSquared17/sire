// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PlaneBase.pypp.hpp"

namespace bp = boost::python;

#include "Siren/stream.h"

#include "plane.h"

#include "plane.h"

#include "Siren/logger.h"

Siren::Primitive<SireMaths::Plane> __copy__(const Siren::Primitive<SireMaths::Plane> &other){ return Siren::Primitive<SireMaths::Plane>(other); }

const char* pvt_get_name(const Siren::Primitive<SireMaths::Plane>&){ return "Siren::Primitive<SireMaths::Plane>";}

void register_PlaneBase_class(){

    { //::Siren::Primitive< SireMaths::Plane >
        typedef bp::class_< Siren::Primitive< SireMaths::Plane > > PlaneBase_exposer_t;
        PlaneBase_exposer_t PlaneBase_exposer = PlaneBase_exposer_t( "PlaneBase" );
        bp::scope PlaneBase_scope( PlaneBase_exposer );
        PlaneBase_exposer.def( bp::init< >() );
        { //::Siren::Primitive< SireMaths::Plane >::getClass
        
            typedef Siren::Primitive< SireMaths::Plane > exported_class_t;
            typedef ::Siren::Class ( ::Siren::Primitive< SireMaths::Plane >::*getClass_function_type )(  ) const;
            getClass_function_type getClass_function_value( &::Siren::Primitive< SireMaths::Plane >::getClass );
            
            PlaneBase_exposer.def( 
                "getClass"
                , getClass_function_value );
        
        }
        { //::Siren::Primitive< SireMaths::Plane >::hashCode
        
            typedef Siren::Primitive< SireMaths::Plane > exported_class_t;
            typedef ::uint ( ::Siren::Primitive< SireMaths::Plane >::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::Siren::Primitive< SireMaths::Plane >::hashCode );
            
            PlaneBase_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::Siren::Primitive< SireMaths::Plane >::load
        
            typedef Siren::Primitive< SireMaths::Plane > exported_class_t;
            typedef void ( ::Siren::Primitive< SireMaths::Plane >::*load_function_type )( ::Siren::Stream & ) ;
            load_function_type load_function_value( &::Siren::Primitive< SireMaths::Plane >::load );
            
            PlaneBase_exposer.def( 
                "load"
                , load_function_value
                , ( bp::arg("s") ) );
        
        }
        PlaneBase_exposer.def( bp::self != bp::other< SireMaths::Plane >() );
        { //::Siren::Primitive< SireMaths::Plane >::save
        
            typedef Siren::Primitive< SireMaths::Plane > exported_class_t;
            typedef void ( ::Siren::Primitive< SireMaths::Plane >::*save_function_type )( ::Siren::Stream & ) const;
            save_function_type save_function_value( &::Siren::Primitive< SireMaths::Plane >::save );
            
            PlaneBase_exposer.def( 
                "save"
                , save_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::Siren::Primitive< SireMaths::Plane >::stream
        
            typedef Siren::Primitive< SireMaths::Plane > exported_class_t;
            typedef void ( ::Siren::Primitive< SireMaths::Plane >::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::Siren::Primitive< SireMaths::Plane >::stream );
            
            PlaneBase_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::Siren::Primitive< SireMaths::Plane >::test
        
            typedef Siren::Primitive< SireMaths::Plane > exported_class_t;
            typedef bool ( ::Siren::Primitive< SireMaths::Plane >::*test_function_type )(  ) const;
            test_function_type test_function_value( &::Siren::Primitive< SireMaths::Plane >::test );
            
            PlaneBase_exposer.def( 
                "test"
                , test_function_value );
        
        }
        { //::Siren::Primitive< SireMaths::Plane >::test
        
            typedef Siren::Primitive< SireMaths::Plane > exported_class_t;
            typedef bool ( ::Siren::Primitive< SireMaths::Plane >::*test_function_type )( ::Siren::Logger & ) const;
            test_function_type test_function_value( &::Siren::Primitive< SireMaths::Plane >::test );
            
            PlaneBase_exposer.def( 
                "test"
                , test_function_value
                , ( bp::arg("logger") ) );
        
        }
        { //::Siren::Primitive< SireMaths::Plane >::typeName
        
            typedef Siren::Primitive< SireMaths::Plane > exported_class_t;
            typedef ::QString ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Siren::Primitive< SireMaths::Plane >::typeName );
            
            PlaneBase_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::Siren::Primitive< SireMaths::Plane >::what
        
            typedef Siren::Primitive< SireMaths::Plane > exported_class_t;
            typedef ::QString ( ::Siren::Primitive< SireMaths::Plane >::*what_function_type )(  ) const;
            what_function_type what_function_value( &::Siren::Primitive< SireMaths::Plane >::what );
            
            PlaneBase_exposer.def( 
                "what"
                , what_function_value );
        
        }
        PlaneBase_exposer.staticmethod( "typeName" );
        PlaneBase_exposer.def( "__copy__", &__copy__);
        PlaneBase_exposer.def( "__deepcopy__", &__copy__);
        PlaneBase_exposer.def( "clone", &__copy__);
        PlaneBase_exposer.def( "__str__", &pvt_get_name);
        PlaneBase_exposer.def( "__repr__", &pvt_get_name);
    }

}
