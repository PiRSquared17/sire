// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ComplexBase.pypp.hpp"

namespace bp = boost::python;

#include "complex.h"

#include "rational.h"

#include <boost/math/complex.hpp>

#include "complex.h"

#include "Siren/logger.h"

Siren::Primitive<SireMaths::Complex> __copy__(const Siren::Primitive<SireMaths::Complex> &other){ return Siren::Primitive<SireMaths::Complex>(other); }

const char* pvt_get_name(const Siren::Primitive<SireMaths::Complex>&){ return "Siren::Primitive<SireMaths::Complex>";}

void register_ComplexBase_class(){

    { //::Siren::Primitive< SireMaths::Complex >
        typedef bp::class_< Siren::Primitive< SireMaths::Complex > > ComplexBase_exposer_t;
        ComplexBase_exposer_t ComplexBase_exposer = ComplexBase_exposer_t( "ComplexBase" );
        bp::scope ComplexBase_scope( ComplexBase_exposer );
        ComplexBase_exposer.def( bp::init< >() );
        { //::Siren::Primitive< SireMaths::Complex >::getClass
        
            typedef Siren::Primitive< SireMaths::Complex > exported_class_t;
            typedef ::Siren::Class ( ::Siren::Primitive< SireMaths::Complex >::*getClass_function_type )(  ) const;
            getClass_function_type getClass_function_value( &::Siren::Primitive< SireMaths::Complex >::getClass );
            
            ComplexBase_exposer.def( 
                "getClass"
                , getClass_function_value );
        
        }
        { //::Siren::Primitive< SireMaths::Complex >::hashCode
        
            typedef Siren::Primitive< SireMaths::Complex > exported_class_t;
            typedef ::uint ( ::Siren::Primitive< SireMaths::Complex >::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::Siren::Primitive< SireMaths::Complex >::hashCode );
            
            ComplexBase_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::Siren::Primitive< SireMaths::Complex >::load
        
            typedef Siren::Primitive< SireMaths::Complex > exported_class_t;
            typedef void ( ::Siren::Primitive< SireMaths::Complex >::*load_function_type )( ::Siren::Stream & ) ;
            load_function_type load_function_value( &::Siren::Primitive< SireMaths::Complex >::load );
            
            ComplexBase_exposer.def( 
                "load"
                , load_function_value
                , ( bp::arg("s") ) );
        
        }
        ComplexBase_exposer.def( bp::self != bp::other< SireMaths::Complex >() );
        { //::Siren::Primitive< SireMaths::Complex >::save
        
            typedef Siren::Primitive< SireMaths::Complex > exported_class_t;
            typedef void ( ::Siren::Primitive< SireMaths::Complex >::*save_function_type )( ::Siren::Stream & ) const;
            save_function_type save_function_value( &::Siren::Primitive< SireMaths::Complex >::save );
            
            ComplexBase_exposer.def( 
                "save"
                , save_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::Siren::Primitive< SireMaths::Complex >::stream
        
            typedef Siren::Primitive< SireMaths::Complex > exported_class_t;
            typedef void ( ::Siren::Primitive< SireMaths::Complex >::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::Siren::Primitive< SireMaths::Complex >::stream );
            
            ComplexBase_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::Siren::Primitive< SireMaths::Complex >::test
        
            typedef Siren::Primitive< SireMaths::Complex > exported_class_t;
            typedef bool ( ::Siren::Primitive< SireMaths::Complex >::*test_function_type )(  ) const;
            test_function_type test_function_value( &::Siren::Primitive< SireMaths::Complex >::test );
            
            ComplexBase_exposer.def( 
                "test"
                , test_function_value );
        
        }
        { //::Siren::Primitive< SireMaths::Complex >::test
        
            typedef Siren::Primitive< SireMaths::Complex > exported_class_t;
            typedef bool ( ::Siren::Primitive< SireMaths::Complex >::*test_function_type )( ::Siren::Logger & ) const;
            test_function_type test_function_value( &::Siren::Primitive< SireMaths::Complex >::test );
            
            ComplexBase_exposer.def( 
                "test"
                , test_function_value
                , ( bp::arg("logger") ) );
        
        }
        { //::Siren::Primitive< SireMaths::Complex >::typeName
        
            typedef Siren::Primitive< SireMaths::Complex > exported_class_t;
            typedef ::QString ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Siren::Primitive< SireMaths::Complex >::typeName );
            
            ComplexBase_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::Siren::Primitive< SireMaths::Complex >::what
        
            typedef Siren::Primitive< SireMaths::Complex > exported_class_t;
            typedef ::QString ( ::Siren::Primitive< SireMaths::Complex >::*what_function_type )(  ) const;
            what_function_type what_function_value( &::Siren::Primitive< SireMaths::Complex >::what );
            
            ComplexBase_exposer.def( 
                "what"
                , what_function_value );
        
        }
        ComplexBase_exposer.staticmethod( "typeName" );
        ComplexBase_exposer.def( "__copy__", &__copy__);
        ComplexBase_exposer.def( "__deepcopy__", &__copy__);
        ComplexBase_exposer.def( "clone", &__copy__);
        ComplexBase_exposer.def( "__str__", &pvt_get_name);
        ComplexBase_exposer.def( "__repr__", &pvt_get_name);
    }

}
