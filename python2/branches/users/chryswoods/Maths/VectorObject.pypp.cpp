// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "VectorObject.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/quickcopy.hpp"

#include "SireUnits/units.h"

#include "Siren/errors.h"

#include "Siren/logger.h"

#include "Siren/stream.h"

#include "Siren/tester.h"

#include "matrix.h"

#include "quaternion.h"

#include "vector.h"

#include <QRegExp>

#include <QString>

#include <cmath>

#include "vector.h"

#include "Siren/logger.h"

Siren::PrimitiveObject<SireMaths::Vector> __copy__(const Siren::PrimitiveObject<SireMaths::Vector> &other){ return Siren::PrimitiveObject<SireMaths::Vector>(other); }

#include "Siren/str.hpp"

void register_VectorObject_class(){

    { //::Siren::PrimitiveObject< SireMaths::Vector >
        typedef bp::class_< Siren::PrimitiveObject< SireMaths::Vector >, bp::bases< Siren::Object > > VectorObject_exposer_t;
        VectorObject_exposer_t VectorObject_exposer = VectorObject_exposer_t( "VectorObject" );
        bp::scope VectorObject_scope( VectorObject_exposer );
        VectorObject_exposer.def( bp::init< >() );
        VectorObject_exposer.def( bp::init< SireMaths::Vector const & >(( bp::arg("primitive") )) );
        VectorObject_exposer.def( bp::init< Siren::PrimitiveObject< SireMaths::Vector > const & >(( bp::arg("other") )) );
        { //::Siren::PrimitiveObject< SireMaths::Vector >::hashCode
        
            typedef Siren::PrimitiveObject< SireMaths::Vector > exported_class_t;
            typedef ::uint ( ::Siren::PrimitiveObject< SireMaths::Vector >::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::Siren::PrimitiveObject< SireMaths::Vector >::hashCode );
            
            VectorObject_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        VectorObject_exposer.def( bp::self != bp::self );
        VectorObject_exposer.def( bp::self != bp::other< SireMaths::Vector >() );
        { //::Siren::PrimitiveObject< SireMaths::Vector >::operator=
        
            typedef Siren::PrimitiveObject< SireMaths::Vector > exported_class_t;
            typedef ::Siren::PrimitiveObject< SireMaths::Vector > & ( ::Siren::PrimitiveObject< SireMaths::Vector >::*assign_function_type )( ::Siren::PrimitiveObject< SireMaths::Vector > const & ) ;
            assign_function_type assign_function_value( &::Siren::PrimitiveObject< SireMaths::Vector >::operator= );
            
            VectorObject_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::Siren::PrimitiveObject< SireMaths::Vector >::operator=
        
            typedef Siren::PrimitiveObject< SireMaths::Vector > exported_class_t;
            typedef ::Siren::PrimitiveObject< SireMaths::Vector > & ( ::Siren::PrimitiveObject< SireMaths::Vector >::*assign_function_type )( ::SireMaths::Vector const & ) ;
            assign_function_type assign_function_value( &::Siren::PrimitiveObject< SireMaths::Vector >::operator= );
            
            VectorObject_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("primitive") )
                , bp::return_self< >() );
        
        }
        VectorObject_exposer.def( bp::self == bp::self );
        VectorObject_exposer.def( bp::self == bp::other< SireMaths::Vector >() );
        { //::Siren::PrimitiveObject< SireMaths::Vector >::stream
        
            typedef Siren::PrimitiveObject< SireMaths::Vector > exported_class_t;
            typedef void ( ::Siren::PrimitiveObject< SireMaths::Vector >::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::Siren::PrimitiveObject< SireMaths::Vector >::stream );
            
            VectorObject_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::Siren::PrimitiveObject< SireMaths::Vector >::test
        
            typedef Siren::PrimitiveObject< SireMaths::Vector > exported_class_t;
            typedef bool ( ::Siren::PrimitiveObject< SireMaths::Vector >::*test_function_type )( ::Siren::Logger & ) const;
            test_function_type test_function_value( &::Siren::PrimitiveObject< SireMaths::Vector >::test );
            
            VectorObject_exposer.def( 
                "test"
                , test_function_value
                , ( bp::arg("logger") ) );
        
        }
        { //::Siren::PrimitiveObject< SireMaths::Vector >::toString
        
            typedef Siren::PrimitiveObject< SireMaths::Vector > exported_class_t;
            typedef ::QString ( ::Siren::PrimitiveObject< SireMaths::Vector >::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::Siren::PrimitiveObject< SireMaths::Vector >::toString );
            
            VectorObject_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        VectorObject_exposer.def( "__copy__", &__copy__);
        VectorObject_exposer.def( "__deepcopy__", &__copy__);
        VectorObject_exposer.def( "clone", &__copy__);
        VectorObject_exposer.def( "__str__", &__str__< ::Siren::PrimitiveObject<SireMaths::Vector> > );
        VectorObject_exposer.def( "__repr__", &__str__< ::Siren::PrimitiveObject<SireMaths::Vector> > );
    }

}
