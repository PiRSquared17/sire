// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Camera.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireUnits/units.h"

#include "camera.h"

#include <QGLWidget>

#include "camera.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Camera_class(){

    { //::Spier::Camera
        typedef bp::class_< Spier::Camera, bp::bases< SireBase::Property >, boost::noncopyable > Camera_exposer_t;
        Camera_exposer_t Camera_exposer = Camera_exposer_t( "Camera", bp::no_init );
        bp::scope Camera_scope( Camera_exposer );
        { //::Spier::Camera::depthCue
        
            typedef bool ( ::Spier::Camera::*depthCue_function_type )(  ) const;
            depthCue_function_type depthCue_function_value( &::Spier::Camera::depthCue );
            
            Camera_exposer.def( 
                "depthCue"
                , depthCue_function_value );
        
        }
        { //::Spier::Camera::differentView
        
            typedef bool ( ::Spier::Camera::*differentView_function_type )( ::Spier::Camera const & ) const;
            differentView_function_type differentView_function_value( &::Spier::Camera::differentView );
            
            Camera_exposer.def( 
                "differentView"
                , differentView_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::Spier::Camera::location
        
            typedef ::SireMaths::Vector const & ( ::Spier::Camera::*location_function_type )(  ) const;
            location_function_type location_function_value( &::Spier::Camera::location );
            
            Camera_exposer.def( 
                "location"
                , location_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Spier::Camera::lookAt
        
            typedef void ( ::Spier::Camera::*lookAt_function_type )( ::SireMaths::Vector const & ) ;
            lookAt_function_type lookAt_function_value( &::Spier::Camera::lookAt );
            
            Camera_exposer.def( 
                "lookAt"
                , lookAt_function_value
                , ( bp::arg("lookat") ) );
        
        }
        { //::Spier::Camera::null
        
            typedef ::Spier::NullCamera const & ( *null_function_type )(  );
            null_function_type null_function_value( &::Spier::Camera::null );
            
            Camera_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Spier::Camera::range
        
            typedef float ( ::Spier::Camera::*range_function_type )(  ) const;
            range_function_type range_function_value( &::Spier::Camera::range );
            
            Camera_exposer.def( 
                "range"
                , range_function_value );
        
        }
        { //::Spier::Camera::reset
        
            typedef void ( ::Spier::Camera::*reset_function_type )(  ) ;
            reset_function_type reset_function_value( &::Spier::Camera::reset );
            
            Camera_exposer.def( 
                "reset"
                , reset_function_value );
        
        }
        { //::Spier::Camera::rotate
        
            typedef void ( ::Spier::Camera::*rotate_function_type )( ::SireUnits::Dimension::Angle const &,::SireUnits::Dimension::Angle const & ) ;
            rotate_function_type rotate_function_value( &::Spier::Camera::rotate );
            
            Camera_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("dx"), bp::arg("dy") ) );
        
        }
        { //::Spier::Camera::sideVector
        
            typedef ::SireMaths::Vector const & ( ::Spier::Camera::*sideVector_function_type )(  ) const;
            sideVector_function_type sideVector_function_value( &::Spier::Camera::sideVector );
            
            Camera_exposer.def( 
                "sideVector"
                , sideVector_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Spier::Camera::spin
        
            typedef void ( ::Spier::Camera::*spin_function_type )( ::SireUnits::Dimension::Angle const & ) ;
            spin_function_type spin_function_value( &::Spier::Camera::spin );
            
            Camera_exposer.def( 
                "spin"
                , spin_function_value
                , ( bp::arg("delta") ) );
        
        }
        { //::Spier::Camera::translate
        
            typedef void ( ::Spier::Camera::*translate_function_type )( double,double ) ;
            translate_function_type translate_function_value( &::Spier::Camera::translate );
            
            Camera_exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("dx"), bp::arg("dy") ) );
        
        }
        { //::Spier::Camera::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Spier::Camera::typeName );
            
            Camera_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::Spier::Camera::upVector
        
            typedef ::SireMaths::Vector const & ( ::Spier::Camera::*upVector_function_type )(  ) const;
            upVector_function_type upVector_function_value( &::Spier::Camera::upVector );
            
            Camera_exposer.def( 
                "upVector"
                , upVector_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Spier::Camera::viewAngle
        
            typedef ::SireUnits::Dimension::Angle const & ( ::Spier::Camera::*viewAngle_function_type )(  ) const;
            viewAngle_function_type viewAngle_function_value( &::Spier::Camera::viewAngle );
            
            Camera_exposer.def( 
                "viewAngle"
                , viewAngle_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Spier::Camera::viewVector
        
            typedef ::SireMaths::Vector const & ( ::Spier::Camera::*viewVector_function_type )(  ) const;
            viewVector_function_type viewVector_function_value( &::Spier::Camera::viewVector );
            
            Camera_exposer.def( 
                "viewVector"
                , viewVector_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Spier::Camera::zoom
        
            typedef void ( ::Spier::Camera::*zoom_function_type )( double ) ;
            zoom_function_type zoom_function_value( &::Spier::Camera::zoom );
            
            Camera_exposer.def( 
                "zoom"
                , zoom_function_value
                , ( bp::arg("delta") ) );
        
        }
        Camera_exposer.staticmethod( "null" );
        Camera_exposer.staticmethod( "typeName" );
        Camera_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Spier::Camera >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Camera_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Spier::Camera >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Camera_exposer.def( "__str__", &__str__< ::Spier::Camera > );
        Camera_exposer.def( "__repr__", &__str__< ::Spier::Camera > );
    }

}
