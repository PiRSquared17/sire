// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Plane.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "plane.h"

#include "plane.h"

SireMaths::Plane __copy__(const SireMaths::Plane &other){ return SireMaths::Plane(other); }

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMaths::Plane&){ return "SireMaths::Plane";}

void register_Plane_class(){

    { //::SireMaths::Plane
        typedef bp::class_< SireMaths::Plane > Plane_exposer_t;
        Plane_exposer_t Plane_exposer = Plane_exposer_t( "Plane" );
        bp::scope Plane_scope( Plane_exposer );
        Plane_exposer.def( bp::init< >() );
        Plane_exposer.def( bp::init< SireMaths::Vector const &, double const & >(( bp::arg("normal"), bp::arg("distance") )) );
        Plane_exposer.def( bp::init< double const &, double const &, double const &, double const & >(( bp::arg("a"), bp::arg("b"), bp::arg("c"), bp::arg("d") )) );
        Plane_exposer.def( bp::init< SireMaths::Vector const &, SireMaths::Vector const & >(( bp::arg("normal"), bp::arg("contains_point") )) );
        Plane_exposer.def( bp::init< SireMaths::Plane const & >(( bp::arg("other") )) );
        { //::SireMaths::Plane::a
        
            typedef double ( ::SireMaths::Plane::*a_function_type )(  ) const;
            a_function_type a_function_value( &::SireMaths::Plane::a );
            
            Plane_exposer.def( 
                "a"
                , a_function_value );
        
        }
        { //::SireMaths::Plane::b
        
            typedef double ( ::SireMaths::Plane::*b_function_type )(  ) const;
            b_function_type b_function_value( &::SireMaths::Plane::b );
            
            Plane_exposer.def( 
                "b"
                , b_function_value );
        
        }
        { //::SireMaths::Plane::c
        
            typedef double ( ::SireMaths::Plane::*c_function_type )(  ) const;
            c_function_type c_function_value( &::SireMaths::Plane::c );
            
            Plane_exposer.def( 
                "c"
                , c_function_value );
        
        }
        { //::SireMaths::Plane::d
        
            typedef double ( ::SireMaths::Plane::*d_function_type )(  ) const;
            d_function_type d_function_value( &::SireMaths::Plane::d );
            
            Plane_exposer.def( 
                "d"
                , d_function_value );
        
        }
        { //::SireMaths::Plane::distance
        
            typedef double ( ::SireMaths::Plane::*distance_function_type )( ::SireMaths::Vector const & ) const;
            distance_function_type distance_function_value( &::SireMaths::Plane::distance );
            
            Plane_exposer.def( 
                "distance"
                , distance_function_value
                , ( bp::arg("point") ) );
        
        }
        { //::SireMaths::Plane::distanceFromOrigin
        
            typedef double const & ( ::SireMaths::Plane::*distanceFromOrigin_function_type )(  ) const;
            distanceFromOrigin_function_type distanceFromOrigin_function_value( &::SireMaths::Plane::distanceFromOrigin );
            
            Plane_exposer.def( 
                "distanceFromOrigin"
                , distanceFromOrigin_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMaths::Plane::normal
        
            typedef ::SireMaths::Vector const & ( ::SireMaths::Plane::*normal_function_type )(  ) const;
            normal_function_type normal_function_value( &::SireMaths::Plane::normal );
            
            Plane_exposer.def( 
                "normal"
                , normal_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMaths::Plane::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::Plane::typeName );
            
            Plane_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMaths::Plane::what
        
            typedef char const * ( ::SireMaths::Plane::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMaths::Plane::what );
            
            Plane_exposer.def( 
                "what"
                , what_function_value );
        
        }
        Plane_exposer.staticmethod( "typeName" );
        Plane_exposer.def( "__copy__", &__copy__);
        Plane_exposer.def( "__deepcopy__", &__copy__);
        Plane_exposer.def( "clone", &__copy__);
        Plane_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::Plane >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Plane_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::Plane >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Plane_exposer.def( "__str__", &pvt_get_name);
        Plane_exposer.def( "__repr__", &pvt_get_name);
    }

}
