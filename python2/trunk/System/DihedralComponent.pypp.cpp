// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "DihedralComponent.pypp.hpp"

namespace bp = boost::python;

#include "SireID/index.h"

#include "SireMaths/torsion.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/space.h"

#include "delta.h"

#include "dihedralcomponent.h"

#include "dihedralcomponent.h"

SireSystem::DihedralComponent __copy__(const SireSystem::DihedralComponent &other){ return SireSystem::DihedralComponent(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_DihedralComponent_class(){

    { //::SireSystem::DihedralComponent
        typedef bp::class_< SireSystem::DihedralComponent, bp::bases< SireSystem::GeometryComponent, SireSystem::Constraint, SireBase::Property > > DihedralComponent_exposer_t;
        DihedralComponent_exposer_t DihedralComponent_exposer = DihedralComponent_exposer_t( "DihedralComponent", bp::init< >() );
        bp::scope DihedralComponent_scope( DihedralComponent_exposer );
        DihedralComponent_exposer.def( bp::init< SireCAS::Symbol const &, SireFF::PointRef const &, SireFF::PointRef const &, SireFF::PointRef const &, SireFF::PointRef const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("constrained_symbol"), bp::arg("point0"), bp::arg("point1"), bp::arg("point2"), bp::arg("point3"), bp::arg("map")=SireBase::PropertyMap() )) );
        DihedralComponent_exposer.def( bp::init< SireCAS::Symbol const &, SireFF::PointRef const &, SireFF::PointRef const &, SireFF::PointRef const &, SireFF::PointRef const &, SireCAS::Expression const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("constrained_symbol"), bp::arg("point0"), bp::arg("point1"), bp::arg("point2"), bp::arg("point3"), bp::arg("geometry_expression"), bp::arg("map")=SireBase::PropertyMap() )) );
        DihedralComponent_exposer.def( bp::init< SireSystem::DihedralComponent const & >(( bp::arg("other") )) );
        { //::SireSystem::DihedralComponent::nPoints
        
            typedef int ( ::SireSystem::DihedralComponent::*nPoints_function_type )(  ) const;
            nPoints_function_type nPoints_function_value( &::SireSystem::DihedralComponent::nPoints );
            
            DihedralComponent_exposer.def( 
                "nPoints"
                , nPoints_function_value );
        
        }
        DihedralComponent_exposer.def( bp::self != bp::self );
        { //::SireSystem::DihedralComponent::operator=
        
            typedef ::SireSystem::DihedralComponent & ( ::SireSystem::DihedralComponent::*assign_function_type )( ::SireSystem::DihedralComponent const & ) ;
            assign_function_type assign_function_value( &::SireSystem::DihedralComponent::operator= );
            
            DihedralComponent_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        DihedralComponent_exposer.def( bp::self == bp::self );
        { //::SireSystem::DihedralComponent::phi
        
            typedef ::SireCAS::Symbol const & ( *phi_function_type )(  );
            phi_function_type phi_function_value( &::SireSystem::DihedralComponent::phi );
            
            DihedralComponent_exposer.def( 
                "phi"
                , phi_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::point
        
            typedef ::SireFF::Point const & ( ::SireSystem::DihedralComponent::*point_function_type )( int ) const;
            point_function_type point_function_value( &::SireSystem::DihedralComponent::point );
            
            DihedralComponent_exposer.def( 
                "point"
                , point_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::point0
        
            typedef ::SireFF::Point const & ( ::SireSystem::DihedralComponent::*point0_function_type )(  ) const;
            point0_function_type point0_function_value( &::SireSystem::DihedralComponent::point0 );
            
            DihedralComponent_exposer.def( 
                "point0"
                , point0_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::point1
        
            typedef ::SireFF::Point const & ( ::SireSystem::DihedralComponent::*point1_function_type )(  ) const;
            point1_function_type point1_function_value( &::SireSystem::DihedralComponent::point1 );
            
            DihedralComponent_exposer.def( 
                "point1"
                , point1_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::point2
        
            typedef ::SireFF::Point const & ( ::SireSystem::DihedralComponent::*point2_function_type )(  ) const;
            point2_function_type point2_function_value( &::SireSystem::DihedralComponent::point2 );
            
            DihedralComponent_exposer.def( 
                "point2"
                , point2_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::point3
        
            typedef ::SireFF::Point const & ( ::SireSystem::DihedralComponent::*point3_function_type )(  ) const;
            point3_function_type point3_function_value( &::SireSystem::DihedralComponent::point3 );
            
            DihedralComponent_exposer.def( 
                "point3"
                , point3_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::r01
        
            typedef ::SireCAS::Symbol const & ( *r01_function_type )(  );
            r01_function_type r01_function_value( &::SireSystem::DihedralComponent::r01 );
            
            DihedralComponent_exposer.def( 
                "r01"
                , r01_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::r03
        
            typedef ::SireCAS::Symbol const & ( *r03_function_type )(  );
            r03_function_type r03_function_value( &::SireSystem::DihedralComponent::r03 );
            
            DihedralComponent_exposer.def( 
                "r03"
                , r03_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::r12
        
            typedef ::SireCAS::Symbol const & ( *r12_function_type )(  );
            r12_function_type r12_function_value( &::SireSystem::DihedralComponent::r12 );
            
            DihedralComponent_exposer.def( 
                "r12"
                , r12_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::r23
        
            typedef ::SireCAS::Symbol const & ( *r23_function_type )(  );
            r23_function_type r23_function_value( &::SireSystem::DihedralComponent::r23 );
            
            DihedralComponent_exposer.def( 
                "r23"
                , r23_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::theta012
        
            typedef ::SireCAS::Symbol const & ( *theta012_function_type )(  );
            theta012_function_type theta012_function_value( &::SireSystem::DihedralComponent::theta012 );
            
            DihedralComponent_exposer.def( 
                "theta012"
                , theta012_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::theta123
        
            typedef ::SireCAS::Symbol const & ( *theta123_function_type )(  );
            theta123_function_type theta123_function_value( &::SireSystem::DihedralComponent::theta123 );
            
            DihedralComponent_exposer.def( 
                "theta123"
                , theta123_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::DihedralComponent::toString
        
            typedef ::QString ( ::SireSystem::DihedralComponent::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireSystem::DihedralComponent::toString );
            
            DihedralComponent_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireSystem::DihedralComponent::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::DihedralComponent::typeName );
            
            DihedralComponent_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        DihedralComponent_exposer.staticmethod( "phi" );
        DihedralComponent_exposer.staticmethod( "r01" );
        DihedralComponent_exposer.staticmethod( "r03" );
        DihedralComponent_exposer.staticmethod( "r12" );
        DihedralComponent_exposer.staticmethod( "r23" );
        DihedralComponent_exposer.staticmethod( "theta012" );
        DihedralComponent_exposer.staticmethod( "theta123" );
        DihedralComponent_exposer.staticmethod( "typeName" );
        DihedralComponent_exposer.def( "__copy__", &__copy__);
        DihedralComponent_exposer.def( "__deepcopy__", &__copy__);
        DihedralComponent_exposer.def( "clone", &__copy__);
        DihedralComponent_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::DihedralComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DihedralComponent_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::DihedralComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DihedralComponent_exposer.def( "__str__", &__str__< ::SireSystem::DihedralComponent > );
        DihedralComponent_exposer.def( "__repr__", &__str__< ::SireSystem::DihedralComponent > );
    }

}
