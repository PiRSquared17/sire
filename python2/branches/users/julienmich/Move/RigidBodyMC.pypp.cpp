// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "RigidBodyMC.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/quaternion.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "SireVol/space.h"

#include "ensemble.h"

#include "rigidbodymc.h"

#include "uniformsampler.h"

#include <QDebug>

#include <QTime>

#include "rigidbodymc.h"

SireMove::RigidBodyMC __copy__(const SireMove::RigidBodyMC &other){ return SireMove::RigidBodyMC(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_RigidBodyMC_class(){

    { //::SireMove::RigidBodyMC
        typedef bp::class_< SireMove::RigidBodyMC, bp::bases< SireMove::MonteCarlo, SireMove::Move, SireBase::Property > > RigidBodyMC_exposer_t;
        RigidBodyMC_exposer_t RigidBodyMC_exposer = RigidBodyMC_exposer_t( "RigidBodyMC" );
        bp::scope RigidBodyMC_scope( RigidBodyMC_exposer );
        RigidBodyMC_exposer.def( bp::init< >() );
        RigidBodyMC_exposer.def( bp::init< SireMol::MoleculeGroup const & >(( bp::arg("molgroup") )) );
        RigidBodyMC_exposer.def( bp::init< SireMove::Sampler const & >(( bp::arg("sampler") )) );
        RigidBodyMC_exposer.def( bp::init< SireMove::RigidBodyMC const & >(( bp::arg("other") )) );
        { //::SireMove::RigidBodyMC::maximumRotation
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::RigidBodyMC::*maximumRotation_function_type )(  ) const;
            maximumRotation_function_type maximumRotation_function_value( &::SireMove::RigidBodyMC::maximumRotation );
            
            RigidBodyMC_exposer.def( 
                "maximumRotation"
                , maximumRotation_function_value );
        
        }
        { //::SireMove::RigidBodyMC::maximumTranslation
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::RigidBodyMC::*maximumTranslation_function_type )(  ) const;
            maximumTranslation_function_type maximumTranslation_function_value( &::SireMove::RigidBodyMC::maximumTranslation );
            
            RigidBodyMC_exposer.def( 
                "maximumTranslation"
                , maximumTranslation_function_value );
        
        }
        { //::SireMove::RigidBodyMC::moleculeGroup
        
            typedef ::SireMol::MoleculeGroup const & ( ::SireMove::RigidBodyMC::*moleculeGroup_function_type )(  ) const;
            moleculeGroup_function_type moleculeGroup_function_value( &::SireMove::RigidBodyMC::moleculeGroup );
            
            RigidBodyMC_exposer.def( 
                "moleculeGroup"
                , moleculeGroup_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMove::RigidBodyMC::move
        
            typedef void ( ::SireMove::RigidBodyMC::*move_function_type )( ::SireSystem::System &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::RigidBodyMC::move );
            
            RigidBodyMC_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) ) );
        
        }
        RigidBodyMC_exposer.def( bp::self != bp::self );
        { //::SireMove::RigidBodyMC::operator=
        
            typedef ::SireMove::RigidBodyMC & ( ::SireMove::RigidBodyMC::*assign_function_type )( ::SireMove::RigidBodyMC const & ) ;
            assign_function_type assign_function_value( &::SireMove::RigidBodyMC::operator= );
            
            RigidBodyMC_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        RigidBodyMC_exposer.def( bp::self == bp::self );
        { //::SireMove::RigidBodyMC::sampler
        
            typedef ::SireMove::Sampler const & ( ::SireMove::RigidBodyMC::*sampler_function_type )(  ) const;
            sampler_function_type sampler_function_value( &::SireMove::RigidBodyMC::sampler );
            
            RigidBodyMC_exposer.def( 
                "sampler"
                , sampler_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMove::RigidBodyMC::setGenerator
        
            typedef void ( ::SireMove::RigidBodyMC::*setGenerator_function_type )( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::RigidBodyMC::setGenerator );
            
            RigidBodyMC_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("rangenerator") ) );
        
        }
        { //::SireMove::RigidBodyMC::setMaximumRotation
        
            typedef void ( ::SireMove::RigidBodyMC::*setMaximumRotation_function_type )( ::SireUnits::Dimension::Angle ) ;
            setMaximumRotation_function_type setMaximumRotation_function_value( &::SireMove::RigidBodyMC::setMaximumRotation );
            
            RigidBodyMC_exposer.def( 
                "setMaximumRotation"
                , setMaximumRotation_function_value
                , ( bp::arg("max_rotation") ) );
        
        }
        { //::SireMove::RigidBodyMC::setMaximumTranslation
        
            typedef void ( ::SireMove::RigidBodyMC::*setMaximumTranslation_function_type )( ::SireUnits::Dimension::Length ) ;
            setMaximumTranslation_function_type setMaximumTranslation_function_value( &::SireMove::RigidBodyMC::setMaximumTranslation );
            
            RigidBodyMC_exposer.def( 
                "setMaximumTranslation"
                , setMaximumTranslation_function_value
                , ( bp::arg("max_translation") ) );
        
        }
        { //::SireMove::RigidBodyMC::setSampler
        
            typedef void ( ::SireMove::RigidBodyMC::*setSampler_function_type )( ::SireMove::Sampler const & ) ;
            setSampler_function_type setSampler_function_value( &::SireMove::RigidBodyMC::setSampler );
            
            RigidBodyMC_exposer.def( 
                "setSampler"
                , setSampler_function_value
                , ( bp::arg("sampler") ) );
        
        }
        { //::SireMove::RigidBodyMC::setSampler
        
            typedef void ( ::SireMove::RigidBodyMC::*setSampler_function_type )( ::SireMol::MoleculeGroup const & ) ;
            setSampler_function_type setSampler_function_value( &::SireMove::RigidBodyMC::setSampler );
            
            RigidBodyMC_exposer.def( 
                "setSampler"
                , setSampler_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireMove::RigidBodyMC::setSharedRotationCenter
        
            typedef void ( ::SireMove::RigidBodyMC::*setSharedRotationCenter_function_type )( bool ) ;
            setSharedRotationCenter_function_type setSharedRotationCenter_function_value( &::SireMove::RigidBodyMC::setSharedRotationCenter );
            
            RigidBodyMC_exposer.def( 
                "setSharedRotationCenter"
                , setSharedRotationCenter_function_value
                , ( bp::arg("on") ) );
        
        }
        { //::SireMove::RigidBodyMC::setSynchronisedRotation
        
            typedef void ( ::SireMove::RigidBodyMC::*setSynchronisedRotation_function_type )( bool ) ;
            setSynchronisedRotation_function_type setSynchronisedRotation_function_value( &::SireMove::RigidBodyMC::setSynchronisedRotation );
            
            RigidBodyMC_exposer.def( 
                "setSynchronisedRotation"
                , setSynchronisedRotation_function_value
                , ( bp::arg("on") ) );
        
        }
        { //::SireMove::RigidBodyMC::setSynchronisedTranslation
        
            typedef void ( ::SireMove::RigidBodyMC::*setSynchronisedTranslation_function_type )( bool ) ;
            setSynchronisedTranslation_function_type setSynchronisedTranslation_function_value( &::SireMove::RigidBodyMC::setSynchronisedTranslation );
            
            RigidBodyMC_exposer.def( 
                "setSynchronisedTranslation"
                , setSynchronisedTranslation_function_value
                , ( bp::arg("on") ) );
        
        }
        { //::SireMove::RigidBodyMC::sharedRotationCenter
        
            typedef bool ( ::SireMove::RigidBodyMC::*sharedRotationCenter_function_type )(  ) const;
            sharedRotationCenter_function_type sharedRotationCenter_function_value( &::SireMove::RigidBodyMC::sharedRotationCenter );
            
            RigidBodyMC_exposer.def( 
                "sharedRotationCenter"
                , sharedRotationCenter_function_value );
        
        }
        { //::SireMove::RigidBodyMC::synchronisedRotation
        
            typedef bool ( ::SireMove::RigidBodyMC::*synchronisedRotation_function_type )(  ) const;
            synchronisedRotation_function_type synchronisedRotation_function_value( &::SireMove::RigidBodyMC::synchronisedRotation );
            
            RigidBodyMC_exposer.def( 
                "synchronisedRotation"
                , synchronisedRotation_function_value );
        
        }
        { //::SireMove::RigidBodyMC::synchronisedTranslation
        
            typedef bool ( ::SireMove::RigidBodyMC::*synchronisedTranslation_function_type )(  ) const;
            synchronisedTranslation_function_type synchronisedTranslation_function_value( &::SireMove::RigidBodyMC::synchronisedTranslation );
            
            RigidBodyMC_exposer.def( 
                "synchronisedTranslation"
                , synchronisedTranslation_function_value );
        
        }
        { //::SireMove::RigidBodyMC::toString
        
            typedef ::QString ( ::SireMove::RigidBodyMC::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMove::RigidBodyMC::toString );
            
            RigidBodyMC_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMove::RigidBodyMC::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::RigidBodyMC::typeName );
            
            RigidBodyMC_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        RigidBodyMC_exposer.staticmethod( "typeName" );
        RigidBodyMC_exposer.def( "__copy__", &__copy__);
        RigidBodyMC_exposer.def( "__deepcopy__", &__copy__);
        RigidBodyMC_exposer.def( "clone", &__copy__);
        RigidBodyMC_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::RigidBodyMC >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RigidBodyMC_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::RigidBodyMC >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RigidBodyMC_exposer.def( "__str__", &__str__< ::SireMove::RigidBodyMC > );
        RigidBodyMC_exposer.def( "__repr__", &__str__< ::SireMove::RigidBodyMC > );
    }

}
