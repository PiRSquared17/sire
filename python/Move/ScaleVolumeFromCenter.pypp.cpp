// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "ScaleVolumeFromCenter.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/moleculegroup.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/space.h"

#include "math.h"

#include "sire_config.h"

#include "volumechanger.h"

#include <cmath>

#include "volumechanger.h"

SireMove::ScaleVolumeFromCenter __copy__(const SireMove::ScaleVolumeFromCenter &other){ return SireMove::ScaleVolumeFromCenter(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ScaleVolumeFromCenter_class(){

    { //::SireMove::ScaleVolumeFromCenter
        typedef bp::class_< SireMove::ScaleVolumeFromCenter, bp::bases< SireMove::VolumeChanger, SireBase::Property > > ScaleVolumeFromCenter_exposer_t;
        ScaleVolumeFromCenter_exposer_t ScaleVolumeFromCenter_exposer = ScaleVolumeFromCenter_exposer_t( "ScaleVolumeFromCenter", bp::init< >() );
        bp::scope ScaleVolumeFromCenter_scope( ScaleVolumeFromCenter_exposer );
        ScaleVolumeFromCenter_exposer.def( bp::init< SireMol::MGID const & >(( bp::arg("mgid") )) );
        ScaleVolumeFromCenter_exposer.def( bp::init< SireMol::MoleculeGroup const & >(( bp::arg("molgroup") )) );
        ScaleVolumeFromCenter_exposer.def( bp::init< SireMol::MGID const &, SireFF::PointRef const & >(( bp::arg("mgid"), bp::arg("point") )) );
        ScaleVolumeFromCenter_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireFF::PointRef const & >(( bp::arg("molgroup"), bp::arg("point") )) );
        ScaleVolumeFromCenter_exposer.def( bp::init< SireMove::ScaleVolumeFromCenter const & >(( bp::arg("other") )) );
        { //::SireMove::ScaleVolumeFromCenter::center
        
            typedef ::SireFF::Point const & ( ::SireMove::ScaleVolumeFromCenter::*center_function_type )(  ) const;
            center_function_type center_function_value( &::SireMove::ScaleVolumeFromCenter::center );
            
            ScaleVolumeFromCenter_exposer.def( 
                "center"
                , center_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        ScaleVolumeFromCenter_exposer.def( bp::self != bp::self );
        { //::SireMove::ScaleVolumeFromCenter::operator=
        
            typedef ::SireMove::ScaleVolumeFromCenter & ( ::SireMove::ScaleVolumeFromCenter::*assign_function_type )( ::SireMove::ScaleVolumeFromCenter const & ) ;
            assign_function_type assign_function_value( &::SireMove::ScaleVolumeFromCenter::operator= );
            
            ScaleVolumeFromCenter_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        ScaleVolumeFromCenter_exposer.def( bp::self == bp::self );
        { //::SireMove::ScaleVolumeFromCenter::setCenter
        
            typedef void ( ::SireMove::ScaleVolumeFromCenter::*setCenter_function_type )( ::SireFF::PointRef const & ) ;
            setCenter_function_type setCenter_function_value( &::SireMove::ScaleVolumeFromCenter::setCenter );
            
            ScaleVolumeFromCenter_exposer.def( 
                "setCenter"
                , setCenter_function_value
                , ( bp::arg("center") ) );
        
        }
        { //::SireMove::ScaleVolumeFromCenter::setVolume
        
            typedef int ( ::SireMove::ScaleVolumeFromCenter::*setVolume_function_type )( ::SireSystem::System &,::SireUnits::Dimension::Volume const &,::SireBase::PropertyMap const & ) const;
            setVolume_function_type setVolume_function_value( &::SireMove::ScaleVolumeFromCenter::setVolume );
            
            ScaleVolumeFromCenter_exposer.def( 
                "setVolume"
                , setVolume_function_value
                , ( bp::arg("system"), bp::arg("volume"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        { //::SireMove::ScaleVolumeFromCenter::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::ScaleVolumeFromCenter::typeName );
            
            ScaleVolumeFromCenter_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        ScaleVolumeFromCenter_exposer.staticmethod( "typeName" );
        ScaleVolumeFromCenter_exposer.def( "__copy__", &__copy__);
        ScaleVolumeFromCenter_exposer.def( "__deepcopy__", &__copy__);
        ScaleVolumeFromCenter_exposer.def( "clone", &__copy__);
        ScaleVolumeFromCenter_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::ScaleVolumeFromCenter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ScaleVolumeFromCenter_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::ScaleVolumeFromCenter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ScaleVolumeFromCenter_exposer.def( "__str__", &__str__< ::SireMove::ScaleVolumeFromCenter > );
        ScaleVolumeFromCenter_exposer.def( "__repr__", &__str__< ::SireMove::ScaleVolumeFromCenter > );
    }

}
