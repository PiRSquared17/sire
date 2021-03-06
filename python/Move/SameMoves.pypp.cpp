// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "SameMoves.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/rangenerator.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/constraints.h"

#include "SireSystem/errors.h"

#include "SireSystem/system.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "ensemble.h"

#include "moves.h"

#include "simulation.h"

#include <QDebug>

#include <QMutex>

#include "moves.h"

SireMove::SameMoves __copy__(const SireMove::SameMoves &other){ return SireMove::SameMoves(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_SameMoves_class(){

    { //::SireMove::SameMoves
        typedef bp::class_< SireMove::SameMoves, bp::bases< SireMove::Moves, SireBase::Property > > SameMoves_exposer_t;
        SameMoves_exposer_t SameMoves_exposer = SameMoves_exposer_t( "SameMoves", bp::init< >() );
        bp::scope SameMoves_scope( SameMoves_exposer );
        SameMoves_exposer.def( bp::init< SireMove::Move const & >(( bp::arg("move") )) );
        SameMoves_exposer.def( bp::init< SireMove::SameMoves const & >(( bp::arg("other") )) );
        { //::SireMove::SameMoves::clearStatistics
        
            typedef void ( ::SireMove::SameMoves::*clearStatistics_function_type )(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireMove::SameMoves::clearStatistics );
            
            SameMoves_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value );
        
        }
        { //::SireMove::SameMoves::energyComponent
        
            typedef ::SireCAS::Symbol const & ( ::SireMove::SameMoves::*energyComponent_function_type )(  ) const;
            energyComponent_function_type energyComponent_function_value( &::SireMove::SameMoves::energyComponent );
            
            SameMoves_exposer.def( 
                "energyComponent"
                , energyComponent_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMove::SameMoves::move
        
            typedef ::SireSystem::System ( ::SireMove::SameMoves::*move_function_type )( ::SireSystem::System const &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::SameMoves::move );
            
            SameMoves_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves")=(int)(1), bp::arg("record_stats")=(bool)(false) ) );
        
        }
        { //::SireMove::SameMoves::moves
        
            typedef ::QList< SireBase::PropPtr< SireMove::Move > > ( ::SireMove::SameMoves::*moves_function_type )(  ) const;
            moves_function_type moves_function_value( &::SireMove::SameMoves::moves );
            
            SameMoves_exposer.def( 
                "moves"
                , moves_function_value );
        
        }
        SameMoves_exposer.def( bp::self != bp::self );
        { //::SireMove::SameMoves::operator=
        
            typedef ::SireMove::SameMoves & ( ::SireMove::SameMoves::*assign_function_type )( ::SireMove::SameMoves const & ) ;
            assign_function_type assign_function_value( &::SireMove::SameMoves::operator= );
            
            SameMoves_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        SameMoves_exposer.def( bp::self == bp::self );
        { //::SireMove::SameMoves::setEnergyComponent
        
            typedef void ( ::SireMove::SameMoves::*setEnergyComponent_function_type )( ::SireCAS::Symbol const & ) ;
            setEnergyComponent_function_type setEnergyComponent_function_value( &::SireMove::SameMoves::setEnergyComponent );
            
            SameMoves_exposer.def( 
                "setEnergyComponent"
                , setEnergyComponent_function_value
                , ( bp::arg("component") ) );
        
        }
        { //::SireMove::SameMoves::setGenerator
        
            typedef void ( ::SireMove::SameMoves::*setGenerator_function_type )( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::SameMoves::setGenerator );
            
            SameMoves_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("rangenerator") ) );
        
        }
        { //::SireMove::SameMoves::setSpaceProperty
        
            typedef void ( ::SireMove::SameMoves::*setSpaceProperty_function_type )( ::SireBase::PropertyName const & ) ;
            setSpaceProperty_function_type setSpaceProperty_function_value( &::SireMove::SameMoves::setSpaceProperty );
            
            SameMoves_exposer.def( 
                "setSpaceProperty"
                , setSpaceProperty_function_value
                , ( bp::arg("spaceproperty") ) );
        
        }
        { //::SireMove::SameMoves::spaceProperty
        
            typedef ::SireBase::PropertyName const & ( ::SireMove::SameMoves::*spaceProperty_function_type )(  ) const;
            spaceProperty_function_type spaceProperty_function_value( &::SireMove::SameMoves::spaceProperty );
            
            SameMoves_exposer.def( 
                "spaceProperty"
                , spaceProperty_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::SameMoves::toString
        
            typedef ::QString ( ::SireMove::SameMoves::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMove::SameMoves::toString );
            
            SameMoves_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMove::SameMoves::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::SameMoves::typeName );
            
            SameMoves_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        SameMoves_exposer.staticmethod( "typeName" );
        SameMoves_exposer.def( "__copy__", &__copy__);
        SameMoves_exposer.def( "__deepcopy__", &__copy__);
        SameMoves_exposer.def( "clone", &__copy__);
        SameMoves_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::SameMoves >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SameMoves_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::SameMoves >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SameMoves_exposer.def( "__str__", &__str__< ::SireMove::SameMoves > );
        SameMoves_exposer.def( "__repr__", &__str__< ::SireMove::SameMoves > );
        SameMoves_exposer.def( "__len__", &__len_size< ::SireMove::SameMoves > );
    }

}
