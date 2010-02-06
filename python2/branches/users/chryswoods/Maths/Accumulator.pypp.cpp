// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Accumulator.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/maths.h"

#include "Siren/errors.h"

#include "Siren/objref.h"

#include "Siren/stream.h"

#include "Siren/streamqt.h"

#include "accumulator.h"

#include <QDebug>

#include <cmath>

#include "accumulator.h"

#include "Siren/logger.h"

#include "Siren/str.hpp"

void register_Accumulator_class(){

    { //::SireMaths::Accumulator
        typedef bp::class_< SireMaths::Accumulator, bp::bases< Siren::Mutable, Siren::Object >, boost::noncopyable > Accumulator_exposer_t;
        Accumulator_exposer_t Accumulator_exposer = Accumulator_exposer_t( "Accumulator", bp::no_init );
        bp::scope Accumulator_scope( Accumulator_exposer );
        { //::SireMaths::Accumulator::accumulate
        
            typedef void ( ::SireMaths::Accumulator::*accumulate_function_type )( double ) ;
            accumulate_function_type accumulate_function_value( &::SireMaths::Accumulator::accumulate );
            
            Accumulator_exposer.def( 
                "accumulate"
                , accumulate_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMaths::Accumulator::accumulate
        
            typedef void ( ::SireMaths::Accumulator::*accumulate_function_type )( ::QVector< double > const & ) ;
            accumulate_function_type accumulate_function_value( &::SireMaths::Accumulator::accumulate );
            
            Accumulator_exposer.def( 
                "accumulate"
                , accumulate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMaths::Accumulator::accumulate
        
            typedef void ( ::SireMaths::Accumulator::*accumulate_function_type )( ::QList< double > const & ) ;
            accumulate_function_type accumulate_function_value( &::SireMaths::Accumulator::accumulate );
            
            Accumulator_exposer.def( 
                "accumulate"
                , accumulate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMaths::Accumulator::clear
        
            typedef void ( ::SireMaths::Accumulator::*clear_function_type )(  ) ;
            clear_function_type clear_function_value( &::SireMaths::Accumulator::clear );
            
            Accumulator_exposer.def( 
                "clear"
                , clear_function_value );
        
        }
        { //::SireMaths::Accumulator::hashCode
        
            typedef ::uint ( ::SireMaths::Accumulator::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::SireMaths::Accumulator::hashCode );
            
            Accumulator_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::SireMaths::Accumulator::nSamples
        
            typedef int ( ::SireMaths::Accumulator::*nSamples_function_type )(  ) const;
            nSamples_function_type nSamples_function_value( &::SireMaths::Accumulator::nSamples );
            
            Accumulator_exposer.def( 
                "nSamples"
                , nSamples_function_value );
        
        }
        { //::SireMaths::Accumulator::restoreState
        
            typedef void ( ::SireMaths::Accumulator::*restoreState_function_type )( ::Siren::Object const & ) ;
            restoreState_function_type restoreState_function_value( &::SireMaths::Accumulator::restoreState );
            
            Accumulator_exposer.def( 
                "restoreState"
                , restoreState_function_value
                , ( bp::arg("object") ) );
        
        }
        { //::SireMaths::Accumulator::saveState
        
            typedef ::Siren::ObjRef ( ::SireMaths::Accumulator::*saveState_function_type )(  ) const;
            saveState_function_type saveState_function_value( &::SireMaths::Accumulator::saveState );
            
            Accumulator_exposer.def( 
                "saveState"
                , saveState_function_value );
        
        }
        { //::SireMaths::Accumulator::stream
        
            typedef void ( ::SireMaths::Accumulator::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireMaths::Accumulator::stream );
            
            Accumulator_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::SireMaths::Accumulator::typeName
        
            typedef ::QString ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::Accumulator::typeName );
            
            Accumulator_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMaths::Accumulator::value
        
            typedef double ( ::SireMaths::Accumulator::*value_function_type )(  ) const;
            value_function_type value_function_value( &::SireMaths::Accumulator::value );
            
            Accumulator_exposer.def( 
                "value"
                , value_function_value );
        
        }
        Accumulator_exposer.staticmethod( "typeName" );
        Accumulator_exposer.def( "__str__", &__str__< ::SireMaths::Accumulator > );
        Accumulator_exposer.def( "__repr__", &__str__< ::SireMaths::Accumulator > );
    }

}