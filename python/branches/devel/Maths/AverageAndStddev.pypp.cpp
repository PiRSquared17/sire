// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AverageAndStddev.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/maths.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "accumulator.h"

#include <QDebug>

#include <QMutex>

#include <cmath>

#include "accumulator.h"

SireMaths::AverageAndStddev __copy__(const SireMaths::AverageAndStddev &other){ return SireMaths::AverageAndStddev(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AverageAndStddev_class(){

    { //::SireMaths::AverageAndStddev
        typedef bp::class_< SireMaths::AverageAndStddev, bp::bases< SireMaths::Average, SireMaths::Accumulator, SireBase::Property > > AverageAndStddev_exposer_t;
        AverageAndStddev_exposer_t AverageAndStddev_exposer = AverageAndStddev_exposer_t( "AverageAndStddev", bp::init< >() );
        bp::scope AverageAndStddev_scope( AverageAndStddev_exposer );
        AverageAndStddev_exposer.def( bp::init< SireMaths::AverageAndStddev const & >(( bp::arg("other") )) );
        { //::SireMaths::AverageAndStddev::accumulate
        
            typedef void ( ::SireMaths::AverageAndStddev::*accumulate_function_type )( double ) ;
            accumulate_function_type accumulate_function_value( &::SireMaths::AverageAndStddev::accumulate );
            
            AverageAndStddev_exposer.def( 
                "accumulate"
                , accumulate_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMaths::AverageAndStddev::clear
        
            typedef void ( ::SireMaths::AverageAndStddev::*clear_function_type )(  ) ;
            clear_function_type clear_function_value( &::SireMaths::AverageAndStddev::clear );
            
            AverageAndStddev_exposer.def( 
                "clear"
                , clear_function_value );
        
        }
        { //::SireMaths::AverageAndStddev::meanOfSquares
        
            typedef double ( ::SireMaths::AverageAndStddev::*meanOfSquares_function_type )(  ) const;
            meanOfSquares_function_type meanOfSquares_function_value( &::SireMaths::AverageAndStddev::meanOfSquares );
            
            AverageAndStddev_exposer.def( 
                "meanOfSquares"
                , meanOfSquares_function_value );
        
        }
        AverageAndStddev_exposer.def( bp::self != bp::self );
        AverageAndStddev_exposer.def( bp::self + bp::self );
        { //::SireMaths::AverageAndStddev::operator=
        
            typedef ::SireMaths::AverageAndStddev & ( ::SireMaths::AverageAndStddev::*assign_function_type )( ::SireMaths::AverageAndStddev const & ) ;
            assign_function_type assign_function_value( &::SireMaths::AverageAndStddev::operator= );
            
            AverageAndStddev_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        AverageAndStddev_exposer.def( bp::self == bp::self );
        { //::SireMaths::AverageAndStddev::standardDeviation
        
            typedef double ( ::SireMaths::AverageAndStddev::*standardDeviation_function_type )(  ) const;
            standardDeviation_function_type standardDeviation_function_value( &::SireMaths::AverageAndStddev::standardDeviation );
            
            AverageAndStddev_exposer.def( 
                "standardDeviation"
                , standardDeviation_function_value );
        
        }
        { //::SireMaths::AverageAndStddev::stddev
        
            typedef double ( ::SireMaths::AverageAndStddev::*stddev_function_type )(  ) const;
            stddev_function_type stddev_function_value( &::SireMaths::AverageAndStddev::stddev );
            
            AverageAndStddev_exposer.def( 
                "stddev"
                , stddev_function_value );
        
        }
        { //::SireMaths::AverageAndStddev::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::AverageAndStddev::typeName );
            
            AverageAndStddev_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        AverageAndStddev_exposer.staticmethod( "typeName" );
        AverageAndStddev_exposer.def( "__copy__", &__copy__);
        AverageAndStddev_exposer.def( "__deepcopy__", &__copy__);
        AverageAndStddev_exposer.def( "clone", &__copy__);
        AverageAndStddev_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::AverageAndStddev >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AverageAndStddev_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::AverageAndStddev >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AverageAndStddev_exposer.def( "__str__", &__str__< ::SireMaths::AverageAndStddev > );
        AverageAndStddev_exposer.def( "__repr__", &__str__< ::SireMaths::AverageAndStddev > );
    }

}
