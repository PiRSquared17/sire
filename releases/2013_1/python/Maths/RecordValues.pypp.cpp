// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "RecordValues.pypp.hpp"

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

SireMaths::RecordValues __copy__(const SireMaths::RecordValues &other){ return SireMaths::RecordValues(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_RecordValues_class(){

    { //::SireMaths::RecordValues
        typedef bp::class_< SireMaths::RecordValues, bp::bases< SireMaths::Accumulator, SireBase::Property > > RecordValues_exposer_t;
        RecordValues_exposer_t RecordValues_exposer = RecordValues_exposer_t( "RecordValues", bp::init< >() );
        bp::scope RecordValues_scope( RecordValues_exposer );
        RecordValues_exposer.def( bp::init< SireMaths::RecordValues const & >(( bp::arg("other") )) );
        { //::SireMaths::RecordValues::accumulate
        
            typedef void ( ::SireMaths::RecordValues::*accumulate_function_type )( double ) ;
            accumulate_function_type accumulate_function_value( &::SireMaths::RecordValues::accumulate );
            
            RecordValues_exposer.def( 
                "accumulate"
                , accumulate_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMaths::RecordValues::clear
        
            typedef void ( ::SireMaths::RecordValues::*clear_function_type )(  ) ;
            clear_function_type clear_function_value( &::SireMaths::RecordValues::clear );
            
            RecordValues_exposer.def( 
                "clear"
                , clear_function_value );
        
        }
        { //::SireMaths::RecordValues::count
        
            typedef int ( ::SireMaths::RecordValues::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireMaths::RecordValues::count );
            
            RecordValues_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMaths::RecordValues::max
        
            typedef double ( ::SireMaths::RecordValues::*max_function_type )(  ) const;
            max_function_type max_function_value( &::SireMaths::RecordValues::max );
            
            RecordValues_exposer.def( 
                "max"
                , max_function_value );
        
        }
        { //::SireMaths::RecordValues::maximum
        
            typedef double ( ::SireMaths::RecordValues::*maximum_function_type )(  ) const;
            maximum_function_type maximum_function_value( &::SireMaths::RecordValues::maximum );
            
            RecordValues_exposer.def( 
                "maximum"
                , maximum_function_value );
        
        }
        { //::SireMaths::RecordValues::mean
        
            typedef double ( ::SireMaths::RecordValues::*mean_function_type )(  ) const;
            mean_function_type mean_function_value( &::SireMaths::RecordValues::mean );
            
            RecordValues_exposer.def( 
                "mean"
                , mean_function_value );
        
        }
        { //::SireMaths::RecordValues::meanOfSquares
        
            typedef double ( ::SireMaths::RecordValues::*meanOfSquares_function_type )(  ) const;
            meanOfSquares_function_type meanOfSquares_function_value( &::SireMaths::RecordValues::meanOfSquares );
            
            RecordValues_exposer.def( 
                "meanOfSquares"
                , meanOfSquares_function_value );
        
        }
        { //::SireMaths::RecordValues::median
        
            typedef double ( ::SireMaths::RecordValues::*median_function_type )(  ) const;
            median_function_type median_function_value( &::SireMaths::RecordValues::median );
            
            RecordValues_exposer.def( 
                "median"
                , median_function_value );
        
        }
        { //::SireMaths::RecordValues::min
        
            typedef double ( ::SireMaths::RecordValues::*min_function_type )(  ) const;
            min_function_type min_function_value( &::SireMaths::RecordValues::min );
            
            RecordValues_exposer.def( 
                "min"
                , min_function_value );
        
        }
        { //::SireMaths::RecordValues::minimum
        
            typedef double ( ::SireMaths::RecordValues::*minimum_function_type )(  ) const;
            minimum_function_type minimum_function_value( &::SireMaths::RecordValues::minimum );
            
            RecordValues_exposer.def( 
                "minimum"
                , minimum_function_value );
        
        }
        { //::SireMaths::RecordValues::nValues
        
            typedef int ( ::SireMaths::RecordValues::*nValues_function_type )(  ) const;
            nValues_function_type nValues_function_value( &::SireMaths::RecordValues::nValues );
            
            RecordValues_exposer.def( 
                "nValues"
                , nValues_function_value );
        
        }
        RecordValues_exposer.def( bp::self != bp::self );
        { //::SireMaths::RecordValues::operator=
        
            typedef ::SireMaths::RecordValues & ( ::SireMaths::RecordValues::*assign_function_type )( ::SireMaths::RecordValues const & ) ;
            assign_function_type assign_function_value( &::SireMaths::RecordValues::operator= );
            
            RecordValues_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        RecordValues_exposer.def( bp::self == bp::self );
        { //::SireMaths::RecordValues::size
        
            typedef int ( ::SireMaths::RecordValues::*size_function_type )(  ) const;
            size_function_type size_function_value( &::SireMaths::RecordValues::size );
            
            RecordValues_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMaths::RecordValues::standardDeviation
        
            typedef double ( ::SireMaths::RecordValues::*standardDeviation_function_type )(  ) const;
            standardDeviation_function_type standardDeviation_function_value( &::SireMaths::RecordValues::standardDeviation );
            
            RecordValues_exposer.def( 
                "standardDeviation"
                , standardDeviation_function_value );
        
        }
        { //::SireMaths::RecordValues::stddev
        
            typedef double ( ::SireMaths::RecordValues::*stddev_function_type )(  ) const;
            stddev_function_type stddev_function_value( &::SireMaths::RecordValues::stddev );
            
            RecordValues_exposer.def( 
                "stddev"
                , stddev_function_value );
        
        }
        { //::SireMaths::RecordValues::sum
        
            typedef double ( ::SireMaths::RecordValues::*sum_function_type )(  ) const;
            sum_function_type sum_function_value( &::SireMaths::RecordValues::sum );
            
            RecordValues_exposer.def( 
                "sum"
                , sum_function_value );
        
        }
        { //::SireMaths::RecordValues::sum2
        
            typedef double ( ::SireMaths::RecordValues::*sum2_function_type )(  ) const;
            sum2_function_type sum2_function_value( &::SireMaths::RecordValues::sum2 );
            
            RecordValues_exposer.def( 
                "sum2"
                , sum2_function_value );
        
        }
        { //::SireMaths::RecordValues::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::RecordValues::typeName );
            
            RecordValues_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMaths::RecordValues::values
        
            typedef ::QVector< double > ( ::SireMaths::RecordValues::*values_function_type )(  ) const;
            values_function_type values_function_value( &::SireMaths::RecordValues::values );
            
            RecordValues_exposer.def( 
                "values"
                , values_function_value );
        
        }
        RecordValues_exposer.staticmethod( "typeName" );
        RecordValues_exposer.def( "__copy__", &__copy__);
        RecordValues_exposer.def( "__deepcopy__", &__copy__);
        RecordValues_exposer.def( "clone", &__copy__);
        RecordValues_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::RecordValues >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RecordValues_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::RecordValues >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RecordValues_exposer.def( "__str__", &__str__< ::SireMaths::RecordValues > );
        RecordValues_exposer.def( "__repr__", &__str__< ::SireMaths::RecordValues > );
    }

}
