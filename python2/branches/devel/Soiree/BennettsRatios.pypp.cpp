// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "BennettsRatios.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMaths/maths.h"

#include "SireStream/shareddatastream.h"

#include "bennetts.h"

#include "tostring.h"

#include "bennetts.h"

Soiree::BennettsRatios __copy__(const Soiree::BennettsRatios &other){ return Soiree::BennettsRatios(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_BennettsRatios_class(){

    { //::Soiree::BennettsRatios
        typedef bp::class_< Soiree::BennettsRatios, bp::bases< SireBase::Property > > BennettsRatios_exposer_t;
        BennettsRatios_exposer_t BennettsRatios_exposer = BennettsRatios_exposer_t( "BennettsRatios", bp::init< >() );
        bp::scope BennettsRatios_scope( BennettsRatios_exposer );
        BennettsRatios_exposer.def( bp::init< QList< double > const &, QMap< double, SireMaths::BennettsFreeEnergyAverage > const &, QMap< double, SireMaths::BennettsFreeEnergyAverage > const & >(( bp::arg("windows"), bp::arg("forwards_ratios"), bp::arg("backwards_ratios") )) );
        BennettsRatios_exposer.def( bp::init< Soiree::BennettsRatios const & >(( bp::arg("other") )) );
        { //::Soiree::BennettsRatios::backwardsData
        
            typedef ::QMap< double, SireMaths::BennettsFreeEnergyAverage > ( ::Soiree::BennettsRatios::*backwardsData_function_type )(  ) const;
            backwardsData_function_type backwardsData_function_value( &::Soiree::BennettsRatios::backwardsData );
            
            BennettsRatios_exposer.def( 
                "backwardsData"
                , backwardsData_function_value );
        
        }
        { //::Soiree::BennettsRatios::backwardsRatios
        
            typedef ::QMap< double, SireMaths::BennettsFreeEnergyAverage > ( ::Soiree::BennettsRatios::*backwardsRatios_function_type )(  ) const;
            backwardsRatios_function_type backwardsRatios_function_value( &::Soiree::BennettsRatios::backwardsRatios );
            
            BennettsRatios_exposer.def( 
                "backwardsRatios"
                , backwardsRatios_function_value );
        
        }
        { //::Soiree::BennettsRatios::denominators
        
            typedef ::QVector< Soiree::DataPoint > ( ::Soiree::BennettsRatios::*denominators_function_type )(  ) const;
            denominators_function_type denominators_function_value( &::Soiree::BennettsRatios::denominators );
            
            BennettsRatios_exposer.def( 
                "denominators"
                , denominators_function_value );
        
        }
        { //::Soiree::BennettsRatios::forwardsData
        
            typedef ::QMap< double, SireMaths::BennettsFreeEnergyAverage > ( ::Soiree::BennettsRatios::*forwardsData_function_type )(  ) const;
            forwardsData_function_type forwardsData_function_value( &::Soiree::BennettsRatios::forwardsData );
            
            BennettsRatios_exposer.def( 
                "forwardsData"
                , forwardsData_function_value );
        
        }
        { //::Soiree::BennettsRatios::forwardsRatios
        
            typedef ::QMap< double, SireMaths::BennettsFreeEnergyAverage > ( ::Soiree::BennettsRatios::*forwardsRatios_function_type )(  ) const;
            forwardsRatios_function_type forwardsRatios_function_value( &::Soiree::BennettsRatios::forwardsRatios );
            
            BennettsRatios_exposer.def( 
                "forwardsRatios"
                , forwardsRatios_function_value );
        
        }
        { //::Soiree::BennettsRatios::integrate
        
            typedef ::Soiree::PMF ( ::Soiree::BennettsRatios::*integrate_function_type )(  ) const;
            integrate_function_type integrate_function_value( &::Soiree::BennettsRatios::integrate );
            
            BennettsRatios_exposer.def( 
                "integrate"
                , integrate_function_value );
        
        }
        { //::Soiree::BennettsRatios::isEmpty
        
            typedef bool ( ::Soiree::BennettsRatios::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::Soiree::BennettsRatios::isEmpty );
            
            BennettsRatios_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::Soiree::BennettsRatios::lambdaValues
        
            typedef ::QList< double > ( ::Soiree::BennettsRatios::*lambdaValues_function_type )(  ) const;
            lambdaValues_function_type lambdaValues_function_value( &::Soiree::BennettsRatios::lambdaValues );
            
            BennettsRatios_exposer.def( 
                "lambdaValues"
                , lambdaValues_function_value );
        
        }
        { //::Soiree::BennettsRatios::merge
        
            typedef ::Soiree::BennettsRatios ( *merge_function_type )( ::QList< Soiree::BennettsRatios > const & );
            merge_function_type merge_function_value( &::Soiree::BennettsRatios::merge );
            
            BennettsRatios_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("deltas") ) );
        
        }
        { //::Soiree::BennettsRatios::nLambdaValues
        
            typedef int ( ::Soiree::BennettsRatios::*nLambdaValues_function_type )(  ) const;
            nLambdaValues_function_type nLambdaValues_function_value( &::Soiree::BennettsRatios::nLambdaValues );
            
            BennettsRatios_exposer.def( 
                "nLambdaValues"
                , nLambdaValues_function_value );
        
        }
        { //::Soiree::BennettsRatios::nSamples
        
            typedef ::qint64 ( ::Soiree::BennettsRatios::*nSamples_function_type )(  ) const;
            nSamples_function_type nSamples_function_value( &::Soiree::BennettsRatios::nSamples );
            
            BennettsRatios_exposer.def( 
                "nSamples"
                , nSamples_function_value );
        
        }
        { //::Soiree::BennettsRatios::nWindows
        
            typedef int ( ::Soiree::BennettsRatios::*nWindows_function_type )(  ) const;
            nWindows_function_type nWindows_function_value( &::Soiree::BennettsRatios::nWindows );
            
            BennettsRatios_exposer.def( 
                "nWindows"
                , nWindows_function_value );
        
        }
        { //::Soiree::BennettsRatios::numerators
        
            typedef ::QVector< Soiree::DataPoint > ( ::Soiree::BennettsRatios::*numerators_function_type )(  ) const;
            numerators_function_type numerators_function_value( &::Soiree::BennettsRatios::numerators );
            
            BennettsRatios_exposer.def( 
                "numerators"
                , numerators_function_value );
        
        }
        BennettsRatios_exposer.def( bp::self != bp::self );
        BennettsRatios_exposer.def( bp::self + bp::self );
        { //::Soiree::BennettsRatios::operator=
        
            typedef ::Soiree::BennettsRatios & ( ::Soiree::BennettsRatios::*assign_function_type )( ::Soiree::BennettsRatios const & ) ;
            assign_function_type assign_function_value( &::Soiree::BennettsRatios::operator= );
            
            BennettsRatios_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        BennettsRatios_exposer.def( bp::self == bp::self );
        { //::Soiree::BennettsRatios::sum
        
            typedef ::Soiree::PMF ( ::Soiree::BennettsRatios::*sum_function_type )(  ) const;
            sum_function_type sum_function_value( &::Soiree::BennettsRatios::sum );
            
            BennettsRatios_exposer.def( 
                "sum"
                , sum_function_value );
        
        }
        { //::Soiree::BennettsRatios::temperature
        
            typedef ::SireUnits::Dimension::Temperature ( ::Soiree::BennettsRatios::*temperature_function_type )(  ) const;
            temperature_function_type temperature_function_value( &::Soiree::BennettsRatios::temperature );
            
            BennettsRatios_exposer.def( 
                "temperature"
                , temperature_function_value );
        
        }
        { //::Soiree::BennettsRatios::toString
        
            typedef ::QString ( ::Soiree::BennettsRatios::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::Soiree::BennettsRatios::toString );
            
            BennettsRatios_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::Soiree::BennettsRatios::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Soiree::BennettsRatios::typeName );
            
            BennettsRatios_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::Soiree::BennettsRatios::values
        
            typedef ::QVector< Soiree::DataPoint > ( ::Soiree::BennettsRatios::*values_function_type )(  ) const;
            values_function_type values_function_value( &::Soiree::BennettsRatios::values );
            
            BennettsRatios_exposer.def( 
                "values"
                , values_function_value );
        
        }
        { //::Soiree::BennettsRatios::what
        
            typedef char const * ( ::Soiree::BennettsRatios::*what_function_type )(  ) const;
            what_function_type what_function_value( &::Soiree::BennettsRatios::what );
            
            BennettsRatios_exposer.def( 
                "what"
                , what_function_value );
        
        }
        { //::Soiree::BennettsRatios::windows
        
            typedef ::QList< double > ( ::Soiree::BennettsRatios::*windows_function_type )(  ) const;
            windows_function_type windows_function_value( &::Soiree::BennettsRatios::windows );
            
            BennettsRatios_exposer.def( 
                "windows"
                , windows_function_value );
        
        }
        BennettsRatios_exposer.staticmethod( "merge" );
        BennettsRatios_exposer.staticmethod( "typeName" );
        BennettsRatios_exposer.def( "__copy__", &__copy__);
        BennettsRatios_exposer.def( "__deepcopy__", &__copy__);
        BennettsRatios_exposer.def( "clone", &__copy__);
        BennettsRatios_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Soiree::BennettsRatios >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BennettsRatios_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Soiree::BennettsRatios >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BennettsRatios_exposer.def( "__str__", &__str__< ::Soiree::BennettsRatios > );
        BennettsRatios_exposer.def( "__repr__", &__str__< ::Soiree::BennettsRatios > );
    }

}
