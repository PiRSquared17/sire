// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "FEPDeltas.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMaths/maths.h"

#include "SireStream/shareddatastream.h"

#include "fep.h"

#include "tostring.h"

#include "fep.h"

Soiree::FEPDeltas __copy__(const Soiree::FEPDeltas &other){ return Soiree::FEPDeltas(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_FEPDeltas_class(){

    { //::Soiree::FEPDeltas
        typedef bp::class_< Soiree::FEPDeltas, bp::bases< SireBase::Property > > FEPDeltas_exposer_t;
        FEPDeltas_exposer_t FEPDeltas_exposer = FEPDeltas_exposer_t( "FEPDeltas", bp::init< >() );
        bp::scope FEPDeltas_scope( FEPDeltas_exposer );
        FEPDeltas_exposer.def( bp::init< QList< double > const &, QMap< double, SireMaths::FreeEnergyAverage > const & >(( bp::arg("windows"), bp::arg("deltas") )) );
        FEPDeltas_exposer.def( bp::init< QList< double > const &, QMap< double, SireMaths::FreeEnergyAverage > const &, QMap< double, SireMaths::FreeEnergyAverage > const & >(( bp::arg("windows"), bp::arg("forwards_deltas"), bp::arg("backwards_deltas") )) );
        FEPDeltas_exposer.def( bp::init< Soiree::FEPDeltas const & >(( bp::arg("other") )) );
        { //::Soiree::FEPDeltas::backwardsData
        
            typedef ::QMap< double, SireMaths::FreeEnergyAverage > ( ::Soiree::FEPDeltas::*backwardsData_function_type )(  ) const;
            backwardsData_function_type backwardsData_function_value( &::Soiree::FEPDeltas::backwardsData );
            
            FEPDeltas_exposer.def( 
                "backwardsData"
                , backwardsData_function_value );
        
        }
        { //::Soiree::FEPDeltas::backwardsDeltas
        
            typedef ::QMap< double, SireMaths::FreeEnergyAverage > ( ::Soiree::FEPDeltas::*backwardsDeltas_function_type )(  ) const;
            backwardsDeltas_function_type backwardsDeltas_function_value( &::Soiree::FEPDeltas::backwardsDeltas );
            
            FEPDeltas_exposer.def( 
                "backwardsDeltas"
                , backwardsDeltas_function_value );
        
        }
        { //::Soiree::FEPDeltas::backwardsValues
        
            typedef ::QVector< Soiree::DataPoint > ( ::Soiree::FEPDeltas::*backwardsValues_function_type )(  ) const;
            backwardsValues_function_type backwardsValues_function_value( &::Soiree::FEPDeltas::backwardsValues );
            
            FEPDeltas_exposer.def( 
                "backwardsValues"
                , backwardsValues_function_value );
        
        }
        { //::Soiree::FEPDeltas::forwardsData
        
            typedef ::QMap< double, SireMaths::FreeEnergyAverage > ( ::Soiree::FEPDeltas::*forwardsData_function_type )(  ) const;
            forwardsData_function_type forwardsData_function_value( &::Soiree::FEPDeltas::forwardsData );
            
            FEPDeltas_exposer.def( 
                "forwardsData"
                , forwardsData_function_value );
        
        }
        { //::Soiree::FEPDeltas::forwardsDeltas
        
            typedef ::QMap< double, SireMaths::FreeEnergyAverage > ( ::Soiree::FEPDeltas::*forwardsDeltas_function_type )(  ) const;
            forwardsDeltas_function_type forwardsDeltas_function_value( &::Soiree::FEPDeltas::forwardsDeltas );
            
            FEPDeltas_exposer.def( 
                "forwardsDeltas"
                , forwardsDeltas_function_value );
        
        }
        { //::Soiree::FEPDeltas::forwardsValues
        
            typedef ::QVector< Soiree::DataPoint > ( ::Soiree::FEPDeltas::*forwardsValues_function_type )(  ) const;
            forwardsValues_function_type forwardsValues_function_value( &::Soiree::FEPDeltas::forwardsValues );
            
            FEPDeltas_exposer.def( 
                "forwardsValues"
                , forwardsValues_function_value );
        
        }
        { //::Soiree::FEPDeltas::integrate
        
            typedef ::Soiree::PMF ( ::Soiree::FEPDeltas::*integrate_function_type )(  ) const;
            integrate_function_type integrate_function_value( &::Soiree::FEPDeltas::integrate );
            
            FEPDeltas_exposer.def( 
                "integrate"
                , integrate_function_value );
        
        }
        { //::Soiree::FEPDeltas::isEmpty
        
            typedef bool ( ::Soiree::FEPDeltas::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::Soiree::FEPDeltas::isEmpty );
            
            FEPDeltas_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::Soiree::FEPDeltas::lambdaValues
        
            typedef ::QList< double > ( ::Soiree::FEPDeltas::*lambdaValues_function_type )(  ) const;
            lambdaValues_function_type lambdaValues_function_value( &::Soiree::FEPDeltas::lambdaValues );
            
            FEPDeltas_exposer.def( 
                "lambdaValues"
                , lambdaValues_function_value );
        
        }
        { //::Soiree::FEPDeltas::merge
        
            typedef ::Soiree::FEPDeltas ( *merge_function_type )( ::QList< Soiree::FEPDeltas > const & );
            merge_function_type merge_function_value( &::Soiree::FEPDeltas::merge );
            
            FEPDeltas_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("deltas") ) );
        
        }
        { //::Soiree::FEPDeltas::nLambdaValues
        
            typedef int ( ::Soiree::FEPDeltas::*nLambdaValues_function_type )(  ) const;
            nLambdaValues_function_type nLambdaValues_function_value( &::Soiree::FEPDeltas::nLambdaValues );
            
            FEPDeltas_exposer.def( 
                "nLambdaValues"
                , nLambdaValues_function_value );
        
        }
        { //::Soiree::FEPDeltas::nSamples
        
            typedef ::qint64 ( ::Soiree::FEPDeltas::*nSamples_function_type )(  ) const;
            nSamples_function_type nSamples_function_value( &::Soiree::FEPDeltas::nSamples );
            
            FEPDeltas_exposer.def( 
                "nSamples"
                , nSamples_function_value );
        
        }
        { //::Soiree::FEPDeltas::nWindows
        
            typedef int ( ::Soiree::FEPDeltas::*nWindows_function_type )(  ) const;
            nWindows_function_type nWindows_function_value( &::Soiree::FEPDeltas::nWindows );
            
            FEPDeltas_exposer.def( 
                "nWindows"
                , nWindows_function_value );
        
        }
        FEPDeltas_exposer.def( bp::self != bp::self );
        FEPDeltas_exposer.def( bp::self + bp::self );
        { //::Soiree::FEPDeltas::operator=
        
            typedef ::Soiree::FEPDeltas & ( ::Soiree::FEPDeltas::*assign_function_type )( ::Soiree::FEPDeltas const & ) ;
            assign_function_type assign_function_value( &::Soiree::FEPDeltas::operator= );
            
            FEPDeltas_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        FEPDeltas_exposer.def( bp::self == bp::self );
        { //::Soiree::FEPDeltas::sum
        
            typedef ::Soiree::PMF ( ::Soiree::FEPDeltas::*sum_function_type )(  ) const;
            sum_function_type sum_function_value( &::Soiree::FEPDeltas::sum );
            
            FEPDeltas_exposer.def( 
                "sum"
                , sum_function_value );
        
        }
        { //::Soiree::FEPDeltas::temperature
        
            typedef ::SireUnits::Dimension::Temperature ( ::Soiree::FEPDeltas::*temperature_function_type )(  ) const;
            temperature_function_type temperature_function_value( &::Soiree::FEPDeltas::temperature );
            
            FEPDeltas_exposer.def( 
                "temperature"
                , temperature_function_value );
        
        }
        { //::Soiree::FEPDeltas::toString
        
            typedef ::QString ( ::Soiree::FEPDeltas::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::Soiree::FEPDeltas::toString );
            
            FEPDeltas_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::Soiree::FEPDeltas::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Soiree::FEPDeltas::typeName );
            
            FEPDeltas_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::Soiree::FEPDeltas::values
        
            typedef ::QVector< Soiree::DataPoint > ( ::Soiree::FEPDeltas::*values_function_type )(  ) const;
            values_function_type values_function_value( &::Soiree::FEPDeltas::values );
            
            FEPDeltas_exposer.def( 
                "values"
                , values_function_value );
        
        }
        { //::Soiree::FEPDeltas::what
        
            typedef char const * ( ::Soiree::FEPDeltas::*what_function_type )(  ) const;
            what_function_type what_function_value( &::Soiree::FEPDeltas::what );
            
            FEPDeltas_exposer.def( 
                "what"
                , what_function_value );
        
        }
        { //::Soiree::FEPDeltas::windows
        
            typedef ::QList< double > ( ::Soiree::FEPDeltas::*windows_function_type )(  ) const;
            windows_function_type windows_function_value( &::Soiree::FEPDeltas::windows );
            
            FEPDeltas_exposer.def( 
                "windows"
                , windows_function_value );
        
        }
        FEPDeltas_exposer.staticmethod( "merge" );
        FEPDeltas_exposer.staticmethod( "typeName" );
        FEPDeltas_exposer.def( "__copy__", &__copy__);
        FEPDeltas_exposer.def( "__deepcopy__", &__copy__);
        FEPDeltas_exposer.def( "clone", &__copy__);
        FEPDeltas_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Soiree::FEPDeltas >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FEPDeltas_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Soiree::FEPDeltas >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FEPDeltas_exposer.def( "__str__", &__str__< ::Soiree::FEPDeltas > );
        FEPDeltas_exposer.def( "__repr__", &__str__< ::Soiree::FEPDeltas > );
    }

}
