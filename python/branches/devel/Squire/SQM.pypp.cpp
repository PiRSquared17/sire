// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SQM.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/findexe.h"

#include "SireBase/process.h"

#include "SireBase/tempdir.h"

#include "SireError/errors.h"

#include "SireFF/potentialtable.h"

#include "SireMM/cljprobe.h"

#include "SireMaths/vector.h"

#include "SireMol/element.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "SireVol/grid.h"

#include "latticecharges.h"

#include "qmpotential.h"

#include "sqm.h"

#include "tostring.h"

#include <QDebug>

#include <QTime>

#include "sqm.h"

Squire::SQM __copy__(const Squire::SQM &other){ return Squire::SQM(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SQM_class(){

    { //::Squire::SQM
        typedef bp::class_< Squire::SQM, bp::bases< Squire::QMProgram, SireBase::Property > > SQM_exposer_t;
        SQM_exposer_t SQM_exposer = SQM_exposer_t( "SQM", bp::init< >() );
        bp::scope SQM_scope( SQM_exposer );
        SQM_exposer.def( bp::init< QString const & >(( bp::arg("SQM") )) );
        SQM_exposer.def( bp::init< Squire::SQM const & >(( bp::arg("other") )) );
        { //::Squire::SQM::energyTemplate
        
            typedef ::QString const & ( ::Squire::SQM::*energyTemplate_function_type )(  ) const;
            energyTemplate_function_type energyTemplate_function_value( &::Squire::SQM::energyTemplate );
            
            SQM_exposer.def( 
                "energyTemplate"
                , energyTemplate_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::SQM::environment
        
            typedef ::QHash< QString, QString > const & ( ::Squire::SQM::*environment_function_type )(  ) const;
            environment_function_type environment_function_value( &::Squire::SQM::environment );
            
            SQM_exposer.def( 
                "environment"
                , environment_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::SQM::environment
        
            typedef ::QString ( ::Squire::SQM::*environment_function_type )( ::QString const & ) const;
            environment_function_type environment_function_value( &::Squire::SQM::environment );
            
            SQM_exposer.def( 
                "environment"
                , environment_function_value
                , ( bp::arg("variable") ) );
        
        }
        { //::Squire::SQM::executable
        
            typedef ::QString ( ::Squire::SQM::*executable_function_type )(  ) const;
            executable_function_type executable_function_value( &::Squire::SQM::executable );
            
            SQM_exposer.def( 
                "executable"
                , executable_function_value );
        
        }
        { //::Squire::SQM::forceTemplate
        
            typedef ::QString const & ( ::Squire::SQM::*forceTemplate_function_type )(  ) const;
            forceTemplate_function_type forceTemplate_function_value( &::Squire::SQM::forceTemplate );
            
            SQM_exposer.def( 
                "forceTemplate"
                , forceTemplate_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::SQM::maximumRunTime
        
            typedef int ( ::Squire::SQM::*maximumRunTime_function_type )(  ) const;
            maximumRunTime_function_type maximumRunTime_function_value( &::Squire::SQM::maximumRunTime );
            
            SQM_exposer.def( 
                "maximumRunTime"
                , maximumRunTime_function_value );
        
        }
        { //::Squire::SQM::method
        
            typedef ::QString const & ( ::Squire::SQM::*method_function_type )(  ) const;
            method_function_type method_function_value( &::Squire::SQM::method );
            
            SQM_exposer.def( 
                "method"
                , method_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        SQM_exposer.def( bp::self != bp::self );
        { //::Squire::SQM::operator=
        
            typedef ::Squire::SQM & ( ::Squire::SQM::*assign_function_type )( ::Squire::SQM const & ) ;
            assign_function_type assign_function_value( &::Squire::SQM::operator= );
            
            SQM_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        SQM_exposer.def( bp::self == bp::self );
        { //::Squire::SQM::setEnergyTemplate
        
            typedef void ( ::Squire::SQM::*setEnergyTemplate_function_type )( ::QString const & ) ;
            setEnergyTemplate_function_type setEnergyTemplate_function_value( &::Squire::SQM::setEnergyTemplate );
            
            SQM_exposer.def( 
                "setEnergyTemplate"
                , setEnergyTemplate_function_value
                , ( bp::arg("energy_template") ) );
        
        }
        { //::Squire::SQM::setEnvironment
        
            typedef void ( ::Squire::SQM::*setEnvironment_function_type )( ::QString const &,::QString const & ) ;
            setEnvironment_function_type setEnvironment_function_value( &::Squire::SQM::setEnvironment );
            
            SQM_exposer.def( 
                "setEnvironment"
                , setEnvironment_function_value
                , ( bp::arg("variable"), bp::arg("value") ) );
        
        }
        { //::Squire::SQM::setExecutable
        
            typedef void ( ::Squire::SQM::*setExecutable_function_type )( ::QString const & ) ;
            setExecutable_function_type setExecutable_function_value( &::Squire::SQM::setExecutable );
            
            SQM_exposer.def( 
                "setExecutable"
                , setExecutable_function_value
                , ( bp::arg("SQM_exe") ) );
        
        }
        { //::Squire::SQM::setForceTemplate
        
            typedef void ( ::Squire::SQM::*setForceTemplate_function_type )( ::QString const & ) ;
            setForceTemplate_function_type setForceTemplate_function_value( &::Squire::SQM::setForceTemplate );
            
            SQM_exposer.def( 
                "setForceTemplate"
                , setForceTemplate_function_value
                , ( bp::arg("force_template") ) );
        
        }
        { //::Squire::SQM::setMaximumRunTime
        
            typedef void ( ::Squire::SQM::*setMaximumRunTime_function_type )( int ) ;
            setMaximumRunTime_function_type setMaximumRunTime_function_value( &::Squire::SQM::setMaximumRunTime );
            
            SQM_exposer.def( 
                "setMaximumRunTime"
                , setMaximumRunTime_function_value
                , ( bp::arg("max_runtime") ) );
        
        }
        { //::Squire::SQM::setMethod
        
            typedef void ( ::Squire::SQM::*setMethod_function_type )( ::QString const & ) ;
            setMethod_function_type setMethod_function_value( &::Squire::SQM::setMethod );
            
            SQM_exposer.def( 
                "setMethod"
                , setMethod_function_value
                , ( bp::arg("method") ) );
        
        }
        { //::Squire::SQM::setTotalCharge
        
            typedef void ( ::Squire::SQM::*setTotalCharge_function_type )( int ) ;
            setTotalCharge_function_type setTotalCharge_function_value( &::Squire::SQM::setTotalCharge );
            
            SQM_exposer.def( 
                "setTotalCharge"
                , setTotalCharge_function_value
                , ( bp::arg("charge") ) );
        
        }
        { //::Squire::SQM::supportsLatticeCharges
        
            typedef bool ( ::Squire::SQM::*supportsLatticeCharges_function_type )(  ) const;
            supportsLatticeCharges_function_type supportsLatticeCharges_function_value( &::Squire::SQM::supportsLatticeCharges );
            
            SQM_exposer.def( 
                "supportsLatticeCharges"
                , supportsLatticeCharges_function_value );
        
        }
        { //::Squire::SQM::totalCharge
        
            typedef int ( ::Squire::SQM::*totalCharge_function_type )(  ) const;
            totalCharge_function_type totalCharge_function_value( &::Squire::SQM::totalCharge );
            
            SQM_exposer.def( 
                "totalCharge"
                , totalCharge_function_value );
        
        }
        { //::Squire::SQM::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::SQM::typeName );
            
            SQM_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        SQM_exposer.staticmethod( "typeName" );
        SQM_exposer.def( "__copy__", &__copy__);
        SQM_exposer.def( "__deepcopy__", &__copy__);
        SQM_exposer.def( "clone", &__copy__);
        SQM_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::SQM >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SQM_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::SQM >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SQM_exposer.def( "__str__", &__str__< ::Squire::SQM > );
        SQM_exposer.def( "__repr__", &__str__< ::Squire::SQM > );
    }

}
