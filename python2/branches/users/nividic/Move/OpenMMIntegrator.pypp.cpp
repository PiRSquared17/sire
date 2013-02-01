// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "OpenMMIntegrator.pypp.hpp"

namespace bp = boost::python;

#include "SireFF/forcetable.h"

#include "SireIO/amber.h"

#include "SireMM/atomljs.h"

#include "SireMM/internalff.h"

#include "SireMaths/rangenerator.h"

#include "SireMaths/vector.h"

#include "SireMol/amberparameters.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atommasses.h"

#include "SireMol/bondid.h"

#include "SireMol/connectivity.h"

#include "SireMol/mgname.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleditor.h"

#include "SireMol/partialmolecule.h"

#include "SireMove/flexibility.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/convert.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "ensemble.h"

#include "fastio.h"

#include "openmmintegrator.h"

#include <QDebug>

#include <QElapsedTimer>

#include <iomanip>

#include <iostream>

#include <queue>

#include "openmmintegrator.h"

SireMove::OpenMMIntegrator __copy__(const SireMove::OpenMMIntegrator &other){ return SireMove::OpenMMIntegrator(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_OpenMMIntegrator_class(){

    { //::SireMove::OpenMMIntegrator
        typedef bp::class_< SireMove::OpenMMIntegrator, bp::bases< SireMove::Integrator, SireBase::Property > > OpenMMIntegrator_exposer_t;
        OpenMMIntegrator_exposer_t OpenMMIntegrator_exposer = OpenMMIntegrator_exposer_t( "OpenMMIntegrator", bp::init< bp::optional< bool > >(( bp::arg("frequent_save_velocities")=(bool)(false) )) );
        bp::scope OpenMMIntegrator_scope( OpenMMIntegrator_exposer );
        OpenMMIntegrator_exposer.def( bp::init< SireMove::OpenMMIntegrator const & >(( bp::arg("other") )) );
        { //::SireMove::OpenMMIntegrator::createWorkspace
        
            typedef ::SireMove::IntegratorWorkspacePtr ( ::SireMove::OpenMMIntegrator::*createWorkspace_function_type )( ::SireBase::PropertyMap const & ) const;
            createWorkspace_function_type createWorkspace_function_value( &::SireMove::OpenMMIntegrator::createWorkspace );
            
            OpenMMIntegrator_exposer.def( 
                "createWorkspace"
                , createWorkspace_function_value
                , ( bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        { //::SireMove::OpenMMIntegrator::createWorkspace
        
            typedef ::SireMove::IntegratorWorkspacePtr ( ::SireMove::OpenMMIntegrator::*createWorkspace_function_type )( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) const;
            createWorkspace_function_type createWorkspace_function_value( &::SireMove::OpenMMIntegrator::createWorkspace );
            
            OpenMMIntegrator_exposer.def( 
                "createWorkspace"
                , createWorkspace_function_value
                , ( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        { //::SireMove::OpenMMIntegrator::ensemble
        
            typedef ::SireMove::Ensemble ( ::SireMove::OpenMMIntegrator::*ensemble_function_type )(  ) const;
            ensemble_function_type ensemble_function_value( &::SireMove::OpenMMIntegrator::ensemble );
            
            OpenMMIntegrator_exposer.def( 
                "ensemble"
                , ensemble_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getAlchemical_values
        
            typedef ::QVector< double > ( ::SireMove::OpenMMIntegrator::*getAlchemical_values_function_type )(  ) ;
            getAlchemical_values_function_type getAlchemical_values_function_value( &::SireMove::OpenMMIntegrator::getAlchemical_values );
            
            OpenMMIntegrator_exposer.def( 
                "getAlchemical_values"
                , getAlchemical_values_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getAndersen
        
            typedef bool ( ::SireMove::OpenMMIntegrator::*getAndersen_function_type )(  ) ;
            getAndersen_function_type getAndersen_function_value( &::SireMove::OpenMMIntegrator::getAndersen );
            
            OpenMMIntegrator_exposer.def( 
                "getAndersen"
                , getAndersen_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getAndersen_frequency
        
            typedef double ( ::SireMove::OpenMMIntegrator::*getAndersen_frequency_function_type )(  ) ;
            getAndersen_frequency_function_type getAndersen_frequency_function_value( &::SireMove::OpenMMIntegrator::getAndersen_frequency );
            
            OpenMMIntegrator_exposer.def( 
                "getAndersen_frequency"
                , getAndersen_frequency_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getCMMremoval_frequency
        
            typedef int ( ::SireMove::OpenMMIntegrator::*getCMMremoval_frequency_function_type )(  ) ;
            getCMMremoval_frequency_function_type getCMMremoval_frequency_function_value( &::SireMove::OpenMMIntegrator::getCMMremoval_frequency );
            
            OpenMMIntegrator_exposer.def( 
                "getCMMremoval_frequency"
                , getCMMremoval_frequency_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getConstraintType
        
            typedef ::QString ( ::SireMove::OpenMMIntegrator::*getConstraintType_function_type )(  ) ;
            getConstraintType_function_type getConstraintType_function_value( &::SireMove::OpenMMIntegrator::getConstraintType );
            
            OpenMMIntegrator_exposer.def( 
                "getConstraintType"
                , getConstraintType_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getCutoffType
        
            typedef ::QString ( ::SireMove::OpenMMIntegrator::*getCutoffType_function_type )(  ) ;
            getCutoffType_function_type getCutoffType_function_value( &::SireMove::OpenMMIntegrator::getCutoffType );
            
            OpenMMIntegrator_exposer.def( 
                "getCutoffType"
                , getCutoffType_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getCutoff_distance
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::OpenMMIntegrator::*getCutoff_distance_function_type )(  ) ;
            getCutoff_distance_function_type getCutoff_distance_function_value( &::SireMove::OpenMMIntegrator::getCutoff_distance );
            
            OpenMMIntegrator_exposer.def( 
                "getCutoff_distance"
                , getCutoff_distance_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getField_dielectric
        
            typedef double ( ::SireMove::OpenMMIntegrator::*getField_dielectric_function_type )(  ) ;
            getField_dielectric_function_type getField_dielectric_function_value( &::SireMove::OpenMMIntegrator::getField_dielectric );
            
            OpenMMIntegrator_exposer.def( 
                "getField_dielectric"
                , getField_dielectric_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getMCBarostat
        
            typedef bool ( ::SireMove::OpenMMIntegrator::*getMCBarostat_function_type )(  ) ;
            getMCBarostat_function_type getMCBarostat_function_value( &::SireMove::OpenMMIntegrator::getMCBarostat );
            
            OpenMMIntegrator_exposer.def( 
                "getMCBarostat"
                , getMCBarostat_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getMCBarostat_frequency
        
            typedef int ( ::SireMove::OpenMMIntegrator::*getMCBarostat_frequency_function_type )(  ) ;
            getMCBarostat_frequency_function_type getMCBarostat_frequency_function_value( &::SireMove::OpenMMIntegrator::getMCBarostat_frequency );
            
            OpenMMIntegrator_exposer.def( 
                "getMCBarostat_frequency"
                , getMCBarostat_frequency_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getPlatform
        
            typedef ::QString ( ::SireMove::OpenMMIntegrator::*getPlatform_function_type )(  ) ;
            getPlatform_function_type getPlatform_function_value( &::SireMove::OpenMMIntegrator::getPlatform );
            
            OpenMMIntegrator_exposer.def( 
                "getPlatform"
                , getPlatform_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getPressure
        
            typedef ::SireUnits::Dimension::Pressure ( ::SireMove::OpenMMIntegrator::*getPressure_function_type )(  ) ;
            getPressure_function_type getPressure_function_value( &::SireMove::OpenMMIntegrator::getPressure );
            
            OpenMMIntegrator_exposer.def( 
                "getPressure"
                , getPressure_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getRestraint
        
            typedef bool ( ::SireMove::OpenMMIntegrator::*getRestraint_function_type )(  ) ;
            getRestraint_function_type getRestraint_function_value( &::SireMove::OpenMMIntegrator::getRestraint );
            
            OpenMMIntegrator_exposer.def( 
                "getRestraint"
                , getRestraint_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::getTemperature
        
            typedef ::SireUnits::Dimension::Temperature ( ::SireMove::OpenMMIntegrator::*getTemperature_function_type )(  ) ;
            getTemperature_function_type getTemperature_function_value( &::SireMove::OpenMMIntegrator::getTemperature );
            
            OpenMMIntegrator_exposer.def( 
                "getTemperature"
                , getTemperature_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::integrate
        
            typedef void ( ::SireMove::OpenMMIntegrator::*integrate_function_type )( ::SireMove::IntegratorWorkspace &,::SireCAS::Symbol const &,::SireUnits::Dimension::Time,int,bool ) const;
            integrate_function_type integrate_function_value( &::SireMove::OpenMMIntegrator::integrate );
            
            OpenMMIntegrator_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("workspace"), bp::arg("nrg_component"), bp::arg("timestep"), bp::arg("nmoves"), bp::arg("record_stats") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::isTimeReversible
        
            typedef bool ( ::SireMove::OpenMMIntegrator::*isTimeReversible_function_type )(  ) const;
            isTimeReversible_function_type isTimeReversible_function_value( &::SireMove::OpenMMIntegrator::isTimeReversible );
            
            OpenMMIntegrator_exposer.def( 
                "isTimeReversible"
                , isTimeReversible_function_value );
        
        }
        OpenMMIntegrator_exposer.def( bp::self != bp::self );
        { //::SireMove::OpenMMIntegrator::operator=
        
            typedef ::SireMove::OpenMMIntegrator & ( ::SireMove::OpenMMIntegrator::*assign_function_type )( ::SireMove::OpenMMIntegrator const & ) ;
            assign_function_type assign_function_value( &::SireMove::OpenMMIntegrator::operator= );
            
            OpenMMIntegrator_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        OpenMMIntegrator_exposer.def( bp::self == bp::self );
        { //::SireMove::OpenMMIntegrator::setAlchemical_values
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setAlchemical_values_function_type )( ::QVector< double > ) ;
            setAlchemical_values_function_type setAlchemical_values_function_value( &::SireMove::OpenMMIntegrator::setAlchemical_values );
            
            OpenMMIntegrator_exposer.def( 
                "setAlchemical_values"
                , setAlchemical_values_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setAndersen
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setAndersen_function_type )( bool ) ;
            setAndersen_function_type setAndersen_function_value( &::SireMove::OpenMMIntegrator::setAndersen );
            
            OpenMMIntegrator_exposer.def( 
                "setAndersen"
                , setAndersen_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setAndersen_frequency
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setAndersen_frequency_function_type )( double ) ;
            setAndersen_frequency_function_type setAndersen_frequency_function_value( &::SireMove::OpenMMIntegrator::setAndersen_frequency );
            
            OpenMMIntegrator_exposer.def( 
                "setAndersen_frequency"
                , setAndersen_frequency_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setCMMremoval_frequency
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setCMMremoval_frequency_function_type )( int ) ;
            setCMMremoval_frequency_function_type setCMMremoval_frequency_function_value( &::SireMove::OpenMMIntegrator::setCMMremoval_frequency );
            
            OpenMMIntegrator_exposer.def( 
                "setCMMremoval_frequency"
                , setCMMremoval_frequency_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setConstraintType
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setConstraintType_function_type )( ::QString ) ;
            setConstraintType_function_type setConstraintType_function_value( &::SireMove::OpenMMIntegrator::setConstraintType );
            
            OpenMMIntegrator_exposer.def( 
                "setConstraintType"
                , setConstraintType_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setCutoffType
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setCutoffType_function_type )( ::QString ) ;
            setCutoffType_function_type setCutoffType_function_value( &::SireMove::OpenMMIntegrator::setCutoffType );
            
            OpenMMIntegrator_exposer.def( 
                "setCutoffType"
                , setCutoffType_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setCutoff_distance
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setCutoff_distance_function_type )( ::SireUnits::Dimension::Length ) ;
            setCutoff_distance_function_type setCutoff_distance_function_value( &::SireMove::OpenMMIntegrator::setCutoff_distance );
            
            OpenMMIntegrator_exposer.def( 
                "setCutoff_distance"
                , setCutoff_distance_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setField_dielectric
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setField_dielectric_function_type )( double ) ;
            setField_dielectric_function_type setField_dielectric_function_value( &::SireMove::OpenMMIntegrator::setField_dielectric );
            
            OpenMMIntegrator_exposer.def( 
                "setField_dielectric"
                , setField_dielectric_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setMCBarostat
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setMCBarostat_function_type )( bool ) ;
            setMCBarostat_function_type setMCBarostat_function_value( &::SireMove::OpenMMIntegrator::setMCBarostat );
            
            OpenMMIntegrator_exposer.def( 
                "setMCBarostat"
                , setMCBarostat_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setMCBarostat_frequency
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setMCBarostat_frequency_function_type )( int ) ;
            setMCBarostat_frequency_function_type setMCBarostat_frequency_function_value( &::SireMove::OpenMMIntegrator::setMCBarostat_frequency );
            
            OpenMMIntegrator_exposer.def( 
                "setMCBarostat_frequency"
                , setMCBarostat_frequency_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setPlatform
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setPlatform_function_type )( ::QString ) ;
            setPlatform_function_type setPlatform_function_value( &::SireMove::OpenMMIntegrator::setPlatform );
            
            OpenMMIntegrator_exposer.def( 
                "setPlatform"
                , setPlatform_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setPressure
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setPressure_function_type )( ::SireUnits::Dimension::Pressure ) ;
            setPressure_function_type setPressure_function_value( &::SireMove::OpenMMIntegrator::setPressure );
            
            OpenMMIntegrator_exposer.def( 
                "setPressure"
                , setPressure_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setRestraint
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setRestraint_function_type )( bool ) ;
            setRestraint_function_type setRestraint_function_value( &::SireMove::OpenMMIntegrator::setRestraint );
            
            OpenMMIntegrator_exposer.def( 
                "setRestraint"
                , setRestraint_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::setTemperature
        
            typedef void ( ::SireMove::OpenMMIntegrator::*setTemperature_function_type )( ::SireUnits::Dimension::Temperature ) ;
            setTemperature_function_type setTemperature_function_value( &::SireMove::OpenMMIntegrator::setTemperature );
            
            OpenMMIntegrator_exposer.def( 
                "setTemperature"
                , setTemperature_function_value
                , ( bp::arg("arg0") ) );
        
        }
        { //::SireMove::OpenMMIntegrator::toString
        
            typedef ::QString ( ::SireMove::OpenMMIntegrator::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMove::OpenMMIntegrator::toString );
            
            OpenMMIntegrator_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMove::OpenMMIntegrator::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::OpenMMIntegrator::typeName );
            
            OpenMMIntegrator_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        OpenMMIntegrator_exposer.staticmethod( "typeName" );
        OpenMMIntegrator_exposer.def( "__copy__", &__copy__);
        OpenMMIntegrator_exposer.def( "__deepcopy__", &__copy__);
        OpenMMIntegrator_exposer.def( "clone", &__copy__);
        OpenMMIntegrator_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::OpenMMIntegrator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        OpenMMIntegrator_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::OpenMMIntegrator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        OpenMMIntegrator_exposer.def( "__str__", &__str__< ::SireMove::OpenMMIntegrator > );
        OpenMMIntegrator_exposer.def( "__repr__", &__str__< ::SireMove::OpenMMIntegrator > );
    }

}
