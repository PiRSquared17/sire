// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "SystemWideDeleter.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "moldeleter.h"

#include "uniformsampler.h"

#include "moldeleter.h"

SireMove::SystemWideDeleter __copy__(const SireMove::SystemWideDeleter &other){ return SireMove::SystemWideDeleter(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SystemWideDeleter_class(){

    { //::SireMove::SystemWideDeleter
        typedef bp::class_< SireMove::SystemWideDeleter, bp::bases< SireMove::MolDeleter, SireBase::Property > > SystemWideDeleter_exposer_t;
        SystemWideDeleter_exposer_t SystemWideDeleter_exposer = SystemWideDeleter_exposer_t( "SystemWideDeleter" );
        bp::scope SystemWideDeleter_scope( SystemWideDeleter_exposer );
        SystemWideDeleter_exposer.def( bp::init< >() );
        SystemWideDeleter_exposer.def( bp::init< SireMol::MoleculeGroup const & >(( bp::arg("molgroup") )) );
        SystemWideDeleter_exposer.def( bp::init< SireMove::Sampler const & >(( bp::arg("sampler") )) );
        SystemWideDeleter_exposer.def( bp::init< SireMove::SystemWideDeleter const & >(( bp::arg("other") )) );
        { //::SireMove::SystemWideDeleter::deleteFrom
        
            typedef ::boost::tuples::tuple< SireMol::Molecule, double, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireMove::SystemWideDeleter::*deleteFrom_function_type )( ::SireSystem::System & ) ;
            deleteFrom_function_type deleteFrom_function_value( &::SireMove::SystemWideDeleter::deleteFrom );
            
            SystemWideDeleter_exposer.def( 
                "deleteFrom"
                , deleteFrom_function_value
                , ( bp::arg("system") ) );
        
        }
        { //::SireMove::SystemWideDeleter::generator
        
            typedef ::SireMaths::RanGenerator const & ( ::SireMove::SystemWideDeleter::*generator_function_type )(  ) const;
            generator_function_type generator_function_value( &::SireMove::SystemWideDeleter::generator );
            
            SystemWideDeleter_exposer.def( 
                "generator"
                , generator_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::SystemWideDeleter::group
        
            typedef ::SireMol::MoleculeGroup const & ( ::SireMove::SystemWideDeleter::*group_function_type )(  ) const;
            group_function_type group_function_value( &::SireMove::SystemWideDeleter::group );
            
            SystemWideDeleter_exposer.def( 
                "group"
                , group_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        SystemWideDeleter_exposer.def( bp::self != bp::self );
        { //::SireMove::SystemWideDeleter::operator=
        
            typedef ::SireMove::SystemWideDeleter & ( ::SireMove::SystemWideDeleter::*assign_function_type )( ::SireMove::SystemWideDeleter const & ) ;
            assign_function_type assign_function_value( &::SireMove::SystemWideDeleter::operator= );
            
            SystemWideDeleter_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        SystemWideDeleter_exposer.def( bp::self == bp::self );
        { //::SireMove::SystemWideDeleter::sampler
        
            typedef ::SireMove::Sampler const & ( ::SireMove::SystemWideDeleter::*sampler_function_type )(  ) const;
            sampler_function_type sampler_function_value( &::SireMove::SystemWideDeleter::sampler );
            
            SystemWideDeleter_exposer.def( 
                "sampler"
                , sampler_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMove::SystemWideDeleter::setGenerator
        
            typedef void ( ::SireMove::SystemWideDeleter::*setGenerator_function_type )( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::SystemWideDeleter::setGenerator );
            
            SystemWideDeleter_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("generator") ) );
        
        }
        { //::SireMove::SystemWideDeleter::setGroup
        
            typedef void ( ::SireMove::SystemWideDeleter::*setGroup_function_type )( ::SireMol::MoleculeGroup const & ) ;
            setGroup_function_type setGroup_function_value( &::SireMove::SystemWideDeleter::setGroup );
            
            SystemWideDeleter_exposer.def( 
                "setGroup"
                , setGroup_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireMove::SystemWideDeleter::setSampler
        
            typedef void ( ::SireMove::SystemWideDeleter::*setSampler_function_type )( ::SireMove::Sampler const & ) ;
            setSampler_function_type setSampler_function_value( &::SireMove::SystemWideDeleter::setSampler );
            
            SystemWideDeleter_exposer.def( 
                "setSampler"
                , setSampler_function_value
                , ( bp::arg("sampler") ) );
        
        }
        { //::SireMove::SystemWideDeleter::setSampler
        
            typedef void ( ::SireMove::SystemWideDeleter::*setSampler_function_type )( ::SireMol::MoleculeGroup const & ) ;
            setSampler_function_type setSampler_function_value( &::SireMove::SystemWideDeleter::setSampler );
            
            SystemWideDeleter_exposer.def( 
                "setSampler"
                , setSampler_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireMove::SystemWideDeleter::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::SystemWideDeleter::typeName );
            
            SystemWideDeleter_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        SystemWideDeleter_exposer.staticmethod( "typeName" );
        SystemWideDeleter_exposer.def( "__copy__", &__copy__);
        SystemWideDeleter_exposer.def( "__deepcopy__", &__copy__);
        SystemWideDeleter_exposer.def( "clone", &__copy__);
        SystemWideDeleter_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::SystemWideDeleter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SystemWideDeleter_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::SystemWideDeleter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SystemWideDeleter_exposer.def( "__str__", &__str__< ::SireMove::SystemWideDeleter > );
        SystemWideDeleter_exposer.def( "__repr__", &__str__< ::SireMove::SystemWideDeleter > );
    }

}