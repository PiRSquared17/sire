// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Sampler.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/molecule.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "sampler.h"

#include "uniformsampler.h"

#include "sampler.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Sampler_class(){

    { //::SireMove::Sampler
        typedef bp::class_< SireMove::Sampler, bp::bases< SireBase::Property >, boost::noncopyable > Sampler_exposer_t;
        Sampler_exposer_t Sampler_exposer = Sampler_exposer_t( "Sampler", bp::no_init );
        bp::scope Sampler_scope( Sampler_exposer );
        { //::SireMove::Sampler::generator
        
            typedef ::SireMaths::RanGenerator const & ( ::SireMove::Sampler::*generator_function_type )(  ) const;
            generator_function_type generator_function_value( &::SireMove::Sampler::generator );
            
            Sampler_exposer.def( 
                "generator"
                , generator_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::Sampler::group
        
            typedef ::SireMol::MoleculeGroup const & ( ::SireMove::Sampler::*group_function_type )(  ) const;
            group_function_type group_function_value( &::SireMove::Sampler::group );
            
            Sampler_exposer.def( 
                "group"
                , group_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMove::Sampler::isBiased
        
            typedef bool ( ::SireMove::Sampler::*isBiased_function_type )(  ) const;
            isBiased_function_type isBiased_function_value( &::SireMove::Sampler::isBiased );
            
            Sampler_exposer.def( 
                "isBiased"
                , isBiased_function_value );
        
        }
        { //::SireMove::Sampler::null
        
            typedef ::SireMove::UniformSampler const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMove::Sampler::null );
            
            Sampler_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMove::Sampler::probabilityOf
        
            typedef double ( ::SireMove::Sampler::*probabilityOf_function_type )( ::SireMol::PartialMolecule const & ) const;
            probabilityOf_function_type probabilityOf_function_value( &::SireMove::Sampler::probabilityOf );
            
            Sampler_exposer.def( 
                "probabilityOf"
                , probabilityOf_function_value
                , ( bp::arg("molecule") ) );
        
        }
        { //::SireMove::Sampler::probabilityOfMolecule
        
            typedef double ( ::SireMove::Sampler::*probabilityOfMolecule_function_type )( ::SireMol::Molecule const & ) const;
            probabilityOfMolecule_function_type probabilityOfMolecule_function_value( &::SireMove::Sampler::probabilityOfMolecule );
            
            Sampler_exposer.def( 
                "probabilityOfMolecule"
                , probabilityOfMolecule_function_value
                , ( bp::arg("molecule") ) );
        
        }
        { //::SireMove::Sampler::sample
        
            typedef ::boost::tuples::tuple< SireMol::PartialMolecule, double, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireMove::Sampler::*sample_function_type )(  ) const;
            sample_function_type sample_function_value( &::SireMove::Sampler::sample );
            
            Sampler_exposer.def( 
                "sample"
                , sample_function_value );
        
        }
        { //::SireMove::Sampler::sampleMolecule
        
            typedef ::boost::tuples::tuple< SireMol::Molecule, double, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireMove::Sampler::*sampleMolecule_function_type )(  ) const;
            sampleMolecule_function_type sampleMolecule_function_value( &::SireMove::Sampler::sampleMolecule );
            
            Sampler_exposer.def( 
                "sampleMolecule"
                , sampleMolecule_function_value );
        
        }
        { //::SireMove::Sampler::setGenerator
        
            typedef void ( ::SireMove::Sampler::*setGenerator_function_type )( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::Sampler::setGenerator );
            
            Sampler_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("generator") ) );
        
        }
        { //::SireMove::Sampler::setGroup
        
            typedef void ( ::SireMove::Sampler::*setGroup_function_type )( ::SireMol::MoleculeGroup const & ) ;
            setGroup_function_type setGroup_function_value( &::SireMove::Sampler::setGroup );
            
            Sampler_exposer.def( 
                "setGroup"
                , setGroup_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireMove::Sampler::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::Sampler::typeName );
            
            Sampler_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMove::Sampler::updateFrom
        
            typedef void ( ::SireMove::Sampler::*updateFrom_function_type )( ::SireSystem::System const & ) ;
            updateFrom_function_type updateFrom_function_value( &::SireMove::Sampler::updateFrom );
            
            Sampler_exposer.def( 
                "updateFrom"
                , updateFrom_function_value
                , ( bp::arg("system") ) );
        
        }
        Sampler_exposer.staticmethod( "null" );
        Sampler_exposer.staticmethod( "typeName" );
        Sampler_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::Sampler >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Sampler_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::Sampler >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Sampler_exposer.def( "__str__", &__str__< ::SireMove::Sampler > );
        Sampler_exposer.def( "__repr__", &__str__< ::SireMove::Sampler > );
    }

}
