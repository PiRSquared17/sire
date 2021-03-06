// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "MolInserter.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/quaternion.h"

#include "SireMol/molecule.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/units.h"

#include "SireVol/space.h"

#include "molinserter.h"

#include "molinserter.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_MolInserter_class(){

    { //::SireMove::MolInserter
        typedef bp::class_< SireMove::MolInserter, bp::bases< SireBase::Property >, boost::noncopyable > MolInserter_exposer_t;
        MolInserter_exposer_t MolInserter_exposer = MolInserter_exposer_t( "MolInserter", bp::no_init );
        bp::scope MolInserter_scope( MolInserter_exposer );
        { //::SireMove::MolInserter::generator
        
            typedef ::SireMaths::RanGenerator const & ( ::SireMove::MolInserter::*generator_function_type )(  ) const;
            generator_function_type generator_function_value( &::SireMove::MolInserter::generator );
            
            MolInserter_exposer.def( 
                "generator"
                , generator_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::MolInserter::groups
        
            typedef ::SireMol::MGIDsAndMaps const & ( ::SireMove::MolInserter::*groups_function_type )(  ) const;
            groups_function_type groups_function_value( &::SireMove::MolInserter::groups );
            
            MolInserter_exposer.def( 
                "groups"
                , groups_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::MolInserter::insert
        
            typedef double ( ::SireMove::MolInserter::*insert_function_type )( ::SireMol::Molecule const &,::SireSystem::System &,::SireVol::Space const & ) ;
            insert_function_type insert_function_value( &::SireMove::MolInserter::insert );
            
            MolInserter_exposer.def( 
                "insert"
                , insert_function_value
                , ( bp::arg("molecule"), bp::arg("system"), bp::arg("space") ) );
        
        }
        { //::SireMove::MolInserter::insert
        
            typedef double ( ::SireMove::MolInserter::*insert_function_type )( ::SireMol::PartialMolecule const &,::SireSystem::System &,::SireVol::Space const & ) ;
            insert_function_type insert_function_value( &::SireMove::MolInserter::insert );
            
            MolInserter_exposer.def( 
                "insert"
                , insert_function_value
                , ( bp::arg("molecule"), bp::arg("system"), bp::arg("space") ) );
        
        }
        { //::SireMove::MolInserter::null
        
            typedef ::SireMove::NullInserter const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMove::MolInserter::null );
            
            MolInserter_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::MolInserter::setGenerator
        
            typedef void ( ::SireMove::MolInserter::*setGenerator_function_type )( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::MolInserter::setGenerator );
            
            MolInserter_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("generator") ) );
        
        }
        { //::SireMove::MolInserter::setGroups
        
            typedef void ( ::SireMove::MolInserter::*setGroups_function_type )( ::SireMol::MGIDsAndMaps const & ) ;
            setGroups_function_type setGroups_function_value( &::SireMove::MolInserter::setGroups );
            
            MolInserter_exposer.def( 
                "setGroups"
                , setGroups_function_value
                , ( bp::arg("mgids") ) );
        
        }
        { //::SireMove::MolInserter::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::MolInserter::typeName );
            
            MolInserter_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        MolInserter_exposer.staticmethod( "null" );
        MolInserter_exposer.staticmethod( "typeName" );
        MolInserter_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::MolInserter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolInserter_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::MolInserter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolInserter_exposer.def( "__str__", &__str__< ::SireMove::MolInserter > );
        MolInserter_exposer.def( "__repr__", &__str__< ::SireMove::MolInserter > );
    }

}
