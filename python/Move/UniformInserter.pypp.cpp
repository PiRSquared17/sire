// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "UniformInserter.pypp.hpp"

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

SireMove::UniformInserter __copy__(const SireMove::UniformInserter &other){ return SireMove::UniformInserter(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_UniformInserter_class(){

    { //::SireMove::UniformInserter
        typedef bp::class_< SireMove::UniformInserter, bp::bases< SireMove::MolInserter, SireBase::Property > > UniformInserter_exposer_t;
        UniformInserter_exposer_t UniformInserter_exposer = UniformInserter_exposer_t( "UniformInserter", bp::init< >() );
        bp::scope UniformInserter_scope( UniformInserter_exposer );
        UniformInserter_exposer.def( bp::init< SireMol::MGIDsAndMaps const & >(( bp::arg("mgids") )) );
        UniformInserter_exposer.def( bp::init< SireMove::UniformInserter const & >(( bp::arg("other") )) );
        { //::SireMove::UniformInserter::insert
        
            typedef double ( ::SireMove::UniformInserter::*insert_function_type )( ::SireMol::Molecule const &,::SireSystem::System &,::SireVol::Space const & ) ;
            insert_function_type insert_function_value( &::SireMove::UniformInserter::insert );
            
            UniformInserter_exposer.def( 
                "insert"
                , insert_function_value
                , ( bp::arg("molecule"), bp::arg("system"), bp::arg("space") ) );
        
        }
        { //::SireMove::UniformInserter::insert
        
            typedef double ( ::SireMove::UniformInserter::*insert_function_type )( ::SireMol::PartialMolecule const &,::SireSystem::System &,::SireVol::Space const & ) ;
            insert_function_type insert_function_value( &::SireMove::UniformInserter::insert );
            
            UniformInserter_exposer.def( 
                "insert"
                , insert_function_value
                , ( bp::arg("molecule"), bp::arg("system"), bp::arg("space") ) );
        
        }
        UniformInserter_exposer.def( bp::self != bp::self );
        { //::SireMove::UniformInserter::operator=
        
            typedef ::SireMove::UniformInserter & ( ::SireMove::UniformInserter::*assign_function_type )( ::SireMove::UniformInserter const & ) ;
            assign_function_type assign_function_value( &::SireMove::UniformInserter::operator= );
            
            UniformInserter_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        UniformInserter_exposer.def( bp::self == bp::self );
        { //::SireMove::UniformInserter::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::UniformInserter::typeName );
            
            UniformInserter_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        UniformInserter_exposer.staticmethod( "typeName" );
        UniformInserter_exposer.def( "__copy__", &__copy__);
        UniformInserter_exposer.def( "__deepcopy__", &__copy__);
        UniformInserter_exposer.def( "clone", &__copy__);
        UniformInserter_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::UniformInserter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        UniformInserter_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::UniformInserter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        UniformInserter_exposer.def( "__str__", &__str__< ::SireMove::UniformInserter > );
        UniformInserter_exposer.def( "__repr__", &__str__< ::SireMove::UniformInserter > );
    }

}
