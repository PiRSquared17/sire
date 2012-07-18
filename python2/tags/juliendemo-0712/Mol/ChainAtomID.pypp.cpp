// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ChainAtomID.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "atomidentifier.h"

#include "cgidentifier.h"

#include "chainidentifier.h"

#include "groupatomids.h"

#include "residentifier.h"

#include "segidentifier.h"

#include "groupatomids.h"

SireMol::GroupAtomID<SireMol::ChainID, SireMol::AtomID> __copy__(const SireMol::GroupAtomID<SireMol::ChainID, SireMol::AtomID> &other){ return SireMol::GroupAtomID<SireMol::ChainID, SireMol::AtomID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ChainAtomID_class(){

    { //::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >
        typedef bp::class_< SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >, bp::bases< SireMol::GroupAtomIDBase, SireMol::AtomID, SireID::ID > > ChainAtomID_exposer_t;
        ChainAtomID_exposer_t ChainAtomID_exposer = ChainAtomID_exposer_t( "ChainAtomID" );
        bp::scope ChainAtomID_scope( ChainAtomID_exposer );
        ChainAtomID_exposer.def( bp::init< >() );
        ChainAtomID_exposer.def( bp::init< SireMol::ChainID const &, SireMol::AtomID const & >(( bp::arg("group"), bp::arg("atom") )) );
        ChainAtomID_exposer.def( bp::init< SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID > const & >(( bp::arg("other") )) );
        { //::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::hash
        
            typedef SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID > exported_class_t;
            typedef ::uint ( ::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::*hash_function_type )(  ) const;
            hash_function_type hash_function_value( &::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::hash );
            
            ChainAtomID_exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::isNull
        
            typedef SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID > exported_class_t;
            typedef bool ( ::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::isNull );
            
            ChainAtomID_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::map
        
            typedef SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID > exported_class_t;
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::*map_function_type )( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::map );
            
            ChainAtomID_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        ChainAtomID_exposer.def( bp::self != bp::self );
        ChainAtomID_exposer.def( bp::self == bp::self );
        ChainAtomID_exposer.def( bp::self == bp::other< SireID::ID >() );
        { //::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::toString
        
            typedef SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID > exported_class_t;
            typedef ::QString ( ::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::toString );
            
            ChainAtomID_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::typeName
        
            typedef SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::typeName );
            
            ChainAtomID_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::what
        
            typedef SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID > exported_class_t;
            typedef char const * ( ::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMol::GroupAtomID< SireMol::ChainID, SireMol::AtomID >::what );
            
            ChainAtomID_exposer.def( 
                "what"
                , what_function_value );
        
        }
        ChainAtomID_exposer.staticmethod( "typeName" );
        ChainAtomID_exposer.def( "__copy__", &__copy__);
        ChainAtomID_exposer.def( "__deepcopy__", &__copy__);
        ChainAtomID_exposer.def( "clone", &__copy__);
        ChainAtomID_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::GroupAtomID<SireMol::ChainID, SireMol::AtomID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainAtomID_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::GroupAtomID<SireMol::ChainID, SireMol::AtomID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainAtomID_exposer.def( "__str__", &__str__< ::SireMol::GroupAtomID<SireMol::ChainID, SireMol::AtomID> > );
        ChainAtomID_exposer.def( "__repr__", &__str__< ::SireMol::GroupAtomID<SireMol::ChainID, SireMol::AtomID> > );
    }

}
