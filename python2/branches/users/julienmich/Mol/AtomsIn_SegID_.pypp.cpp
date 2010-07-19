// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomsIn_SegID_.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "atom.h"

#include "chain.h"

#include "cutgroup.h"

#include "editor.hpp"

#include "groupatomids.h"

#include "groupgroupids.h"

#include "moleculegroup.h"

#include "moleculegroups.h"

#include "molecules.h"

#include "molinfo.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "residue.h"

#include "segid.h"

#include "segidentifier.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include "segid.h"

SireMol::AtomsIn<SireMol::SegID> __copy__(const SireMol::AtomsIn<SireMol::SegID> &other){ return SireMol::AtomsIn<SireMol::SegID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomsIn_SegID__class(){

    { //::SireMol::AtomsIn< SireMol::SegID >
        typedef bp::class_< SireMol::AtomsIn< SireMol::SegID >, bp::bases< SireMol::AtomID, SireID::ID > > AtomsIn_SegID__exposer_t;
        AtomsIn_SegID__exposer_t AtomsIn_SegID__exposer = AtomsIn_SegID__exposer_t( "AtomsIn_SegID_" );
        bp::scope AtomsIn_SegID__scope( AtomsIn_SegID__exposer );
        AtomsIn_SegID__exposer.def( bp::init< >() );
        AtomsIn_SegID__exposer.def( bp::init< SireMol::SegID const & >(( bp::arg("id") )) );
        AtomsIn_SegID__exposer.def( bp::init< SireMol::SegID const &, qint32 >(( bp::arg("id"), bp::arg("i") )) );
        AtomsIn_SegID__exposer.def( bp::init< SireMol::SegID const &, qint32, qint32 >(( bp::arg("id"), bp::arg("i"), bp::arg("j") )) );
        AtomsIn_SegID__exposer.def( bp::init< SireMol::AtomsIn< SireMol::SegID > const & >(( bp::arg("other") )) );
        { //::SireMol::AtomsIn< SireMol::SegID >::hash
        
            typedef SireMol::AtomsIn< SireMol::SegID > exported_class_t;
            typedef ::uint ( ::SireMol::AtomsIn< SireMol::SegID >::*hash_function_type )(  ) const;
            hash_function_type hash_function_value( &::SireMol::AtomsIn< SireMol::SegID >::hash );
            
            AtomsIn_SegID__exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireMol::AtomsIn< SireMol::SegID >::isNull
        
            typedef SireMol::AtomsIn< SireMol::SegID > exported_class_t;
            typedef bool ( ::SireMol::AtomsIn< SireMol::SegID >::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::AtomsIn< SireMol::SegID >::isNull );
            
            AtomsIn_SegID__exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireMol::AtomsIn< SireMol::SegID >::map
        
            typedef SireMol::AtomsIn< SireMol::SegID > exported_class_t;
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::AtomsIn< SireMol::SegID >::*map_function_type )( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::AtomsIn< SireMol::SegID >::map );
            
            AtomsIn_SegID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        AtomsIn_SegID__exposer.def( bp::self != bp::self );
        AtomsIn_SegID__exposer.def( bp::self != bp::other< SireID::ID >() );
        { //::SireMol::AtomsIn< SireMol::SegID >::operator=
        
            typedef SireMol::AtomsIn< SireMol::SegID > exported_class_t;
            typedef ::SireMol::AtomsIn< SireMol::SegID > & ( ::SireMol::AtomsIn< SireMol::SegID >::*assign_function_type )( ::SireMol::AtomsIn< SireMol::SegID > const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomsIn< SireMol::SegID >::operator= );
            
            AtomsIn_SegID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        AtomsIn_SegID__exposer.def( bp::self == bp::self );
        AtomsIn_SegID__exposer.def( bp::self == bp::other< SireID::ID >() );
        { //::SireMol::AtomsIn< SireMol::SegID >::toString
        
            typedef SireMol::AtomsIn< SireMol::SegID > exported_class_t;
            typedef ::QString ( ::SireMol::AtomsIn< SireMol::SegID >::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomsIn< SireMol::SegID >::toString );
            
            AtomsIn_SegID__exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::AtomsIn< SireMol::SegID >::typeName
        
            typedef SireMol::AtomsIn< SireMol::SegID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomsIn< SireMol::SegID >::typeName );
            
            AtomsIn_SegID__exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::AtomsIn< SireMol::SegID >::what
        
            typedef SireMol::AtomsIn< SireMol::SegID > exported_class_t;
            typedef char const * ( ::SireMol::AtomsIn< SireMol::SegID >::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMol::AtomsIn< SireMol::SegID >::what );
            
            AtomsIn_SegID__exposer.def( 
                "what"
                , what_function_value );
        
        }
        AtomsIn_SegID__exposer.staticmethod( "typeName" );
        AtomsIn_SegID__exposer.def( "__copy__", &__copy__);
        AtomsIn_SegID__exposer.def( "__deepcopy__", &__copy__);
        AtomsIn_SegID__exposer.def( "clone", &__copy__);
        AtomsIn_SegID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomsIn<SireMol::SegID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomsIn_SegID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomsIn<SireMol::SegID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomsIn_SegID__exposer.def( "__str__", &__str__< ::SireMol::AtomsIn<SireMol::SegID> > );
        AtomsIn_SegID__exposer.def( "__repr__", &__str__< ::SireMol::AtomsIn<SireMol::SegID> > );
    }

}