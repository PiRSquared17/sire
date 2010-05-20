// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomsIn_ResID_.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "atom.h"

#include "chain.h"

#include "chainresid.h"

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

#include "resid.h"

#include "residentifier.h"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include "resid.h"

SireMol::AtomsIn<SireMol::ResID> __copy__(const SireMol::AtomsIn<SireMol::ResID> &other){ return SireMol::AtomsIn<SireMol::ResID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomsIn_ResID__class(){

    { //::SireMol::AtomsIn< SireMol::ResID >
        typedef bp::class_< SireMol::AtomsIn< SireMol::ResID >, bp::bases< SireMol::AtomID, SireID::ID > > AtomsIn_ResID__exposer_t;
        AtomsIn_ResID__exposer_t AtomsIn_ResID__exposer = AtomsIn_ResID__exposer_t( "AtomsIn_ResID_" );
        bp::scope AtomsIn_ResID__scope( AtomsIn_ResID__exposer );
        AtomsIn_ResID__exposer.def( bp::init< >() );
        AtomsIn_ResID__exposer.def( bp::init< SireMol::ResID const & >(( bp::arg("id") )) );
        AtomsIn_ResID__exposer.def( bp::init< SireMol::ResID const &, qint32 >(( bp::arg("id"), bp::arg("i") )) );
        AtomsIn_ResID__exposer.def( bp::init< SireMol::ResID const &, qint32, qint32 >(( bp::arg("id"), bp::arg("i"), bp::arg("j") )) );
        AtomsIn_ResID__exposer.def( bp::init< SireMol::AtomsIn< SireMol::ResID > const & >(( bp::arg("other") )) );
        { //::SireMol::AtomsIn< SireMol::ResID >::hash
        
            typedef SireMol::AtomsIn< SireMol::ResID > exported_class_t;
            typedef ::uint ( ::SireMol::AtomsIn< SireMol::ResID >::*hash_function_type )(  ) const;
            hash_function_type hash_function_value( &::SireMol::AtomsIn< SireMol::ResID >::hash );
            
            AtomsIn_ResID__exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireMol::AtomsIn< SireMol::ResID >::isNull
        
            typedef SireMol::AtomsIn< SireMol::ResID > exported_class_t;
            typedef bool ( ::SireMol::AtomsIn< SireMol::ResID >::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::AtomsIn< SireMol::ResID >::isNull );
            
            AtomsIn_ResID__exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireMol::AtomsIn< SireMol::ResID >::map
        
            typedef SireMol::AtomsIn< SireMol::ResID > exported_class_t;
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::AtomsIn< SireMol::ResID >::*map_function_type )( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::AtomsIn< SireMol::ResID >::map );
            
            AtomsIn_ResID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        AtomsIn_ResID__exposer.def( bp::self != bp::self );
        AtomsIn_ResID__exposer.def( bp::self != bp::other< SireID::ID >() );
        { //::SireMol::AtomsIn< SireMol::ResID >::operator=
        
            typedef SireMol::AtomsIn< SireMol::ResID > exported_class_t;
            typedef ::SireMol::AtomsIn< SireMol::ResID > & ( ::SireMol::AtomsIn< SireMol::ResID >::*assign_function_type )( ::SireMol::AtomsIn< SireMol::ResID > const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomsIn< SireMol::ResID >::operator= );
            
            AtomsIn_ResID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        AtomsIn_ResID__exposer.def( bp::self == bp::self );
        AtomsIn_ResID__exposer.def( bp::self == bp::other< SireID::ID >() );
        { //::SireMol::AtomsIn< SireMol::ResID >::toString
        
            typedef SireMol::AtomsIn< SireMol::ResID > exported_class_t;
            typedef ::QString ( ::SireMol::AtomsIn< SireMol::ResID >::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomsIn< SireMol::ResID >::toString );
            
            AtomsIn_ResID__exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::AtomsIn< SireMol::ResID >::typeName
        
            typedef SireMol::AtomsIn< SireMol::ResID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomsIn< SireMol::ResID >::typeName );
            
            AtomsIn_ResID__exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::AtomsIn< SireMol::ResID >::what
        
            typedef SireMol::AtomsIn< SireMol::ResID > exported_class_t;
            typedef char const * ( ::SireMol::AtomsIn< SireMol::ResID >::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMol::AtomsIn< SireMol::ResID >::what );
            
            AtomsIn_ResID__exposer.def( 
                "what"
                , what_function_value );
        
        }
        AtomsIn_ResID__exposer.staticmethod( "typeName" );
        AtomsIn_ResID__exposer.def( "__copy__", &__copy__);
        AtomsIn_ResID__exposer.def( "__deepcopy__", &__copy__);
        AtomsIn_ResID__exposer.def( "clone", &__copy__);
        AtomsIn_ResID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomsIn<SireMol::ResID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomsIn_ResID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomsIn<SireMol::ResID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomsIn_ResID__exposer.def( "__str__", &__str__< ::SireMol::AtomsIn<SireMol::ResID> > );
        AtomsIn_ResID__exposer.def( "__repr__", &__str__< ::SireMol::AtomsIn<SireMol::ResID> > );
    }

}
