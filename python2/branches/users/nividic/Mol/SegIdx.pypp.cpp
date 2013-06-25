// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SegIdx.pypp.hpp"

namespace bp = boost::python;

#include "segidx.h"

#include "segidx.h"

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

SireMol::SegIdx __copy__(const SireMol::SegIdx &other){ return SireMol::SegIdx(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SegIdx_class(){

    { //::SireMol::SegIdx
        typedef bp::class_< SireMol::SegIdx, bp::bases< SireMol::SegID, SireID::ID, SireID::IndexBase > > SegIdx_exposer_t;
        SegIdx_exposer_t SegIdx_exposer = SegIdx_exposer_t( "SegIdx" );
        bp::scope SegIdx_scope( SegIdx_exposer );
        SegIdx_exposer.def( bp::init< >() );
        SegIdx_exposer.def( bp::init< quint32 >(( bp::arg("idx") )) );
        SegIdx_exposer.def( bp::init< SireMol::SegIdx const & >(( bp::arg("other") )) );
        { //::SireMol::SegIdx::hash
        
            typedef ::uint ( ::SireMol::SegIdx::*hash_function_type )(  ) const;
            hash_function_type hash_function_value( &::SireMol::SegIdx::hash );
            
            SegIdx_exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireMol::SegIdx::isNull
        
            typedef bool ( ::SireMol::SegIdx::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::SegIdx::isNull );
            
            SegIdx_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireMol::SegIdx::map
        
            typedef ::QList< SireMol::SegIdx > ( ::SireMol::SegIdx::*map_function_type )( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::SegIdx::map );
            
            SegIdx_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::SegIdx::null
        
            typedef ::SireMol::SegIdx ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMol::SegIdx::null );
            
            SegIdx_exposer.def( 
                "null"
                , null_function_value );
        
        }
        { //::SireMol::SegIdx::operator=
        
            typedef ::SireMol::SegIdx & ( ::SireMol::SegIdx::*assign_function_type )( ::SireMol::SegIdx const & ) ;
            assign_function_type assign_function_value( &::SireMol::SegIdx::operator= );
            
            SegIdx_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::SegIdx::toString
        
            typedef ::QString ( ::SireMol::SegIdx::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::SegIdx::toString );
            
            SegIdx_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::SegIdx::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::SegIdx::typeName );
            
            SegIdx_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::SegIdx::what
        
            typedef char const * ( ::SireMol::SegIdx::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMol::SegIdx::what );
            
            SegIdx_exposer.def( 
                "what"
                , what_function_value );
        
        }
        SegIdx_exposer.staticmethod( "null" );
        SegIdx_exposer.staticmethod( "typeName" );
        SegIdx_exposer.def( "__copy__", &__copy__);
        SegIdx_exposer.def( "__deepcopy__", &__copy__);
        SegIdx_exposer.def( "clone", &__copy__);
        SegIdx_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::SegIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SegIdx_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::SegIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SegIdx_exposer.def( "__str__", &__str__< ::SireMol::SegIdx > );
        SegIdx_exposer.def( "__repr__", &__str__< ::SireMol::SegIdx > );
    }

}
