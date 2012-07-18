// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ResNum.pypp.hpp"

namespace bp = boost::python;

#include "resnum.h"

#include "resnum.h"

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

#include "withres.h"

#include "resid.h"

SireMol::ResNum __copy__(const SireMol::ResNum &other){ return SireMol::ResNum(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ResNum_class(){

    { //::SireMol::ResNum
        typedef bp::class_< SireMol::ResNum, bp::bases< SireMol::ResID, SireID::ID, SireID::Number > > ResNum_exposer_t;
        ResNum_exposer_t ResNum_exposer = ResNum_exposer_t( "ResNum" );
        bp::scope ResNum_scope( ResNum_exposer );
        ResNum_exposer.def( bp::init< >() );
        ResNum_exposer.def( bp::init< quint32 >(( bp::arg("num") )) );
        ResNum_exposer.def( bp::init< SireMol::ResNum const & >(( bp::arg("other") )) );
        { //::SireMol::ResNum::hash
        
            typedef ::uint ( ::SireMol::ResNum::*hash_function_type )(  ) const;
            hash_function_type hash_function_value( &::SireMol::ResNum::hash );
            
            ResNum_exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireMol::ResNum::isNull
        
            typedef bool ( ::SireMol::ResNum::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::ResNum::isNull );
            
            ResNum_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireMol::ResNum::map
        
            typedef ::QList< SireMol::ResIdx > ( ::SireMol::ResNum::*map_function_type )( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::ResNum::map );
            
            ResNum_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        ResNum_exposer.def( bp::self != bp::self );
        ResNum_exposer.def( bp::self < bp::self );
        ResNum_exposer.def( bp::self <= bp::self );
        { //::SireMol::ResNum::operator=
        
            typedef ::SireMol::ResNum & ( ::SireMol::ResNum::*assign_function_type )( ::SireMol::ResNum const & ) ;
            assign_function_type assign_function_value( &::SireMol::ResNum::operator= );
            
            ResNum_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        ResNum_exposer.def( bp::self == bp::other< SireID::ID >() );
        ResNum_exposer.def( bp::self == bp::self );
        ResNum_exposer.def( bp::self > bp::self );
        ResNum_exposer.def( bp::self >= bp::self );
        { //::SireMol::ResNum::toString
        
            typedef ::QString ( ::SireMol::ResNum::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::ResNum::toString );
            
            ResNum_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::ResNum::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ResNum::typeName );
            
            ResNum_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::ResNum::what
        
            typedef char const * ( ::SireMol::ResNum::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMol::ResNum::what );
            
            ResNum_exposer.def( 
                "what"
                , what_function_value );
        
        }
        ResNum_exposer.staticmethod( "typeName" );
        ResNum_exposer.def( "__copy__", &__copy__);
        ResNum_exposer.def( "__deepcopy__", &__copy__);
        ResNum_exposer.def( "clone", &__copy__);
        ResNum_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ResNum >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResNum_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ResNum >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResNum_exposer.def( "__str__", &__str__< ::SireMol::ResNum > );
        ResNum_exposer.def( "__repr__", &__str__< ::SireMol::ResNum > );
    }

}
