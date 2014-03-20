// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomIdxMatcher.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "atomidentifier.h"

#include "atomidx.h"

#include "atommatcher.h"

#include "atomname.h"

#include "atomselection.h"

#include "evaluator.h"

#include "moleculeinfodata.h"

#include "moleculeview.h"

#include "atommatcher.h"

SireMol::AtomIdxMatcher __copy__(const SireMol::AtomIdxMatcher &other){ return SireMol::AtomIdxMatcher(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomIdxMatcher_class(){

    { //::SireMol::AtomIdxMatcher
        typedef bp::class_< SireMol::AtomIdxMatcher, bp::bases< SireMol::AtomMatcher, SireBase::Property > > AtomIdxMatcher_exposer_t;
        AtomIdxMatcher_exposer_t AtomIdxMatcher_exposer = AtomIdxMatcher_exposer_t( "AtomIdxMatcher", bp::init< >() );
        bp::scope AtomIdxMatcher_scope( AtomIdxMatcher_exposer );
        AtomIdxMatcher_exposer.def( bp::init< SireMol::AtomIdxMatcher const & >(( bp::arg("arg0") )) );
        { //::SireMol::AtomIdxMatcher::match
        
            typedef ::QHash< SireMol::AtomIdx, SireMol::AtomIdx > ( ::SireMol::AtomIdxMatcher::*match_function_type )( ::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const & ) const;
            match_function_type match_function_value( &::SireMol::AtomIdxMatcher::match );
            
            AtomIdxMatcher_exposer.def( 
                "match"
                , match_function_value
                , ( bp::arg("molinfo0"), bp::arg("molinfo1") ) );
        
        }
        { //::SireMol::AtomIdxMatcher::match
        
            typedef ::QHash< SireMol::AtomIdx, SireMol::AtomIdx > ( ::SireMol::AtomIdxMatcher::*match_function_type )( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const &,::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            match_function_type match_function_value( &::SireMol::AtomIdxMatcher::match );
            
            AtomIdxMatcher_exposer.def( 
                "match"
                , match_function_value
                , ( bp::arg("molview0"), bp::arg("map0"), bp::arg("molview1"), bp::arg("map1") ) );
        
        }
        AtomIdxMatcher_exposer.def( bp::self != bp::self );
        { //::SireMol::AtomIdxMatcher::operator=
        
            typedef ::SireMol::AtomIdxMatcher & ( ::SireMol::AtomIdxMatcher::*assign_function_type )( ::SireMol::AtomIdxMatcher const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomIdxMatcher::operator= );
            
            AtomIdxMatcher_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        AtomIdxMatcher_exposer.def( bp::self == bp::self );
        { //::SireMol::AtomIdxMatcher::toString
        
            typedef ::QString ( ::SireMol::AtomIdxMatcher::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomIdxMatcher::toString );
            
            AtomIdxMatcher_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::AtomIdxMatcher::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomIdxMatcher::typeName );
            
            AtomIdxMatcher_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::AtomIdxMatcher::what
        
            typedef char const * ( ::SireMol::AtomIdxMatcher::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMol::AtomIdxMatcher::what );
            
            AtomIdxMatcher_exposer.def( 
                "what"
                , what_function_value );
        
        }
        AtomIdxMatcher_exposer.staticmethod( "typeName" );
        AtomIdxMatcher_exposer.def( "__copy__", &__copy__);
        AtomIdxMatcher_exposer.def( "__deepcopy__", &__copy__);
        AtomIdxMatcher_exposer.def( "clone", &__copy__);
        AtomIdxMatcher_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomIdxMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomIdxMatcher_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomIdxMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomIdxMatcher_exposer.def( "__str__", &__str__< ::SireMol::AtomIdxMatcher > );
        AtomIdxMatcher_exposer.def( "__repr__", &__str__< ::SireMol::AtomIdxMatcher > );
    }

}
