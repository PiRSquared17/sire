// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomMatcher.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "atomidx.h"

#include "atommatcher.h"

#include "atomname.h"

#include "moleculeinfodata.h"

#include "moleculeview.h"

#include "atommatcher.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomMatcher_class(){

    { //::SireMol::AtomMatcher
        typedef bp::class_< SireMol::AtomMatcher, bp::bases< SireBase::Property >, boost::noncopyable > AtomMatcher_exposer_t;
        AtomMatcher_exposer_t AtomMatcher_exposer = AtomMatcher_exposer_t( "AtomMatcher", bp::no_init );
        bp::scope AtomMatcher_scope( AtomMatcher_exposer );
        { //::SireMol::AtomMatcher::match
        
            typedef ::QHash< SireMol::AtomIdx, SireMol::AtomIdx > ( ::SireMol::AtomMatcher::*match_function_type )( ::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const & ) const;
            match_function_type match_function_value( &::SireMol::AtomMatcher::match );
            
            AtomMatcher_exposer.def( 
                "match"
                , match_function_value
                , ( bp::arg("molinfo0"), bp::arg("molinfo1") ) );
        
        }
        { //::SireMol::AtomMatcher::match
        
            typedef ::QHash< SireMol::AtomIdx, SireMol::AtomIdx > ( ::SireMol::AtomMatcher::*match_function_type )( ::SireMol::MoleculeView const &,::SireMol::MoleculeView const & ) const;
            match_function_type match_function_value( &::SireMol::AtomMatcher::match );
            
            AtomMatcher_exposer.def( 
                "match"
                , match_function_value
                , ( bp::arg("molview0"), bp::arg("molview1") ) );
        
        }
        { //::SireMol::AtomMatcher::toString
        
            typedef ::QString ( ::SireMol::AtomMatcher::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomMatcher::toString );
            
            AtomMatcher_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::AtomMatcher::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomMatcher::typeName );
            
            AtomMatcher_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::AtomMatcher::unchangedAtomOrder
        
            typedef bool ( ::SireMol::AtomMatcher::*unchangedAtomOrder_function_type )( ::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const & ) const;
            unchangedAtomOrder_function_type unchangedAtomOrder_function_value( &::SireMol::AtomMatcher::unchangedAtomOrder );
            
            AtomMatcher_exposer.def( 
                "unchangedAtomOrder"
                , unchangedAtomOrder_function_value
                , ( bp::arg("molinfo0"), bp::arg("molinfo1") ) );
        
        }
        AtomMatcher_exposer.staticmethod( "typeName" );
        AtomMatcher_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomMatcher_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomMatcher_exposer.def( "__str__", &__str__< ::SireMol::AtomMatcher > );
        AtomMatcher_exposer.def( "__repr__", &__str__< ::SireMol::AtomMatcher > );
    }

}