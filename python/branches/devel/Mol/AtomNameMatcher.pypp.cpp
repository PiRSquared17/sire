// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomNameMatcher.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireUnits/units.h"

#include "atomidentifier.h"

#include "atomidx.h"

#include "atommatcher.h"

#include "atommatchers.h"

#include "atomname.h"

#include "atomselection.h"

#include "evaluator.h"

#include "moleculeinfodata.h"

#include "moleculeview.h"

#include "tostring.h"

#include "atommatchers.h"

SireMol::AtomNameMatcher __copy__(const SireMol::AtomNameMatcher &other){ return SireMol::AtomNameMatcher(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomNameMatcher_class(){

    { //::SireMol::AtomNameMatcher
        typedef bp::class_< SireMol::AtomNameMatcher, bp::bases< SireMol::AtomMatcher, SireBase::Property > > AtomNameMatcher_exposer_t;
        AtomNameMatcher_exposer_t AtomNameMatcher_exposer = AtomNameMatcher_exposer_t( "AtomNameMatcher", bp::init< >() );
        bp::scope AtomNameMatcher_scope( AtomNameMatcher_exposer );
        AtomNameMatcher_exposer.def( bp::init< SireMol::AtomNameMatcher const & >(( bp::arg("arg0") )) );
        { //::SireMol::AtomNameMatcher::match
        
            typedef ::QHash< SireMol::AtomIdx, SireMol::AtomIdx > ( ::SireMol::AtomNameMatcher::*match_function_type )( ::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const & ) const;
            match_function_type match_function_value( &::SireMol::AtomNameMatcher::match );
            
            AtomNameMatcher_exposer.def( 
                "match"
                , match_function_value
                , ( bp::arg("molinfo0"), bp::arg("molinfo1") ) );
        
        }
        { //::SireMol::AtomNameMatcher::match
        
            typedef ::QHash< SireMol::AtomIdx, SireMol::AtomIdx > ( ::SireMol::AtomNameMatcher::*match_function_type )( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const &,::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            match_function_type match_function_value( &::SireMol::AtomNameMatcher::match );
            
            AtomNameMatcher_exposer.def( 
                "match"
                , match_function_value
                , ( bp::arg("molview0"), bp::arg("map0"), bp::arg("molview1"), bp::arg("map1") ) );
        
        }
        AtomNameMatcher_exposer.def( bp::self != bp::self );
        { //::SireMol::AtomNameMatcher::operator=
        
            typedef ::SireMol::AtomNameMatcher & ( ::SireMol::AtomNameMatcher::*assign_function_type )( ::SireMol::AtomNameMatcher const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomNameMatcher::operator= );
            
            AtomNameMatcher_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        AtomNameMatcher_exposer.def( bp::self == bp::self );
        { //::SireMol::AtomNameMatcher::toString
        
            typedef ::QString ( ::SireMol::AtomNameMatcher::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomNameMatcher::toString );
            
            AtomNameMatcher_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::AtomNameMatcher::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomNameMatcher::typeName );
            
            AtomNameMatcher_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::AtomNameMatcher::what
        
            typedef char const * ( ::SireMol::AtomNameMatcher::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMol::AtomNameMatcher::what );
            
            AtomNameMatcher_exposer.def( 
                "what"
                , what_function_value );
        
        }
        AtomNameMatcher_exposer.staticmethod( "typeName" );
        AtomNameMatcher_exposer.def( "__copy__", &__copy__);
        AtomNameMatcher_exposer.def( "__deepcopy__", &__copy__);
        AtomNameMatcher_exposer.def( "clone", &__copy__);
        AtomNameMatcher_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomNameMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomNameMatcher_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomNameMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomNameMatcher_exposer.def( "__str__", &__str__< ::SireMol::AtomNameMatcher > );
        AtomNameMatcher_exposer.def( "__repr__", &__str__< ::SireMol::AtomNameMatcher > );
    }

}
