// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MatchMol.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

SireDB::MatchMol __copy__(const SireDB::MatchMol &other){ return SireDB::MatchMol(other); }

#include "SirePy/str.hpp"

void register_MatchMol_class(){

    bp::class_< SireDB::MatchMol >( "MatchMol" )    
        .def( bp::init< >() )    
        .def( 
            "group"
            , (::QString ( ::SireDB::MatchMol::* )(  ) const)( &::SireDB::MatchMol::group ) )    
        .def( 
            "hasGroupMatch"
            , (bool ( ::SireDB::MatchMol::* )(  ) const)( &::SireDB::MatchMol::hasGroupMatch ) )    
        .def( 
            "hasNameMatch"
            , (bool ( ::SireDB::MatchMol::* )(  ) const)( &::SireDB::MatchMol::hasNameMatch ) )    
        .def( 
            "hasStateMatch"
            , (bool ( ::SireDB::MatchMol::* )(  ) const)( &::SireDB::MatchMol::hasStateMatch ) )    
        .def( 
            "isEmpty"
            , (bool ( ::SireDB::MatchMol::* )(  ) const)( &::SireDB::MatchMol::isEmpty ) )    
        .def( 
            "name"
            , (::QString ( ::SireDB::MatchMol::* )(  ) const)( &::SireDB::MatchMol::name ) )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "setGroup"
            , (void ( ::SireDB::MatchMol::* )( ::QString const & ) )( &::SireDB::MatchMol::setGroup )
            , ( bp::arg("group") ) )    
        .def( 
            "setName"
            , (void ( ::SireDB::MatchMol::* )( ::QString const & ) )( &::SireDB::MatchMol::setName )
            , ( bp::arg("name") ) )    
        .def( 
            "setState"
            , (void ( ::SireDB::MatchMol::* )( ::QString const & ) )( &::SireDB::MatchMol::setState )
            , ( bp::arg("state") ) )    
        .def( 
            "state"
            , (::QString ( ::SireDB::MatchMol::* )(  ) const)( &::SireDB::MatchMol::state ) )    
        .def( 
            "toString"
            , (::QString ( ::SireDB::MatchMol::* )(  ) const)( &::SireDB::MatchMol::toString ) )    
        .def( "__copy__", &__copy__)    
        .def( "__str__", &SirePy::__str__< ::SireDB::MatchMol > );

}
