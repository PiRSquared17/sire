// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MatchMRData.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

void register_MatchMRData_class(){

    bp::class_< SireDB::MatchMRData >( "MatchMRData" )    
        .def( bp::init< >() )    
        .def( bp::init< SireDB::MatchMol const & >(( bp::arg("moldata") )) )    
        .def( bp::init< SireDB::MatchRes const & >(( bp::arg("resdata") )) )    
        .def( bp::init< SireDB::MatchMol const &, SireDB::MatchRes const & >(( bp::arg("moldata"), bp::arg("resdata") )) )    
        .def( bp::init< SireDB::MatchMR const & >(( bp::arg("matchmr") )) )    
        .def( bp::init< uint, SireDB::MatchRes const & >(( bp::arg("atm"), bp::arg("resdata") )) )    
        .def( 
            "getMatch"
            , (::SireDB::MatchMR ( ::SireDB::MatchMRData::* )( ::uint ) const)( &::SireDB::MatchMRData::getMatch )
            , ( bp::arg("atm") ) )    
        .def( 
            "hasMoleculeCriteria"
            , (bool ( ::SireDB::MatchMRData::* )(  ) const)( &::SireDB::MatchMRData::hasMoleculeCriteria ) )    
        .def( 
            "hasResidueCriteria"
            , (bool ( ::SireDB::MatchMRData::* )(  ) const)( &::SireDB::MatchMRData::hasResidueCriteria ) )    
        .def( 
            "isEmpty"
            , (bool ( ::SireDB::MatchMRData::* )(  ) const)( &::SireDB::MatchMRData::isEmpty ) )    
        .def( 
            "nMatches"
            , (int ( ::SireDB::MatchMRData::* )(  ) const)( &::SireDB::MatchMRData::nMatches ) )    
        .def( 
            "set"
            , (void ( ::SireDB::MatchMRData::* )( ::uint,::SireDB::MatchMR const & ) )( &::SireDB::MatchMRData::set )
            , ( bp::arg("atm"), bp::arg("match") ) )    
        .def( 
            "simplify"
            , (void ( ::SireDB::MatchMRData::* )( ::uint ) )( &::SireDB::MatchMRData::simplify )
            , ( bp::arg("maxatm") ) )    
        .def( 
            "toString"
            , (::QString ( ::SireDB::MatchMRData::* )(  ) const)( &::SireDB::MatchMRData::toString ) );

}
