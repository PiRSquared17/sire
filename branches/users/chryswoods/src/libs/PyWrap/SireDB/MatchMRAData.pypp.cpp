// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MatchMRAData.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

void register_MatchMRAData_class(){

    bp::class_< SireDB::MatchMRAData, bp::bases< SireDB::MatchData > >( "MatchMRAData" )    
        .def( bp::init< >() )    
        .def( bp::init< SireDB::MatchMol const & >(( bp::arg("moldata") )) )    
        .def( bp::init< SireDB::MatchRes const & >(( bp::arg("resdata") )) )    
        .def( bp::init< SireDB::MatchAtom const & >(( bp::arg("atmdata") )) )    
        .def( bp::init< SireDB::MatchMR const & >(( bp::arg("matchmr") )) )    
        .def( bp::init< SireDB::MatchMRData const & >(( bp::arg("molresdata") )) )    
        .def( bp::init< SireDB::MatchMRA const & >(( bp::arg("match") )) )    
        .def( bp::init< uint, SireDB::MatchRes const & >(( bp::arg("atm"), bp::arg("resdata") )) )    
        .def( bp::init< uint, SireDB::MatchAtom const & >(( bp::arg("atm"), bp::arg("atmdata") )) )    
        .def( 
            "addMirroredTo"
            , (::SireDB::RelateID ( ::SireDB::MatchMRAData::* )( ::SireDB::ParameterDB &,::uint ) const)( &::SireDB::MatchMRAData::addMirroredTo )
            , ( bp::arg("db"), bp::arg("n") ) )    
        .def( 
            "addTo"
            , (::SireDB::RelateID ( ::SireDB::MatchMRAData::* )( ::SireDB::ParameterDB &,::uint ) const)( &::SireDB::MatchMRAData::addTo )
            , ( bp::arg("db"), bp::arg("n") ) )    
        .def( 
            "getMatch"
            , (::SireDB::MatchMRA ( ::SireDB::MatchMRAData::* )( ::uint ) const)( &::SireDB::MatchMRAData::getMatch )
            , ( bp::arg("atm") ) )    
        .def( 
            "hasAtomCriteria"
            , (bool ( ::SireDB::MatchMRAData::* )(  ) const)( &::SireDB::MatchMRAData::hasAtomCriteria ) )    
        .def( 
            "hasMoleculeCriteria"
            , (bool ( ::SireDB::MatchMRAData::* )(  ) const)( &::SireDB::MatchMRAData::hasMoleculeCriteria ) )    
        .def( 
            "hasResidueCriteria"
            , (bool ( ::SireDB::MatchMRAData::* )(  ) const)( &::SireDB::MatchMRAData::hasResidueCriteria ) )    
        .def( 
            "isEmpty"
            , (bool ( ::SireDB::MatchMRAData::* )(  ) const)( &::SireDB::MatchMRAData::isEmpty ) )    
        .def( 
            "nMatches"
            , (int ( ::SireDB::MatchMRAData::* )(  ) const)( &::SireDB::MatchMRAData::nMatches ) )    
        .def( 
            "set"
            , (void ( ::SireDB::MatchMRAData::* )( ::uint,::SireDB::MatchMRA const & ) )( &::SireDB::MatchMRAData::set )
            , ( bp::arg("atm"), bp::arg("match") ) )    
        .def( 
            "simplify"
            , (void ( ::SireDB::MatchMRAData::* )( ::uint ) )( &::SireDB::MatchMRAData::simplify )
            , ( bp::arg("maxatm") ) )    
        .def( 
            "toString"
            , (::QString ( ::SireDB::MatchMRAData::* )(  ) const)( &::SireDB::MatchMRAData::toString ) );

}
