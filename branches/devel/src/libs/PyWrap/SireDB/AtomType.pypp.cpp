// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "AtomType.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

void register_AtomType_class(){

    bp::class_< SireDB::AtomType >( "AtomType" )    
        .def( bp::init< >() )    
        .def( bp::init< QString const & >(( bp::arg("id") )) )    
        .def( bp::init< QString const &, SireMol::Element const & >(( bp::arg("id"), bp::arg("element") )) )    
        .def( bp::init< SireMol::Element const & >(( bp::arg("element") )) )    
        .def( bp::init< SireDB::AtomParameter<SireDB::AtomType> const & >(( bp::arg("atomtype") )) )    
        .def( 
            "ID"
            , (::QString ( ::SireDB::AtomType::* )(  ) const)( &::SireDB::AtomType::ID ) )    
        .def( 
            "dummy"
            , (::SireDB::AtomType (*)(  ))( &::SireDB::AtomType::dummy ) )    
        .def( 
            "element"
            , (::SireMol::Element const & ( ::SireDB::AtomType::* )(  ) const)( &::SireDB::AtomType::element )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "isDummy"
            , (bool ( ::SireDB::AtomType::* )(  ) const)( &::SireDB::AtomType::isDummy ) )    
        .def( 
            "isWild"
            , (bool ( ::SireDB::AtomType::* )(  ) const)( &::SireDB::AtomType::isWild ) )    
        .def( 
            "isWild"
            , (bool ( ::SireDB::AtomType::* )( ::SireMol::Element const & ) const)( &::SireDB::AtomType::isWild )
            , ( bp::arg("element") ) )    
        .def( 
            "isWildElement"
            , (bool ( ::SireDB::AtomType::* )(  ) const)( &::SireDB::AtomType::isWildElement ) )    
        .def( 
            "isWildType"
            , (bool ( ::SireDB::AtomType::* )(  ) const)( &::SireDB::AtomType::isWildType ) )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "toString"
            , (::QString ( ::SireDB::AtomType::* )(  ) const)( &::SireDB::AtomType::toString ) )    
        .def( 
            "wild"
            , (::SireDB::AtomType (*)(  ))( &::SireDB::AtomType::wild ) )    
        .def( 
            "wild"
            , (::SireDB::AtomType (*)( ::SireMol::Element const & ))( &::SireDB::AtomType::wild )
            , ( bp::arg("element") ) )    
        .staticmethod( "dummy" )    
        .staticmethod( "wild" );

}
