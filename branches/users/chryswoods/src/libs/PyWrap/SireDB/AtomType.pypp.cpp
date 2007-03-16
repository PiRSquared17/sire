// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "AtomType.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

SireDB::AtomType __copy__(const SireDB::AtomType &other){ return SireDB::AtomType(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

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
        .staticmethod( "wild" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireDB::AtomType >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireDB::AtomType >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireDB::AtomType > );

}
