// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "TableBase.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

void register_TableBase_class(){

    bp::class_< SireDB::TableBase, boost::noncopyable >( "TableBase", bp::no_init )    
        .def( 
            "add"
            , (void ( ::SireDB::TableBase::* )( ::SireDB::TableBase const & ) )( &::SireDB::TableBase::add )
            , ( bp::arg("other") ) )    
        .def( 
            "assertCompatibleWith"
            , (void ( ::SireDB::TableBase::* )( ::SireMol::Molecule const & ) const)( &::SireDB::TableBase::assertCompatibleWith )
            , ( bp::arg("molecule") ) )    
        .def( 
            "info"
            , (::SireMol::MoleculeInfo const & ( ::SireDB::TableBase::* )(  ) const)( &::SireDB::TableBase::info )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "isCompatibleWith"
            , (bool ( ::SireDB::TableBase::* )( ::SireMol::MoleculeInfo const & ) const)( &::SireDB::TableBase::isCompatibleWith )
            , ( bp::arg("molinfo") ) )    
        .def( 
            "isEmpty"
            , (bool ( ::SireDB::TableBase::* )(  ) const)( &::SireDB::TableBase::isEmpty ) )    
        .def( 
            "what"
            , (char const * ( ::SireDB::TableBase::* )(  ) const)( &::SireDB::TableBase::what ) );

}
