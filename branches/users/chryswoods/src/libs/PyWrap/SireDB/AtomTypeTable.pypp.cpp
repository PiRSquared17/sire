// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "AtomTypeTable.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

void register_AtomTypeTable_class(){

    bp::class_< SireDB::AtomTypeTable, bp::bases< SireDB::AtomTableT<SireDB::AtomType> > >( "AtomTypeTable" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") )) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireDB::AtomTypeTable::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireDB::AtomTypeTable::* )(  ) const)( &::SireDB::AtomTypeTable::what ) )    
        .staticmethod( "typeName" );

}
