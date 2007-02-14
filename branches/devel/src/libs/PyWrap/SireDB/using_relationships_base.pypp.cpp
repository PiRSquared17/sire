// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "using_relationships_base.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

void register_using_relationships_base_class(){

    bp::class_< SireDB::using_relationships_base, bp::bases< SireDB::using_database > >( "using_relationships_base" )    
        .def( bp::init< >() )    
        .def( bp::init< QStringList const & >(( bp::arg("typenames") )) )    
        .def( 
            "what"
            , (char const * ( ::SireDB::using_relationships_base::* )(  ) const)( &::SireDB::using_relationships_base::what ) );

}
