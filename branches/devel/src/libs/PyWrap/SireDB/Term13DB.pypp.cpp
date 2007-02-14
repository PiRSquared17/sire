// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Term13DB.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

void register_Term13DB_class(){

    bp::class_< SireDB::Term13DB, bp::bases< SireDB::FuncDB >, boost::noncopyable >( "Term13DB", bp::no_init )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireDB::Term13DB::typeName ) )    
        .staticmethod( "typeName" );

}
