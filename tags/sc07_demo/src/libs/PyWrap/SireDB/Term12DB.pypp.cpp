// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Term12DB.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

#include "SirePy/str.hpp"

void register_Term12DB_class(){

    bp::class_< SireDB::Term12DB, bp::bases< SireDB::FuncDB >, boost::noncopyable >( "Term12DB", bp::no_init )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireDB::Term12DB::typeName ) )    
        .staticmethod( "typeName" )    
        .def( "__str__", &SirePy::__str__< ::SireDB::Term12DB > );

}
