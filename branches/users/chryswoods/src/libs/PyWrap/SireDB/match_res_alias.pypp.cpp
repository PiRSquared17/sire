// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "match_res_alias.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

SireDB::match_res_alias __copy__(const SireDB::match_res_alias &other){ return SireDB::match_res_alias(other); }

const char* pvt_get_name(const SireDB::match_res_alias&){ return "SireDB::match_res_alias";}

void register_match_res_alias_class(){

    bp::class_< SireDB::match_res_alias >( "match_res_alias" )    
        .def( bp::init< >() )    
        .def( bp::init< uint >(( bp::arg("atmid") )) )    
        .def( bp::self == bp::other< QString >() )    
        .def( "__copy__", &__copy__)    
        .def( "__str__", &pvt_get_name);

}
