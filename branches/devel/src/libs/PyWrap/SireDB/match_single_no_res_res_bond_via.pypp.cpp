// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "match_single_no_res_res_bond_via.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

void register_match_single_no_res_res_bond_via_class(){

    bp::class_< SireDB::match_single_no_res_res_bond_via >( "match_single_no_res_res_bond_via" )    
        .def( bp::init< >() )    
        .def( bp::self == bp::other< QString >() );

}
