// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "match_res_num.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

void register_match_res_num_class(){

    bp::class_< SireDB::match_res_num >( "match_res_num" )    
        .def( bp::init< >() )    
        .def( bp::init< uint >(( bp::arg("atmid") )) )    
        .def( bp::self == bp::other< SireMol::ResNum >() );

}
