// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "match_single_atom_element.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

void register_match_single_atom_element_class(){

    bp::class_< SireDB::match_single_atom_element >( "match_single_atom_element" )    
        .def( bp::init< >() )    
        .def( bp::self == bp::other< QString >() )    
        .def( bp::self == bp::other< SireMol::Element >() )    
        .def( bp::self == bp::other< uint >() );

}
