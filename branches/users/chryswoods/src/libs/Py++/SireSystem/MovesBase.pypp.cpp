// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MovesBase.pypp.hpp"
#include "boost/python.hpp"
#include "siresystem_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"

namespace bp = boost::python;

void register_MovesBase_class(){

    bp::class_< SireSystem::MovesBase, boost::noncopyable >( "MovesBase", bp::no_init )    
        .def( 
            "count"
            , &::SireSystem::MovesBase::count )    
        .def( 
            "initialise"
            , &::SireSystem::MovesBase::initialise
            , ( bp::arg("system") ) )    
        .def( 
            "typeName"
            , &::SireSystem::MovesBase::typeName )    
        .def( 
            "what"
            , &::SireSystem::MovesBase::what )    
        .staticmethod( "typeName" );

}
