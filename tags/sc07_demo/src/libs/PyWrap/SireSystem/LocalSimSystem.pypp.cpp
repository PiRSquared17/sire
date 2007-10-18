// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "LocalSimSystem.pypp.hpp"
#include "boost/python.hpp"
#include "siresystem_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireFF/ffgroupid.h"

namespace bp = boost::python;

const char* pvt_get_name(const SireSystem::LocalSimSystem&){ return "SireSystem::LocalSimSystem";}

void register_LocalSimSystem_class(){

    bp::class_< SireSystem::LocalSimSystem, bp::bases< SireSystem::SimSystem >, boost::noncopyable >( "LocalSimSystem", bp::no_init )    
        .def( 
            "rollBack"
            , (void ( ::SireSystem::LocalSimSystem::* )( ::SireSystem::CheckPoint const & ) )( &::SireSystem::LocalSimSystem::rollBack )
            , ( bp::arg("checkpoint") ) )    
        .def( "__str__", &pvt_get_name);

}
