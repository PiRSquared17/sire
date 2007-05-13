// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MovesBase.pypp.hpp"
#include "boost/python.hpp"
#include "siresystem_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireFF/ffgroupid.h"

namespace bp = boost::python;

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireSystem::MovesBase&){ return "SireSystem::MovesBase";}

void register_MovesBase_class(){

    bp::class_< SireSystem::MovesBase, boost::noncopyable >( "MovesBase", bp::no_init )    
        .def( 
            "count"
            , (int ( ::SireSystem::MovesBase::* )(  ) const)( &::SireSystem::MovesBase::count ) )    
        .def( 
            "initialise"
            , (void ( ::SireSystem::MovesBase::* )( ::SireSystem::QuerySystem & ) )( &::SireSystem::MovesBase::initialise )
            , ( bp::arg("system") ) )    
        .def( 
            "moves"
            , (::QList<SireSystem::Move> ( ::SireSystem::MovesBase::* )(  ) const)( &::SireSystem::MovesBase::moves ) )    
        .def( 
            "setEnergyComponent"
            , (void ( ::SireSystem::MovesBase::* )( ::SireCAS::Symbol const & ) )( &::SireSystem::MovesBase::setEnergyComponent )
            , ( bp::arg("symbol") ) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireSystem::MovesBase::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireSystem::MovesBase::* )(  ) const)( &::SireSystem::MovesBase::what ) )    
        .staticmethod( "typeName" )    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireSystem::MovesBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireSystem::MovesBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
