// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "SystemMonitors.pypp.hpp"
#include "boost/python.hpp"
#include "siresystem_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireFF/ffgroupid.h"

namespace bp = boost::python;

SireSystem::SystemMonitors __copy__(const SireSystem::SystemMonitors &other){ return SireSystem::SystemMonitors(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireSystem::SystemMonitors&){ return "SireSystem::SystemMonitors";}

void register_SystemMonitors_class(){

    bp::class_< SireSystem::SystemMonitors >( "SystemMonitors" )    
        .def( bp::init< >() )    
        .def( 
            "assertCompatibleWith"
            , (void ( ::SireSystem::SystemMonitors::* )( ::SireSystem::QuerySystem const & ) const)( &::SireSystem::SystemMonitors::assertCompatibleWith )
            , ( bp::arg("system") ) )    
        .def( 
            "assertContains"
            , (void ( ::SireSystem::SystemMonitors::* )( ::SireCAS::Symbol const & ) const)( &::SireSystem::SystemMonitors::assertContains )
            , ( bp::arg("symbol") ) )    
        .def( 
            "clear"
            , (void ( ::SireSystem::SystemMonitors::* )(  ) )( &::SireSystem::SystemMonitors::clear ) )    
        .def( 
            "contains"
            , (bool ( ::SireSystem::SystemMonitors::* )( ::SireCAS::Symbol const & ) const)( &::SireSystem::SystemMonitors::contains )
            , ( bp::arg("symbol") ) )    
        .def( 
            "count"
            , (bool ( ::SireSystem::SystemMonitors::* )(  ) const)( &::SireSystem::SystemMonitors::count ) )    
        .def( 
            "forceUpdate"
            , (void ( ::SireSystem::SystemMonitors::* )( ::SireSystem::QuerySystem & ) )( &::SireSystem::SystemMonitors::forceUpdate )
            , ( bp::arg("system") ) )    
        .def( 
            "getDelta"
            , (::quint64 ( ::SireSystem::SystemMonitors::* )( ::SireCAS::Symbol const & ) const)( &::SireSystem::SystemMonitors::getDelta )
            , ( bp::arg("symbol") ) )    
        .def( 
            "isEmpty"
            , (bool ( ::SireSystem::SystemMonitors::* )(  ) const)( &::SireSystem::SystemMonitors::isEmpty ) )    
        .def( 
            "monitor"
            , (::SireSystem::SystemMonitor const & ( ::SireSystem::SystemMonitors::* )( ::SireCAS::Symbol const & ) const)( &::SireSystem::SystemMonitors::monitor )
            , ( bp::arg("symbol") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "monitor"
            , (void ( ::SireSystem::SystemMonitors::* )( ::SireSystem::QuerySystem & ) )( &::SireSystem::SystemMonitors::monitor )
            , ( bp::arg("system") ) )    
        .def( 
            "monitors"
            , (::QHash<SireCAS::Symbol, SireSystem::SystemMonitor> const & ( ::SireSystem::SystemMonitors::* )(  ) const)( &::SireSystem::SystemMonitors::monitors )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "nMonitors"
            , (bool ( ::SireSystem::SystemMonitors::* )(  ) const)( &::SireSystem::SystemMonitors::nMonitors ) )    
        .def( 
            "nUpdates"
            , (::quint64 ( ::SireSystem::SystemMonitors::* )(  ) const)( &::SireSystem::SystemMonitors::nUpdates ) )    
        .def( 
            "remove"
            , (void ( ::SireSystem::SystemMonitors::* )( ::SireCAS::Symbol const & ) )( &::SireSystem::SystemMonitors::remove )
            , ( bp::arg("symbol") ) )    
        .def( 
            "resetStatistics"
            , (void ( ::SireSystem::SystemMonitors::* )(  ) )( &::SireSystem::SystemMonitors::resetStatistics ) )    
        .def( 
            "set"
            , (void ( ::SireSystem::SystemMonitors::* )( ::SireCAS::Symbol const &,::SireSystem::SystemMonitor const &,::quint64 ) )( &::SireSystem::SystemMonitors::set )
            , ( bp::arg("symbol"), bp::arg("monitor"), bp::arg("delta")=(long long unsigned int)(1) ) )    
        .def( 
            "setDelta"
            , (void ( ::SireSystem::SystemMonitors::* )( ::SireCAS::Symbol const &,::quint64 ) )( &::SireSystem::SystemMonitors::setDelta )
            , ( bp::arg("symbol"), bp::arg("delta") ) )    
        .def( 
            "take"
            , (::SireSystem::SystemMonitor ( ::SireSystem::SystemMonitors::* )( ::SireCAS::Symbol const & ) )( &::SireSystem::SystemMonitors::take )
            , ( bp::arg("symbol") ) )    
        .def( 
            "update"
            , (void ( ::SireSystem::SystemMonitors::* )( ::SireSystem::QuerySystem & ) )( &::SireSystem::SystemMonitors::update )
            , ( bp::arg("system") ) )    
        .def( 
            "value"
            , (double ( ::SireSystem::SystemMonitors::* )( ::SireCAS::Symbol const & ) const)( &::SireSystem::SystemMonitors::value )
            , ( bp::arg("symbol") ) )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireSystem::SystemMonitors >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireSystem::SystemMonitors >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
