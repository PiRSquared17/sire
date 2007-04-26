// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "SystemData.pypp.hpp"
#include "boost/python.hpp"
#include "siresystem_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireFF/ffgroupid.h"

namespace bp = boost::python;

SireSystem::SystemData __copy__(const SireSystem::SystemData &other){ return SireSystem::SystemData(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireSystem::SystemData&){ return "SireSystem::SystemData";}

void register_SystemData_class(){

    bp::class_< SireSystem::SystemData >( "SystemData" )    
        .def( bp::init< >() )    
        .def( bp::init< QString const & >(( bp::arg("name") )) )    
        .def( 
            "ID"
            , (::SireSystem::SystemID ( ::SireSystem::SystemData::* )(  ) const)( &::SireSystem::SystemData::ID ) )    
        .def( 
            "add"
            , (void ( ::SireSystem::SystemData::* )( ::SireMol::MoleculeGroup const & ) )( &::SireSystem::SystemData::add )
            , ( bp::arg("group") ) )    
        .def( 
            "add"
            , (void ( ::SireSystem::SystemData::* )( ::SireMol::MoleculeGroups const & ) )( &::SireSystem::SystemData::add )
            , ( bp::arg("groups") ) )    
        .def( 
            "add"
            , (::SireMol::PartialMolecule ( ::SireSystem::SystemData::* )( ::SireMol::PartialMolecule const &,::QSet<SireMol::MoleculeGroupID> const & ) )( &::SireSystem::SystemData::add )
            , ( bp::arg("molecule"), bp::arg("groupids") ) )    
        .def( 
            "add"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireSystem::SystemData::* )( ::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> const &,::QSet<SireMol::MoleculeGroupID> const & ) )( &::SireSystem::SystemData::add )
            , ( bp::arg("molecules"), bp::arg("groupids") ) )    
        .def( 
            "add"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireSystem::SystemData::* )( ::QList<SireMol::PartialMolecule> const &,::QSet<SireMol::MoleculeGroupID> const & ) )( &::SireSystem::SystemData::add )
            , ( bp::arg("molecules"), bp::arg("groupids") ) )    
        .def( 
            "change"
            , (::SireMol::PartialMolecule ( ::SireSystem::SystemData::* )( ::SireMol::PartialMolecule const & ) )( &::SireSystem::SystemData::change )
            , ( bp::arg("molecule") ) )    
        .def( 
            "change"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireSystem::SystemData::* )( ::QList<SireMol::PartialMolecule> const & ) )( &::SireSystem::SystemData::change )
            , ( bp::arg("molecules") ) )    
        .def( 
            "change"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireSystem::SystemData::* )( ::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> const & ) )( &::SireSystem::SystemData::change )
            , ( bp::arg("molecules") ) )    
        .def( 
            "group"
            , (::SireMol::MoleculeGroup const & ( ::SireSystem::SystemData::* )( ::SireMol::MoleculeGroupID ) const)( &::SireSystem::SystemData::group )
            , ( bp::arg("id") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "groups"
            , (::SireMol::MoleculeGroups const & ( ::SireSystem::SystemData::* )(  ) const)( &::SireSystem::SystemData::groups )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "incrementMajorVersion"
            , (void ( ::SireSystem::SystemData::* )(  ) )( &::SireSystem::SystemData::incrementMajorVersion ) )    
        .def( 
            "incrementMinorVersion"
            , (void ( ::SireSystem::SystemData::* )(  ) )( &::SireSystem::SystemData::incrementMinorVersion ) )    
        .def( 
            "mapIntoSystemSpace"
            , (::SireMol::PartialMolecule ( ::SireSystem::SystemData::* )( ::SireMol::PartialMolecule const & ) const)( &::SireSystem::SystemData::mapIntoSystemSpace )
            , ( bp::arg("molecule") ) )    
        .def( 
            "mapIntoSystemSpace"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireSystem::SystemData::* )( ::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> const & ) const)( &::SireSystem::SystemData::mapIntoSystemSpace )
            , ( bp::arg("molecules") ) )    
        .def( 
            "mapIntoSystemSpace"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireSystem::SystemData::* )( ::QList<SireMol::PartialMolecule> const & ) const)( &::SireSystem::SystemData::mapIntoSystemSpace )
            , ( bp::arg("molecules") ) )    
        .def( 
            "name"
            , (::QString const & ( ::SireSystem::SystemData::* )(  ) const)( &::SireSystem::SystemData::name )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "refersTo"
            , (bool ( ::SireSystem::SystemData::* )( ::SireMol::MoleculeGroupID ) const)( &::SireSystem::SystemData::refersTo )
            , ( bp::arg("groupid") ) )    
        .def( 
            "remove"
            , (void ( ::SireSystem::SystemData::* )( ::SireMol::MoleculeGroupID ) )( &::SireSystem::SystemData::remove )
            , ( bp::arg("groupid") ) )    
        .def( 
            "remove"
            , (void ( ::SireSystem::SystemData::* )( ::SireMol::MoleculeGroup const & ) )( &::SireSystem::SystemData::remove )
            , ( bp::arg("group") ) )    
        .def( 
            "remove"
            , (void ( ::SireSystem::SystemData::* )( ::SireMol::PartialMolecule const &,::QSet<SireMol::MoleculeGroupID> const & ) )( &::SireSystem::SystemData::remove )
            , ( bp::arg("molecule"), bp::arg("groupids") ) )    
        .def( 
            "remove"
            , (void ( ::SireSystem::SystemData::* )( ::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> const &,::QSet<SireMol::MoleculeGroupID> const & ) )( &::SireSystem::SystemData::remove )
            , ( bp::arg("molecules"), bp::arg("groupids") ) )    
        .def( 
            "remove"
            , (void ( ::SireSystem::SystemData::* )( ::QList<SireMol::PartialMolecule> const &,::QSet<SireMol::MoleculeGroupID> const & ) )( &::SireSystem::SystemData::remove )
            , ( bp::arg("molecules"), bp::arg("groupids") ) )    
        .def( 
            "remove"
            , (void ( ::SireSystem::SystemData::* )( ::SireMol::PartialMolecule const & ) )( &::SireSystem::SystemData::remove )
            , ( bp::arg("molecule") ) )    
        .def( 
            "remove"
            , (void ( ::SireSystem::SystemData::* )( ::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> const & ) )( &::SireSystem::SystemData::remove )
            , ( bp::arg("molecules") ) )    
        .def( 
            "remove"
            , (void ( ::SireSystem::SystemData::* )( ::QList<SireMol::PartialMolecule> const & ) )( &::SireSystem::SystemData::remove )
            , ( bp::arg("molecules") ) )    
        .def( 
            "space"
            , (::SireVol::Space const & ( ::SireSystem::SystemData::* )(  ) const)( &::SireSystem::SystemData::space )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "version"
            , (::SireBase::Version const & ( ::SireSystem::SystemData::* )(  ) const)( &::SireSystem::SystemData::version )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireSystem::SystemData >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireSystem::SystemData >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
