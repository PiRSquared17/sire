// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MoleculeGroups.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

SireMol::MoleculeGroups __copy__(const SireMol::MoleculeGroups &other){ return SireMol::MoleculeGroups(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMol::MoleculeGroups&){ return "SireMol::MoleculeGroups";}

void register_MoleculeGroups_class(){

    bp::class_< SireMol::MoleculeGroups >( "MoleculeGroups" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::MoleculeGroup const & >(( bp::arg("group") )) )    
        .def( 
            "add"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::SireMol::MoleculeGroup const & ) )( &::SireMol::MoleculeGroups::add )
            , ( bp::arg("group") ) )    
        .def( 
            "add"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::SireMol::PartialMolecule const &,::SireMol::MoleculeGroupID ) )( &::SireMol::MoleculeGroups::add )
            , ( bp::arg("molecule"), bp::arg("groupid") ) )    
        .def( 
            "add"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::QVector<SireMol::PartialMolecule> const &,::SireMol::MoleculeGroupID ) )( &::SireMol::MoleculeGroups::add )
            , ( bp::arg("molecules"), bp::arg("groupid") ) )    
        .def( 
            "assertContains"
            , &::SireMol::MoleculeGroups::assertContains
            , ( bp::arg("groupid") ) )    
        .def( 
            "change"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::SireMol::MoleculeGroup const & ) )( &::SireMol::MoleculeGroups::change )
            , ( bp::arg("group") ) )    
        .def( 
            "change"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::SireMol::PartialMolecule const & ) )( &::SireMol::MoleculeGroups::change )
            , ( bp::arg("molecule") ) )    
        .def( 
            "change"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::QVector<SireMol::PartialMolecule> const & ) )( &::SireMol::MoleculeGroups::change )
            , ( bp::arg("molecules") ) )    
        .def( 
            "change"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> const & ) )( &::SireMol::MoleculeGroups::change )
            , ( bp::arg("molecules") ) )    
        .def( 
            "clear"
            , &::SireMol::MoleculeGroups::clear )    
        .def( 
            "contains"
            , &::SireMol::MoleculeGroups::contains
            , ( bp::arg("molecule") ) )    
        .def( 
            "count"
            , &::SireMol::MoleculeGroups::count )    
        .def( 
            "group"
            , (::SireMol::MoleculeGroup const & ( ::SireMol::MoleculeGroups::* )( ::SireMol::MoleculeGroupID ) const)( &::SireMol::MoleculeGroups::group )
            , ( bp::arg("groupid") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "group"
            , (::SireMol::MoleculeGroup const & ( ::SireMol::MoleculeGroups::* )( ::QString const & ) const)( &::SireMol::MoleculeGroups::group )
            , ( bp::arg("name") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "groups"
            , (::QVector<SireMol::MoleculeGroup> ( ::SireMol::MoleculeGroups::* )(  ) const)( &::SireMol::MoleculeGroups::groups ) )    
        .def( 
            "groups"
            , (::QVector<SireMol::MoleculeGroup> ( ::SireMol::MoleculeGroups::* )( ::SireMol::MoleculeID ) const)( &::SireMol::MoleculeGroups::groups )
            , ( bp::arg("molid") ) )    
        .def( 
            "groupsContaining"
            , &::SireMol::MoleculeGroups::groupsContaining
            , ( bp::arg("molecule") ) )    
        .def( 
            "groupsReferringTo"
            , &::SireMol::MoleculeGroups::groupsReferringTo
            , ( bp::arg("molid") ) )    
        .def( 
            "moleculeIDs"
            , &::SireMol::MoleculeGroups::moleculeIDs )    
        .def( 
            "molecules"
            , &::SireMol::MoleculeGroups::molecules )    
        .def( 
            "nMolecules"
            , &::SireMol::MoleculeGroups::nMolecules )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "refersTo"
            , &::SireMol::MoleculeGroups::refersTo
            , ( bp::arg("molid") ) )    
        .def( 
            "remove"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::SireMol::MoleculeGroup const & ) )( &::SireMol::MoleculeGroups::remove )
            , ( bp::arg("group") ) )    
        .def( 
            "remove"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::SireMol::MoleculeGroupID ) )( &::SireMol::MoleculeGroups::remove )
            , ( bp::arg("groupid") ) )    
        .def( 
            "remove"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::QString const & ) )( &::SireMol::MoleculeGroups::remove )
            , ( bp::arg("groupname") ) )    
        .def( 
            "remove"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::SireMol::PartialMolecule const &,::SireMol::MoleculeGroupID ) )( &::SireMol::MoleculeGroups::remove )
            , ( bp::arg("molecule"), bp::arg("groupid") ) )    
        .def( 
            "remove"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::QVector<SireMol::PartialMolecule> const &,::SireMol::MoleculeGroupID ) )( &::SireMol::MoleculeGroups::remove )
            , ( bp::arg("molecules"), bp::arg("groupid") ) )    
        .def( 
            "remove"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::SireMol::PartialMolecule const & ) )( &::SireMol::MoleculeGroups::remove )
            , ( bp::arg("molecule") ) )    
        .def( 
            "remove"
            , (bool ( ::SireMol::MoleculeGroups::* )( ::QVector<SireMol::PartialMolecule> const & ) )( &::SireMol::MoleculeGroups::remove )
            , ( bp::arg("molecules") ) )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::MoleculeGroups >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::MoleculeGroups >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
