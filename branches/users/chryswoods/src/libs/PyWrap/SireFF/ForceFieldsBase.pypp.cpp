// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ForceFieldsBase.pypp.hpp"
#include "boost/python.hpp"
#include "sireff_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"
#include "SireMol/cgatomid.h"
#include "SireMol/cutgroupnum.h"
#include "SireMol/resnumatomid.h"
#include "SireMol/resid.h"
#include "SireMol/moleculeid.h"
#include "SireBase/property.h"

namespace bp = boost::python;

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireFF::ForceFieldsBase&){ return "SireFF::ForceFieldsBase";}

void register_ForceFieldsBase_class(){

    bp::class_< SireFF::ForceFieldsBase, boost::noncopyable >( "ForceFieldsBase", bp::no_init )    
        .def( 
            "add"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::SireMol::PartialMolecule const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceFieldsBase::add )
            , ( bp::arg("ffid"), bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "add"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::QList<SireMol::PartialMolecule> const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceFieldsBase::add )
            , ( bp::arg("ffid"), bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "add"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::ForceFieldID> const &,::SireMol::PartialMolecule const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceFieldsBase::add )
            , ( bp::arg("ffids"), bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "add"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::ForceFieldID> const &,::QList<SireMol::PartialMolecule> const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceFieldsBase::add )
            , ( bp::arg("ffids"), bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "addTo"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::SireFF::FFBase::Group const &,::SireMol::PartialMolecule const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceFieldsBase::addTo )
            , ( bp::arg("ffid"), bp::arg("group"), bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "addTo"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::SireFF::FFBase::Group const &,::QList<SireMol::PartialMolecule> const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceFieldsBase::addTo )
            , ( bp::arg("ffid"), bp::arg("group"), bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "addTo"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::FFGroupID const &,::SireMol::PartialMolecule const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceFieldsBase::addTo )
            , ( bp::arg("ffgroupid"), bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "addTo"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::FFGroupID const &,::QList<SireMol::PartialMolecule> const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceFieldsBase::addTo )
            , ( bp::arg("ffgroupid"), bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "addTo"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::FFGroupID> const &,::SireMol::PartialMolecule const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceFieldsBase::addTo )
            , ( bp::arg("ffgroupids"), bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "addTo"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::FFGroupID> const &,::QList<SireMol::PartialMolecule> const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceFieldsBase::addTo )
            , ( bp::arg("ffgroupids"), bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "assertSameContents"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldsBase const & ) const)( &::SireFF::ForceFieldsBase::assertSameContents )
            , ( bp::arg("other") ) )    
        .def( 
            "assertValidComponents"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::FFExpression const & ) const)( &::SireFF::ForceFieldsBase::assertValidComponents )
            , ( bp::arg("expression") ) )    
        .def( 
            "change"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireMol::PartialMolecule const & ) )( &::SireFF::ForceFieldsBase::change )
            , ( bp::arg("molecule") ) )    
        .def( 
            "change"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QList<SireMol::PartialMolecule> const & ) )( &::SireFF::ForceFieldsBase::change )
            , ( bp::arg("molecules") ) )    
        .def( 
            "change"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> const & ) )( &::SireFF::ForceFieldsBase::change )
            , ( bp::arg("molecules") ) )    
        .def( 
            "contains"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireCAS::Function const & ) const)( &::SireFF::ForceFieldsBase::contains )
            , ( bp::arg("function") ) )    
        .def( 
            "contains"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID const ) const)( &::SireFF::ForceFieldsBase::contains )
            , ( bp::arg("ffid") ) )    
        .def( 
            "contains"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::PartialMolecule const & ) const)( &::SireFF::ForceFieldsBase::contains )
            , ( bp::arg("molecule") ) )    
        .def( 
            "contains"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::PartialMolecule const &,::SireFF::ForceFieldID ) const)( &::SireFF::ForceFieldsBase::contains )
            , ( bp::arg("molecule"), bp::arg("ffid") ) )    
        .def( 
            "contains"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::PartialMolecule const &,::QSet<SireFF::ForceFieldID> const & ) const)( &::SireFF::ForceFieldsBase::contains )
            , ( bp::arg("molecule"), bp::arg("ffids") ) )    
        .def( 
            "contains"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::PartialMolecule const &,::SireFF::ForceFieldID,::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceFieldsBase::contains )
            , ( bp::arg("molecule"), bp::arg("ffid"), bp::arg("group") ) )    
        .def( 
            "contains"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::PartialMolecule const &,::SireFF::FFGroupID const & ) const)( &::SireFF::ForceFieldsBase::contains )
            , ( bp::arg("molecule"), bp::arg("ffgroupid") ) )    
        .def( 
            "contains"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::PartialMolecule const &,::QSet<SireFF::FFGroupID> const & ) const)( &::SireFF::ForceFieldsBase::contains )
            , ( bp::arg("molecule"), bp::arg("ffgroupids") ) )    
        .def( 
            "containsProperty"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::QString const & ) const)( &::SireFF::ForceFieldsBase::containsProperty )
            , ( bp::arg("name") ) )    
        .def( 
            "containsProperty"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::QString const & ) const)( &::SireFF::ForceFieldsBase::containsProperty )
            , ( bp::arg("ffid"), bp::arg("name") ) )    
        .def( 
            "contents"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )(  ) const)( &::SireFF::ForceFieldsBase::contents ) )    
        .def( 
            "contents"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID ) const)( &::SireFF::ForceFieldsBase::contents )
            , ( bp::arg("ffid") ) )    
        .def( 
            "contents"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::ForceFieldID> const & ) const)( &::SireFF::ForceFieldsBase::contents )
            , ( bp::arg("ffids") ) )    
        .def( 
            "contents"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceFieldsBase::contents )
            , ( bp::arg("ffid"), bp::arg("group") ) )    
        .def( 
            "contents"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )( ::SireFF::FFGroupID const & ) const)( &::SireFF::ForceFieldsBase::contents )
            , ( bp::arg("ffgroupid") ) )    
        .def( 
            "contents"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::FFGroupID> const & ) const)( &::SireFF::ForceFieldsBase::contents )
            , ( bp::arg("ffgroupids") ) )    
        .def( 
            "energies"
            , (::SireCAS::Values ( ::SireFF::ForceFieldsBase::* )(  ) )( &::SireFF::ForceFieldsBase::energies ) )    
        .def( 
            "energies"
            , (::SireCAS::Values ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::FFComponent> const & ) )( &::SireFF::ForceFieldsBase::energies )
            , ( bp::arg("components") ) )    
        .def( 
            "energies"
            , (::SireCAS::Values ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireCAS::Function> const & ) )( &::SireFF::ForceFieldsBase::energies )
            , ( bp::arg("components") ) )    
        .def( 
            "energy"
            , (double ( ::SireFF::ForceFieldsBase::* )( ::SireCAS::Expression const & ) )( &::SireFF::ForceFieldsBase::energy )
            , ( bp::arg("expression") ) )    
        .def( 
            "energy"
            , (double ( ::SireFF::ForceFieldsBase::* )( ::SireFF::FFExpression const & ) )( &::SireFF::ForceFieldsBase::energy )
            , ( bp::arg("expression") ) )    
        .def( 
            "energy"
            , (double ( ::SireFF::ForceFieldsBase::* )( ::SireCAS::Function const & ) )( &::SireFF::ForceFieldsBase::energy )
            , ( bp::arg("component") ) )    
        .def( 
            "energy"
            , (double ( ::SireFF::ForceFieldsBase::* )( ::SireFF::FFComponent const & ) )( &::SireFF::ForceFieldsBase::energy )
            , ( bp::arg("component") ) )    
        .def( 
            "energy"
            , (double ( ::SireFF::ForceFieldsBase::* )(  ) )( &::SireFF::ForceFieldsBase::energy ) )    
        .def( 
            "expression"
            , (::SireFF::FFExpression const & ( ::SireFF::ForceFieldsBase::* )( ::SireCAS::Function const & ) const)( &::SireFF::ForceFieldsBase::expression )
            , ( bp::arg("function") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "expressions"
            , (::QHash<SireCAS::Function,SireFF::FFExpression> ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireCAS::Function> const & ) const)( &::SireFF::ForceFieldsBase::expressions )
            , ( bp::arg("functions") ) )    
        .def( 
            "expressions"
            , (::QHash<SireCAS::Function,SireFF::FFExpression> ( ::SireFF::ForceFieldsBase::* )(  ) const)( &::SireFF::ForceFieldsBase::expressions ) )    
        .def( 
            "forceField"
            , (::SireFF::ForceField ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID ) const)( &::SireFF::ForceFieldsBase::forceField )
            , ( bp::arg("ffid") ) )    
        .def( 
            "forceFieldGroupsContaining"
            , (::QSet<SireFF::FFGroupID> ( ::SireFF::ForceFieldsBase::* )( ::SireMol::PartialMolecule const & ) const)( &::SireFF::ForceFieldsBase::forceFieldGroupsContaining )
            , ( bp::arg("molecule") ) )    
        .def( 
            "forceFieldGroupsReferringTo"
            , (::QSet<SireFF::FFGroupID> ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID ) const)( &::SireFF::ForceFieldsBase::forceFieldGroupsReferringTo )
            , ( bp::arg("molid") ) )    
        .def( 
            "forceFieldIDs"
            , (::QSet<SireFF::ForceFieldID> ( ::SireFF::ForceFieldsBase::* )(  ) const)( &::SireFF::ForceFieldsBase::forceFieldIDs ) )    
        .def( 
            "forceFields"
            , (::QHash<SireFF::ForceFieldID, SireFF::ForceField> ( ::SireFF::ForceFieldsBase::* )(  ) const)( &::SireFF::ForceFieldsBase::forceFields ) )    
        .def( 
            "forceFields"
            , (::QHash<SireFF::ForceFieldID, SireFF::ForceField> ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::ForceFieldID> const & ) const)( &::SireFF::ForceFieldsBase::forceFields )
            , ( bp::arg("ffids") ) )    
        .def( 
            "forceFieldsContaining"
            , (::QSet<SireFF::ForceFieldID> ( ::SireFF::ForceFieldsBase::* )( ::SireMol::PartialMolecule const & ) const)( &::SireFF::ForceFieldsBase::forceFieldsContaining )
            , ( bp::arg("molecule") ) )    
        .def( 
            "forceFieldsReferringTo"
            , (::QSet<SireFF::ForceFieldID> ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID ) const)( &::SireFF::ForceFieldsBase::forceFieldsReferringTo )
            , ( bp::arg("molid") ) )    
        .def( 
            "forceFieldsWithProperty"
            , (::QSet<SireFF::ForceFieldID> ( ::SireFF::ForceFieldsBase::* )( ::QString const & ) const)( &::SireFF::ForceFieldsBase::forceFieldsWithProperty )
            , ( bp::arg("name") ) )    
        .def( 
            "getProperty"
            , (::QHash<SireFF::ForceFieldID,SireBase::Property> ( ::SireFF::ForceFieldsBase::* )( ::QString const & ) const)( &::SireFF::ForceFieldsBase::getProperty )
            , ( bp::arg("name") ) )    
        .def( 
            "getProperty"
            , (::SireBase::Property ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::QString const & ) const)( &::SireFF::ForceFieldsBase::getProperty )
            , ( bp::arg("ffid"), bp::arg("name") ) )    
        .def( 
            "getProperty"
            , (::QHash<SireFF::ForceFieldID,SireBase::Property> ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::ForceFieldID> const &,::QString const & ) const)( &::SireFF::ForceFieldsBase::getProperty )
            , ( bp::arg("ffids"), bp::arg("name") ) )    
        .def( 
            "molecule"
            , (::SireMol::PartialMolecule ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID ) const)( &::SireFF::ForceFieldsBase::molecule )
            , ( bp::arg("molid") ) )    
        .def( 
            "molecule"
            , (::SireMol::PartialMolecule ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID,::SireFF::ForceFieldID ) const)( &::SireFF::ForceFieldsBase::molecule )
            , ( bp::arg("molid"), bp::arg("ffid") ) )    
        .def( 
            "molecule"
            , (::SireMol::PartialMolecule ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID,::QSet<SireFF::ForceFieldID> const & ) const)( &::SireFF::ForceFieldsBase::molecule )
            , ( bp::arg("molid"), bp::arg("ffids") ) )    
        .def( 
            "molecule"
            , (::SireMol::PartialMolecule ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID,::SireFF::ForceFieldID,::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceFieldsBase::molecule )
            , ( bp::arg("molid"), bp::arg("ffid"), bp::arg("group") ) )    
        .def( 
            "molecule"
            , (::SireMol::PartialMolecule ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID,::SireFF::FFGroupID const & ) const)( &::SireFF::ForceFieldsBase::molecule )
            , ( bp::arg("molid"), bp::arg("ffgroupid") ) )    
        .def( 
            "molecule"
            , (::SireMol::PartialMolecule ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID,::QSet<SireFF::FFGroupID> const & ) const)( &::SireFF::ForceFieldsBase::molecule )
            , ( bp::arg("molid"), bp::arg("ffgroupids") ) )    
        .def( 
            "moleculeIDs"
            , (::QSet<SireMol::MoleculeID> ( ::SireFF::ForceFieldsBase::* )(  ) const)( &::SireFF::ForceFieldsBase::moleculeIDs ) )    
        .def( 
            "moleculeIDs"
            , (::QSet<SireMol::MoleculeID> ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID ) const)( &::SireFF::ForceFieldsBase::moleculeIDs )
            , ( bp::arg("ffid") ) )    
        .def( 
            "moleculeIDs"
            , (::QSet<SireMol::MoleculeID> ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceFieldsBase::moleculeIDs )
            , ( bp::arg("ffid"), bp::arg("group") ) )    
        .def( 
            "moleculeIDs"
            , (::QSet<SireMol::MoleculeID> ( ::SireFF::ForceFieldsBase::* )( ::SireFF::FFGroupID const & ) const)( &::SireFF::ForceFieldsBase::moleculeIDs )
            , ( bp::arg("ffgroupid") ) )    
        .def( 
            "molecules"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )(  ) const)( &::SireFF::ForceFieldsBase::molecules ) )    
        .def( 
            "molecules"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID ) const)( &::SireFF::ForceFieldsBase::molecules )
            , ( bp::arg("ffid") ) )    
        .def( 
            "molecules"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::ForceFieldID> const & ) const)( &::SireFF::ForceFieldsBase::molecules )
            , ( bp::arg("ffids") ) )    
        .def( 
            "molecules"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceFieldsBase::molecules )
            , ( bp::arg("ffid"), bp::arg("group") ) )    
        .def( 
            "molecules"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )( ::SireFF::FFGroupID const & ) const)( &::SireFF::ForceFieldsBase::molecules )
            , ( bp::arg("ffgroupid") ) )    
        .def( 
            "molecules"
            , (::QHash<SireMol::MoleculeID,SireMol::PartialMolecule> ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::FFGroupID> const & ) const)( &::SireFF::ForceFieldsBase::molecules )
            , ( bp::arg("ffgroupids") ) )    
        .def( 
            "nForceFields"
            , (int ( ::SireFF::ForceFieldsBase::* )(  ) const)( &::SireFF::ForceFieldsBase::nForceFields ) )    
        .def( 
            "nMolecules"
            , (int ( ::SireFF::ForceFieldsBase::* )(  ) const)( &::SireFF::ForceFieldsBase::nMolecules ) )    
        .def( 
            "parameters"
            , (::SireCAS::Values const & ( ::SireFF::ForceFieldsBase::* )(  ) const)( &::SireFF::ForceFieldsBase::parameters )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "properties"
            , (::QHash<QString,QHash<SireFF::ForceFieldID, SireBase::Property> > ( ::SireFF::ForceFieldsBase::* )(  ) const)( &::SireFF::ForceFieldsBase::properties ) )    
        .def( 
            "properties"
            , (::QHash<QString,QHash<SireFF::ForceFieldID, SireBase::Property> > ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::ForceFieldID> const & ) const)( &::SireFF::ForceFieldsBase::properties )
            , ( bp::arg("ffids") ) )    
        .def( 
            "refersTo"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID ) const)( &::SireFF::ForceFieldsBase::refersTo )
            , ( bp::arg("molid") ) )    
        .def( 
            "refersTo"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID,::SireFF::ForceFieldID ) const)( &::SireFF::ForceFieldsBase::refersTo )
            , ( bp::arg("molid"), bp::arg("ffid") ) )    
        .def( 
            "refersTo"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID,::QSet<SireFF::ForceFieldID> const & ) const)( &::SireFF::ForceFieldsBase::refersTo )
            , ( bp::arg("molid"), bp::arg("ffids") ) )    
        .def( 
            "refersTo"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID,::SireFF::ForceFieldID,::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceFieldsBase::refersTo )
            , ( bp::arg("molid"), bp::arg("ffid"), bp::arg("group") ) )    
        .def( 
            "refersTo"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID,::SireFF::FFGroupID const & ) const)( &::SireFF::ForceFieldsBase::refersTo )
            , ( bp::arg("molid"), bp::arg("ffgroupid") ) )    
        .def( 
            "refersTo"
            , (bool ( ::SireFF::ForceFieldsBase::* )( ::SireMol::MoleculeID,::QSet<SireFF::FFGroupID> const & ) const)( &::SireFF::ForceFieldsBase::refersTo )
            , ( bp::arg("molid"), bp::arg("ffgroupids") ) )    
        .def( 
            "remove"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireMol::PartialMolecule const & ) )( &::SireFF::ForceFieldsBase::remove )
            , ( bp::arg("molecule") ) )    
        .def( 
            "remove"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QList<SireMol::PartialMolecule> const & ) )( &::SireFF::ForceFieldsBase::remove )
            , ( bp::arg("molecules") ) )    
        .def( 
            "remove"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::SireMol::PartialMolecule const & ) )( &::SireFF::ForceFieldsBase::remove )
            , ( bp::arg("ffid"), bp::arg("molecule") ) )    
        .def( 
            "remove"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::QList<SireMol::PartialMolecule> const & ) )( &::SireFF::ForceFieldsBase::remove )
            , ( bp::arg("ffid"), bp::arg("molecules") ) )    
        .def( 
            "remove"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::ForceFieldID> const &,::SireMol::PartialMolecule const & ) )( &::SireFF::ForceFieldsBase::remove )
            , ( bp::arg("ffids"), bp::arg("molecule") ) )    
        .def( 
            "remove"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::ForceFieldID> const &,::QList<SireMol::PartialMolecule> const & ) )( &::SireFF::ForceFieldsBase::remove )
            , ( bp::arg("ffids"), bp::arg("molecules") ) )    
        .def( 
            "removeFrom"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::SireFF::FFBase::Group const &,::SireMol::PartialMolecule const & ) )( &::SireFF::ForceFieldsBase::removeFrom )
            , ( bp::arg("ffid"), bp::arg("group"), bp::arg("molecule") ) )    
        .def( 
            "removeFrom"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::SireFF::FFBase::Group const &,::QList<SireMol::PartialMolecule> const & ) )( &::SireFF::ForceFieldsBase::removeFrom )
            , ( bp::arg("ffid"), bp::arg("group"), bp::arg("molecules") ) )    
        .def( 
            "removeFrom"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::FFGroupID const &,::SireMol::PartialMolecule const & ) )( &::SireFF::ForceFieldsBase::removeFrom )
            , ( bp::arg("ffgroupid"), bp::arg("molecule") ) )    
        .def( 
            "removeFrom"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::FFGroupID const &,::QList<SireMol::PartialMolecule> const & ) )( &::SireFF::ForceFieldsBase::removeFrom )
            , ( bp::arg("ffgroupid"), bp::arg("molecules") ) )    
        .def( 
            "removeFrom"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::FFGroupID> const &,::SireMol::PartialMolecule const & ) )( &::SireFF::ForceFieldsBase::removeFrom )
            , ( bp::arg("ffgroupids"), bp::arg("molecule") ) )    
        .def( 
            "removeFrom"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::FFGroupID> const &,::QList<SireMol::PartialMolecule> const & ) )( &::SireFF::ForceFieldsBase::removeFrom )
            , ( bp::arg("ffgroupids"), bp::arg("molecules") ) )    
        .def( 
            "setParameter"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireCAS::Symbol const &,double ) )( &::SireFF::ForceFieldsBase::setParameter )
            , ( bp::arg("param"), bp::arg("value") ) )    
        .def( 
            "setProperty"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QString const &,::SireBase::Property const & ) )( &::SireFF::ForceFieldsBase::setProperty )
            , ( bp::arg("name"), bp::arg("property") ) )    
        .def( 
            "setProperty"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::ForceFieldID,::QString const &,::SireBase::Property const & ) )( &::SireFF::ForceFieldsBase::setProperty )
            , ( bp::arg("ffid"), bp::arg("name"), bp::arg("property") ) )    
        .def( 
            "setProperty"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::QSet<SireFF::ForceFieldID> const &,::QString const &,::SireBase::Property const & ) )( &::SireFF::ForceFieldsBase::setProperty )
            , ( bp::arg("ffids"), bp::arg("name"), bp::arg("property") ) )    
        .def( 
            "setTotal"
            , (void ( ::SireFF::ForceFieldsBase::* )( ::SireFF::FFExpression const & ) )( &::SireFF::ForceFieldsBase::setTotal )
            , ( bp::arg("expression") ) )    
        .def( 
            "total"
            , (::SireFF::FFExpression const & ( ::SireFF::ForceFieldsBase::* )(  ) const)( &::SireFF::ForceFieldsBase::total )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireFF::ForceFieldsBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireFF::ForceFieldsBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
