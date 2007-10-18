// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ForceField.pypp.hpp"
#include "boost/python.hpp"
#include "sireff_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
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

SireFF::ForceField __copy__(const SireFF::ForceField &other){ return SireFF::ForceField(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireFF::ForceField&){ return "SireFF::ForceField";}

void register_ForceField_class(){

    bp::class_< SireFF::ForceField >( "ForceField" )    
        .def( bp::init< >() )    
        .def( bp::init< SireFF::FFBase const & >(( bp::arg("ffbase") )) )    
        .def( bp::init< SireBase::SharedPolyPointer<SireFF::FFBase> const & >(( bp::arg("ffptr") )) )    
        .def( 
            "ID"
            , (::SireFF::ForceFieldID ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::ID ) )    
        .def( 
            "add"
            , (bool ( ::SireFF::ForceField::* )( ::SireMol::PartialMolecule const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceField::add )
            , ( bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "add"
            , (bool ( ::SireFF::ForceField::* )( ::SireMol::Molecules const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceField::add )
            , ( bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "addTo"
            , (bool ( ::SireFF::ForceField::* )( ::SireFF::FFBase::Group const &,::SireMol::PartialMolecule const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceField::addTo )
            , ( bp::arg("group"), bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "addTo"
            , (bool ( ::SireFF::ForceField::* )( ::SireFF::FFBase::Group const &,::SireMol::Molecules const &,::SireFF::ParameterMap const & ) )( &::SireFF::ForceField::addTo )
            , ( bp::arg("group"), bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) )    
        .def( 
            "assertContains"
            , (void ( ::SireFF::ForceField::* )( ::SireFF::FFComponent const & ) const)( &::SireFF::ForceField::assertContains )
            , ( bp::arg("component") ) )    
        .def( 
            "base"
            , (::SireFF::FFBase const & ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::base )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "change"
            , (bool ( ::SireFF::ForceField::* )( ::SireMol::PartialMolecule const & ) )( &::SireFF::ForceField::change )
            , ( bp::arg("molecule") ) )    
        .def( 
            "change"
            , (bool ( ::SireFF::ForceField::* )( ::SireMol::Molecules const & ) )( &::SireFF::ForceField::change )
            , ( bp::arg("molecules") ) )    
        .def( 
            "components"
            , (::SireFF::FFBase::Components const & ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::components )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "contains"
            , (bool ( ::SireFF::ForceField::* )( ::SireMol::PartialMolecule const & ) const)( &::SireFF::ForceField::contains )
            , ( bp::arg("molecule") ) )    
        .def( 
            "contains"
            , (bool ( ::SireFF::ForceField::* )( ::SireMol::PartialMolecule const &,::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceField::contains )
            , ( bp::arg("molecule"), bp::arg("group") ) )    
        .def( 
            "containsProperty"
            , (bool ( ::SireFF::ForceField::* )( ::QString const & ) const)( &::SireFF::ForceField::containsProperty )
            , ( bp::arg("name") ) )    
        .def( 
            "contents"
            , (::SireMol::Molecules ( ::SireFF::ForceField::* )( ::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceField::contents )
            , ( bp::arg("group") ) )    
        .def( 
            "contents"
            , (::SireMol::Molecules ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::contents ) )    
        .def( 
            "energies"
            , (::SireCAS::Values ( ::SireFF::ForceField::* )(  ) )( &::SireFF::ForceField::energies ) )    
        .def( 
            "energies"
            , (::SireCAS::Values ( ::SireFF::ForceField::* )( ::QSet<SireFF::FFComponent> const & ) )( &::SireFF::ForceField::energies )
            , ( bp::arg("components") ) )    
        .def( 
            "energy"
            , (double ( ::SireFF::ForceField::* )(  ) )( &::SireFF::ForceField::energy ) )    
        .def( 
            "energy"
            , (double ( ::SireFF::ForceField::* )( ::SireFF::FFComponent const & ) )( &::SireFF::ForceField::energy )
            , ( bp::arg("component") ) )    
        .def( 
            "energy"
            , (double ( ::SireFF::ForceField::* )( ::SireCAS::Function const & ) )( &::SireFF::ForceField::energy )
            , ( bp::arg("component") ) )    
        .def( 
            "energy"
            , (double ( ::SireFF::ForceField::* )( ::SireCAS::Symbol const & ) )( &::SireFF::ForceField::energy )
            , ( bp::arg("component") ) )    
        .def( 
            "getProperty"
            , (::SireBase::Property ( ::SireFF::ForceField::* )( ::QString const & ) const)( &::SireFF::ForceField::getProperty )
            , ( bp::arg("name") ) )    
        .def( 
            "groups"
            , (::SireFF::FFBase::Groups const & ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::groups )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "groupsReferringTo"
            , (::QSet<SireFF::FFBase::Group> ( ::SireFF::ForceField::* )( ::SireMol::MoleculeID ) const)( &::SireFF::ForceField::groupsReferringTo )
            , ( bp::arg("molid") ) )    
        .def( 
            "isClean"
            , (bool ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::isClean ) )    
        .def( 
            "isDirty"
            , (bool ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::isDirty ) )    
        .def( 
            "molecule"
            , (::SireMol::PartialMolecule ( ::SireFF::ForceField::* )( ::SireMol::MoleculeID ) const)( &::SireFF::ForceField::molecule )
            , ( bp::arg("molid") ) )    
        .def( 
            "molecule"
            , (::SireMol::PartialMolecule ( ::SireFF::ForceField::* )( ::SireMol::MoleculeID,::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceField::molecule )
            , ( bp::arg("molid"), bp::arg("group") ) )    
        .def( 
            "moleculeIDs"
            , (::QSet<SireMol::MoleculeID> ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::moleculeIDs ) )    
        .def( 
            "moleculeIDs"
            , (::QSet<SireMol::MoleculeID> ( ::SireFF::ForceField::* )( ::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceField::moleculeIDs )
            , ( bp::arg("group") ) )    
        .def( 
            "molecules"
            , (::SireMol::Molecules ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::molecules ) )    
        .def( 
            "molecules"
            , (::SireMol::Molecules ( ::SireFF::ForceField::* )( ::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceField::molecules )
            , ( bp::arg("group") ) )    
        .def( 
            "molecules"
            , (::SireMol::Molecules ( ::SireFF::ForceField::* )( ::QSet<SireMol::MoleculeID> const & ) const)( &::SireFF::ForceField::molecules )
            , ( bp::arg("molids") ) )    
        .def( 
            "mustNowRecalculateFromScratch"
            , (void ( ::SireFF::ForceField::* )(  ) )( &::SireFF::ForceField::mustNowRecalculateFromScratch ) )    
        .def( 
            "name"
            , (::QString const & ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::name )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "parameters"
            , (::SireFF::FFBase::Parameters const & ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::parameters )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "properties"
            , (::QHash<QString,SireBase::Property> ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::properties ) )    
        .def( 
            "refersTo"
            , (bool ( ::SireFF::ForceField::* )( ::SireMol::MoleculeID ) const)( &::SireFF::ForceField::refersTo )
            , ( bp::arg("molid") ) )    
        .def( 
            "refersTo"
            , (bool ( ::SireFF::ForceField::* )( ::SireMol::MoleculeID,::SireFF::FFBase::Group const & ) const)( &::SireFF::ForceField::refersTo )
            , ( bp::arg("molid"), bp::arg("group") ) )    
        .def( 
            "remove"
            , (bool ( ::SireFF::ForceField::* )( ::SireMol::PartialMolecule const & ) )( &::SireFF::ForceField::remove )
            , ( bp::arg("molecule") ) )    
        .def( 
            "remove"
            , (bool ( ::SireFF::ForceField::* )( ::SireMol::Molecules const & ) )( &::SireFF::ForceField::remove )
            , ( bp::arg("molecules") ) )    
        .def( 
            "removeFrom"
            , (bool ( ::SireFF::ForceField::* )( ::SireFF::FFBase::Group const &,::SireMol::PartialMolecule const & ) )( &::SireFF::ForceField::removeFrom )
            , ( bp::arg("group"), bp::arg("molecule") ) )    
        .def( 
            "removeFrom"
            , (bool ( ::SireFF::ForceField::* )( ::SireFF::FFBase::Group const &,::SireMol::Molecules const & ) )( &::SireFF::ForceField::removeFrom )
            , ( bp::arg("group"), bp::arg("molecules") ) )    
        .def( 
            "setName"
            , (void ( ::SireFF::ForceField::* )( ::QString const & ) )( &::SireFF::ForceField::setName )
            , ( bp::arg("name") ) )    
        .def( 
            "setProperty"
            , (bool ( ::SireFF::ForceField::* )( ::QString const &,::SireBase::Property const & ) )( &::SireFF::ForceField::setProperty )
            , ( bp::arg("name"), bp::arg("value") ) )    
        .def( 
            "version"
            , (::SireBase::Version const & ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::version )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "what"
            , (char const * ( ::SireFF::ForceField::* )(  ) const)( &::SireFF::ForceField::what ) )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireFF::ForceField >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireFF::ForceField >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
