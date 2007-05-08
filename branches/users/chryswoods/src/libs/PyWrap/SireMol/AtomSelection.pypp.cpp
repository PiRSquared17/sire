// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "AtomSelection.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMol/moleculedata.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"
#include "SireVol/space.h"

namespace bp = boost::python;

SireMol::AtomSelection __copy__(const SireMol::AtomSelection &other){ return SireMol::AtomSelection(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMol::AtomSelection&){ return "SireMol::AtomSelection";}

void register_AtomSelection_class(){

    bp::class_< SireMol::AtomSelection >( "AtomSelection" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molecule") )) )    
        .def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") )) )    
        .def( 
            "applyMask"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::QSet<SireMol::CutGroupID> const & ) const)( &::SireMol::AtomSelection::applyMask )
            , ( bp::arg("cgids") ) )    
        .def( 
            "applyMask"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::QSet<SireMol::ResNum> const & ) const)( &::SireMol::AtomSelection::applyMask )
            , ( bp::arg("resnums") ) )    
        .def( 
            "applyMask"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::applyMask )
            , ( bp::arg("selection") ) )    
        .def( 
            "applyMask"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::applyMask )
            , ( bp::arg("selection") ) )    
        .def( 
            "asSingleAtom"
            , &::SireMol::AtomSelection::asSingleAtom )    
        .def( 
            "asSingleCutGroup"
            , &::SireMol::AtomSelection::asSingleCutGroup )    
        .def( 
            "asSingleResidue"
            , &::SireMol::AtomSelection::asSingleResidue )    
        .def( 
            "assertCompatibleWith"
            , (void ( ::SireMol::AtomSelection::* )( ::SireMol::MoleculeInfo const & ) const)( &::SireMol::AtomSelection::assertCompatibleWith )
            , ( bp::arg("molinfo") ) )    
        .def( 
            "assertCompatibleWith"
            , (void ( ::SireMol::AtomSelection::* )( ::SireMol::Molecule const & ) const)( &::SireMol::AtomSelection::assertCompatibleWith )
            , ( bp::arg("molecule") ) )    
        .def( 
            "assertCompatibleWith"
            , (void ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::assertCompatibleWith )
            , ( bp::arg("other") ) )    
        .def( 
            "assertCompatibleWith"
            , (void ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::assertCompatibleWith )
            , ( bp::arg("other") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::contains )
            , ( bp::arg("selection") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::contains )
            , ( bp::arg("selection") ) )    
        .def( 
            "deselect"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::CutGroupID ) const)( &::SireMol::AtomSelection::deselect )
            , ( bp::arg("cgid") ) )    
        .def( 
            "deselect"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::ResNum ) const)( &::SireMol::AtomSelection::deselect )
            , ( bp::arg("resnum") ) )    
        .def( 
            "deselect"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::CGAtomID const & ) const)( &::SireMol::AtomSelection::deselect )
            , ( bp::arg("cgatomid") ) )    
        .def( 
            "deselect"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::IDMolAtom const & ) const)( &::SireMol::AtomSelection::deselect )
            , ( bp::arg("atomid") ) )    
        .def( 
            "deselect"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::deselect )
            , ( bp::arg("selection") ) )    
        .def( 
            "deselect"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::deselect )
            , ( bp::arg("selection") ) )    
        .def( 
            "deselectAll"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )(  ) const)( &::SireMol::AtomSelection::deselectAll ) )    
        .def( 
            "deselectAll"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::CutGroupID ) const)( &::SireMol::AtomSelection::deselectAll )
            , ( bp::arg("cgid") ) )    
        .def( 
            "deselectAll"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::ResNum ) const)( &::SireMol::AtomSelection::deselectAll )
            , ( bp::arg("resnum") ) )    
        .def( 
            "deselectAll"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::deselectAll )
            , ( bp::arg("selection") ) )    
        .def( 
            "deselectAll"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::deselectAll )
            , ( bp::arg("selection") ) )    
        .def( 
            "info"
            , &::SireMol::AtomSelection::info
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "intersect"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::intersect )
            , ( bp::arg("selection") ) )    
        .def( 
            "intersect"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::intersect )
            , ( bp::arg("selection") ) )    
        .def( 
            "intersects"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::intersects )
            , ( bp::arg("selection") ) )    
        .def( 
            "intersects"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::intersects )
            , ( bp::arg("selection") ) )    
        .def( 
            "invert"
            , &::SireMol::AtomSelection::invert )    
        .def( 
            "isEmpty"
            , &::SireMol::AtomSelection::isEmpty )    
        .def( 
            "isSingleAtom"
            , &::SireMol::AtomSelection::isSingleAtom )    
        .def( 
            "isSingleCutGroup"
            , &::SireMol::AtomSelection::isSingleCutGroup )    
        .def( 
            "isSingleResidue"
            , &::SireMol::AtomSelection::isSingleResidue )    
        .def( 
            "nSelected"
            , (int ( ::SireMol::AtomSelection::* )(  ) const)( &::SireMol::AtomSelection::nSelected ) )    
        .def( 
            "nSelected"
            , (int ( ::SireMol::AtomSelection::* )( ::SireMol::CutGroupID ) const)( &::SireMol::AtomSelection::nSelected )
            , ( bp::arg("cgid") ) )    
        .def( 
            "nSelected"
            , (int ( ::SireMol::AtomSelection::* )( ::SireMol::ResNum ) const)( &::SireMol::AtomSelection::nSelected )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nSelected"
            , (int ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::nSelected )
            , ( bp::arg("selection") ) )    
        .def( 
            "nSelected"
            , (int ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::nSelected )
            , ( bp::arg("selection") ) )    
        .def( 
            "nSelectedCutGroups"
            , &::SireMol::AtomSelection::nSelectedCutGroups )    
        .def( 
            "nSelectedResidues"
            , &::SireMol::AtomSelection::nSelectedResidues )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "select"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::CutGroupID ) const)( &::SireMol::AtomSelection::select )
            , ( bp::arg("cgid") ) )    
        .def( 
            "select"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::ResNum ) const)( &::SireMol::AtomSelection::select )
            , ( bp::arg("resnum") ) )    
        .def( 
            "select"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::CGAtomID const & ) const)( &::SireMol::AtomSelection::select )
            , ( bp::arg("cgatomid") ) )    
        .def( 
            "select"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::IDMolAtom const & ) const)( &::SireMol::AtomSelection::select )
            , ( bp::arg("atomid") ) )    
        .def( 
            "select"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::select )
            , ( bp::arg("selection") ) )    
        .def( 
            "select"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::select )
            , ( bp::arg("selection") ) )    
        .def( 
            "selectAll"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )(  ) const)( &::SireMol::AtomSelection::selectAll ) )    
        .def( 
            "selectAll"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::CutGroupID ) const)( &::SireMol::AtomSelection::selectAll )
            , ( bp::arg("cgid") ) )    
        .def( 
            "selectAll"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::ResNum ) const)( &::SireMol::AtomSelection::selectAll )
            , ( bp::arg("resnum") ) )    
        .def( 
            "selectAll"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::selectAll )
            , ( bp::arg("selection") ) )    
        .def( 
            "selectAll"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::selectAll )
            , ( bp::arg("selection") ) )    
        .def( 
            "selectNone"
            , &::SireMol::AtomSelection::selectNone )    
        .def( 
            "selected"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::CGAtomID const & ) const)( &::SireMol::AtomSelection::selected )
            , ( bp::arg("cgatomid") ) )    
        .def( 
            "selected"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::IDMolAtom const & ) const)( &::SireMol::AtomSelection::selected )
            , ( bp::arg("atomid") ) )    
        .def( 
            "selectedAll"
            , (bool ( ::SireMol::AtomSelection::* )(  ) const)( &::SireMol::AtomSelection::selectedAll ) )    
        .def( 
            "selectedAll"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::CutGroupID ) const)( &::SireMol::AtomSelection::selectedAll )
            , ( bp::arg("cgid") ) )    
        .def( 
            "selectedAll"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::ResNum ) const)( &::SireMol::AtomSelection::selectedAll )
            , ( bp::arg("resnum") ) )    
        .def( 
            "selectedAll"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::selectedAll )
            , ( bp::arg("selection") ) )    
        .def( 
            "selectedAll"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::selectedAll )
            , ( bp::arg("selection") ) )    
        .def( 
            "selectedAllCutGroups"
            , &::SireMol::AtomSelection::selectedAllCutGroups )    
        .def( 
            "selectedAllResidues"
            , &::SireMol::AtomSelection::selectedAllResidues )    
        .def( 
            "selectedAtoms"
            , (::QSet<SireMol::AtomIndex> ( ::SireMol::AtomSelection::* )(  ) const)( &::SireMol::AtomSelection::selectedAtoms ) )    
        .def( 
            "selectedAtoms"
            , (::QSet<SireMol::CGAtomID> ( ::SireMol::AtomSelection::* )( ::SireMol::CutGroupID ) const)( &::SireMol::AtomSelection::selectedAtoms )
            , ( bp::arg("cgid") ) )    
        .def( 
            "selectedAtoms"
            , (::QSet<SireMol::ResNumAtomID> ( ::SireMol::AtomSelection::* )( ::SireMol::ResNum ) const)( &::SireMol::AtomSelection::selectedAtoms )
            , ( bp::arg("resnum") ) )    
        .def( 
            "selectedCutGroups"
            , &::SireMol::AtomSelection::selectedCutGroups )    
        .def( 
            "selectedNone"
            , (bool ( ::SireMol::AtomSelection::* )(  ) const)( &::SireMol::AtomSelection::selectedNone ) )    
        .def( 
            "selectedNone"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::CutGroupID ) const)( &::SireMol::AtomSelection::selectedNone )
            , ( bp::arg("cgid") ) )    
        .def( 
            "selectedNone"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::ResNum ) const)( &::SireMol::AtomSelection::selectedNone )
            , ( bp::arg("resnum") ) )    
        .def( 
            "selectedNone"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::selectedNone )
            , ( bp::arg("selection") ) )    
        .def( 
            "selectedNone"
            , (bool ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::selectedNone )
            , ( bp::arg("selection") ) )    
        .def( 
            "selectedResidues"
            , &::SireMol::AtomSelection::selectedResidues )    
        .def( 
            "setSelection"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::CutGroupID ) const)( &::SireMol::AtomSelection::setSelection )
            , ( bp::arg("cgid") ) )    
        .def( 
            "setSelection"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::ResNum ) const)( &::SireMol::AtomSelection::setSelection )
            , ( bp::arg("resnum") ) )    
        .def( 
            "setSelection"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::CGAtomID const & ) const)( &::SireMol::AtomSelection::setSelection )
            , ( bp::arg("cgatomid") ) )    
        .def( 
            "setSelection"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::IDMolAtom const & ) const)( &::SireMol::AtomSelection::setSelection )
            , ( bp::arg("atomid") ) )    
        .def( 
            "setSelection"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::setSelection )
            , ( bp::arg("selection") ) )    
        .def( 
            "setSelection"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::setSelection )
            , ( bp::arg("selection") ) )    
        .def( 
            "subtract"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::subtract )
            , ( bp::arg("selection") ) )    
        .def( 
            "subtract"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::subtract )
            , ( bp::arg("selection") ) )    
        .def( 
            "unite"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::AtomSelection const & ) const)( &::SireMol::AtomSelection::unite )
            , ( bp::arg("selection") ) )    
        .def( 
            "unite"
            , (::SireMol::AtomSelection ( ::SireMol::AtomSelection::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::AtomSelection::unite )
            , ( bp::arg("selection") ) )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::AtomSelection >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::AtomSelection >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
