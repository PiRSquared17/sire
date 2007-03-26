// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License


#include "boost/python.hpp"

#include "siremol_headers.h"

#include "SireMaths/angle.h"

#include "SireMaths/quaternion.h"

#include "SireMaths/matrix.h"

#include "SireMaths/triangle.h"

#include "SireMaths/line.h"

#include "SireMaths/torsion.h"

#include "Angle.pypp.hpp"

#include "Atom.pypp.hpp"

#include "AtomCutting.pypp.hpp"

#include "AtomID.pypp.hpp"

#include "AtomIDGroup.pypp.hpp"

#include "AtomIndex.pypp.hpp"

#include "AtomInfo.pypp.hpp"

#include "AtomInfoGroup.pypp.hpp"

#include "AtomNum.pypp.hpp"

#include "AtomSelection.pypp.hpp"

#include "AtomicProperties.pypp.hpp"

#include "Bond.pypp.hpp"

#include "BondAddingFunction.pypp.hpp"

#include "CGAtomID.pypp.hpp"

#include "CGNumAtomID.pypp.hpp"

#include "CutGroup.pypp.hpp"

#include "CutGroupID.pypp.hpp"

#include "CutGroupNum.pypp.hpp"

#include "CuttingFunction.pypp.hpp"

#include "CuttingFunctionBase.pypp.hpp"

#include "Dihedral.pypp.hpp"

#include "EditMol.pypp.hpp"

#include "EditRes.pypp.hpp"

#include "Element.pypp.hpp"

#include "GroupID.pypp.hpp"

#include "GroupIndexID.pypp.hpp"

#include "IDBase.pypp.hpp"

#include "IDMolAtom.pypp.hpp"

#include "Improper.pypp.hpp"

#include "Index.pypp.hpp"

#include "MolCGAtomID.pypp.hpp"

#include "MolCutGroupID.pypp.hpp"

#include "Molecule.pypp.hpp"

#include "MoleculeBonds.pypp.hpp"

#include "MoleculeCutting.pypp.hpp"

#include "MoleculeGroup.pypp.hpp"

#include "MoleculeGroupID.pypp.hpp"

#include "MoleculeGroups.pypp.hpp"

#include "MoleculeID.pypp.hpp"

#include "MoleculeInfo.pypp.hpp"

#include "MoleculeProperty.pypp.hpp"

#include "MoleculeVersion.pypp.hpp"

#include "NewAtom.pypp.hpp"

#include "ResID.pypp.hpp"

#include "ResIDAtomID.pypp.hpp"

#include "ResNum.pypp.hpp"

#include "ResNumAtomID.pypp.hpp"

#include "ResNumIndexID.pypp.hpp"

#include "Residue.pypp.hpp"

#include "ResidueBonds.pypp.hpp"

#include "ResidueCutting.pypp.hpp"

#include "ResidueInfo.pypp.hpp"

#include "VariantProperty.pypp.hpp"

#include "WeightFunction.pypp.hpp"

#include <boost/tuple/tuple.hpp>

#include "siremol_containers.h"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(_Mol){
    register_SireMol_containers();

    register_VariantProperty_class();

    register_Angle_class();

    register_IDBase_class();

    bp::implicitly_convertible< SireMol::IDBase, quint32 >();

    register_ResNum_class();

    register_AtomIndex_class();

    register_Element_class();

    register_AtomInfo_class();

    register_Atom_class();

    register_CuttingFunctionBase_class();

    register_AtomCutting_class();

    register_AtomID_class();

    register_AtomIDGroup_class();

    register_AtomInfoGroup_class();

    register_AtomNum_class();

    register_AtomSelection_class();

    register_MoleculeProperty_class();

    register_AtomicProperties_class();

    register_Bond_class();

    register_BondAddingFunction_class();

    register_CutGroupID_class();

    register_CGAtomID_class();

    register_CutGroupNum_class();

    register_CGNumAtomID_class();

    register_CutGroup_class();

    register_CuttingFunction_class();

    register_Dihedral_class();

    register_EditMol_class();

    register_EditRes_class();

    register_GroupID_class();

    register_Index_class();

    register_GroupIndexID_class();

    register_IDMolAtom_class();

    register_Improper_class();

    register_MoleculeID_class();

    register_MolCGAtomID_class();

    register_MolCutGroupID_class();

    register_Molecule_class();

    register_MoleculeBonds_class();

    register_MoleculeCutting_class();

    register_MoleculeGroup_class();

    register_MoleculeGroupID_class();

    register_MoleculeGroups_class();

    register_MoleculeInfo_class();

    register_MoleculeVersion_class();

    register_NewAtom_class();

    bp::implicitly_convertible< SireMol::NewAtom, SireMaths::Vector >();

    bp::implicitly_convertible< SireMol::NewAtom, SireMol::Element >();

    register_ResID_class();

    register_ResIDAtomID_class();

    register_ResNumAtomID_class();

    register_ResNumIndexID_class();

    register_Residue_class();

    register_ResidueBonds_class();

    register_ResidueCutting_class();

    register_ResidueInfo_class();

    bp::implicitly_convertible< boost::tuples::tuple<SireMol::AtomIndex,SireMol::AtomIndex,SireMol::AtomIndex>, SireMol::Angle >();

    bp::implicitly_convertible< boost::tuples::tuple<SireMol::AtomIndex,SireMol::AtomIndex,SireMol::AtomIndex,SireMol::AtomIndex>, SireMol::Dihedral >();

    bp::implicitly_convertible< boost::tuples::tuple<SireMol::AtomIndex,SireMol::AtomIndex,SireMol::AtomIndex,SireMol::AtomIndex>, SireMol::Improper >();

    bp::implicitly_convertible< boost::tuples::tuple<SireMol::AtomIndex,SireMol::AtomIndex>, SireMol::Bond >();

    bp::implicitly_convertible< SireMol::IDBase, quint32 >();

    bp::implicitly_convertible< boost::tuples::tuple<QString,SireMol::ResNum>, SireMol::AtomIndex >();

    bp::implicitly_convertible< QString, SireMol::Element >();

    bp::implicitly_convertible< quint32, SireMol::Element >();

    bp::implicitly_convertible< boost::tuples::tuple<SireMol::CutGroupID,SireMol::AtomID>, SireMol::CGAtomID >();

    bp::implicitly_convertible< boost::tuples::tuple<SireMol::ResNum,SireMol::AtomID>, SireMol::ResNumAtomID >();

    bp::implicitly_convertible< boost::tuples::tuple<SireMol::ResID,SireMol::AtomID>, SireMol::ResIDAtomID >();

    bp::implicitly_convertible< boost::tuples::tuple<SireMol::GroupID,SireMol::Index>, SireMol::GroupIndexID >();

    bp::implicitly_convertible< QVector<SireMol::Atom>, SireMol::CutGroup >();

    bp::implicitly_convertible< QList<SireMol::Atom>, SireMol::CutGroup >();

    bp::implicitly_convertible< const SireMol::CuttingFunctionBase&, SireMol::CuttingFunction >();

    bp::implicitly_convertible< SireMol::CGAtomID, SireMol::IDMolAtom >();

    bp::implicitly_convertible< SireMol::CGNumAtomID, SireMol::IDMolAtom >();

    bp::implicitly_convertible< SireMol::ResNumAtomID, SireMol::IDMolAtom >();

    bp::implicitly_convertible< SireMol::ResIDAtomID, SireMol::IDMolAtom >();

    bp::implicitly_convertible< SireMol::AtomIndex, SireMol::IDMolAtom >();

    bp::implicitly_convertible< SireMol::AtomID, SireMol::IDMolAtom >();

    bp::implicitly_convertible< SireMol::MoleculeGroup, SireMol::MoleculeGroups >();

    bp::implicitly_convertible< SireMol::NewAtom, SireMaths::Vector >();

    bp::implicitly_convertible< SireMol::NewAtom, SireMol::Element >();

    bp::implicitly_convertible< const SireMol::PropertyBase&, SireMol::Property >();

    bp::implicitly_convertible< QVariant, SireMol::VariantProperty >();

    register_WeightFunction_class();
}

