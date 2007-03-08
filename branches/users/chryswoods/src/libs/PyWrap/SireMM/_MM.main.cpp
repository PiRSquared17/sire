// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License


#include "boost/python.hpp"

#include "siremm_headers.h"

#include "SireMol/molecule.h"

#include "SireMol/residue.h"

#include "SireMol/newatom.h"

#include "SireMol/atom.h"

#include "AngleAtomID.pypp.hpp"

#include "AngleDB.pypp.hpp"

#include "AngleGenerator.pypp.hpp"

#include "AngleResID.pypp.hpp"

#include "AngleTable.pypp.hpp"

#include "AngleTable_Expression_.pypp.hpp"

#include "ArithmeticCombiningRules.pypp.hpp"

#include "AtomTable_ChargeParameter_.pypp.hpp"

#include "AtomTable_LJParameter_.pypp.hpp"

#include "AtomicCharges.pypp.hpp"

#include "AtomicLJs.pypp.hpp"

#include "BondAtomID.pypp.hpp"

#include "BondDB.pypp.hpp"

#include "BondGenerator.pypp.hpp"

#include "BondResID.pypp.hpp"

#include "BondTable.pypp.hpp"

#include "BondTable_Expression_.pypp.hpp"

#include "CLJFF.pypp.hpp"

#include "CLJPair.pypp.hpp"

#include "CLJParameter.pypp.hpp"

#include "ChargeDB.pypp.hpp"

#include "ChargeParameter.pypp.hpp"

#include "ChargeTable.pypp.hpp"

#include "CombiningRuleBase.pypp.hpp"

#include "CombiningRules.pypp.hpp"

#include "DihedralAtomID.pypp.hpp"

#include "DihedralDB.pypp.hpp"

#include "DihedralGenerator.pypp.hpp"

#include "DihedralResID.pypp.hpp"

#include "DihedralTable.pypp.hpp"

#include "DihedralTable_Expression_.pypp.hpp"

#include "GeometricCombiningRules.pypp.hpp"

#include "HarmonicSwitchingFunction.pypp.hpp"

#include "InterCLJFF.pypp.hpp"

#include "InterLJFF.pypp.hpp"

#include "InternalDB_MolAngleInfo_.pypp.hpp"

#include "InternalDB_MolBondInfo_.pypp.hpp"

#include "InternalDB_MolDihedralInfo_.pypp.hpp"

#include "InternalGenerator_MolAngleInfo_.pypp.hpp"

#include "InternalGenerator_MolBondInfo_.pypp.hpp"

#include "InternalGenerator_MolDihedralInfo_.pypp.hpp"

#include "InternalInfo_Angle_.pypp.hpp"

#include "InternalInfo_Bond_.pypp.hpp"

#include "InternalInfo_Dihedral_.pypp.hpp"

#include "InternalTableBase_MolAngleInfo_.pypp.hpp"

#include "InternalTableBase_MolBondInfo_.pypp.hpp"

#include "InternalTableBase_MolDihedralInfo_.pypp.hpp"

#include "InternalTable_MolAngleInfo_Expression_.pypp.hpp"

#include "InternalTable_MolBondInfo_Expression_.pypp.hpp"

#include "InternalTable_MolDihedralInfo_Expression_.pypp.hpp"

#include "LJDB.pypp.hpp"

#include "LJFF.pypp.hpp"

#include "LJPair.pypp.hpp"

#include "LJParameter.pypp.hpp"

#include "LJTable.pypp.hpp"

#include "MolAngleInfo.pypp.hpp"

#include "MolBondInfo.pypp.hpp"

#include "MolDihedralInfo.pypp.hpp"

#include "MolInternalInfo_Angle_.pypp.hpp"

#include "MolInternalInfo_Bond_.pypp.hpp"

#include "MolInternalInfo_Dihedral_.pypp.hpp"

#include "NoCutoff.pypp.hpp"

#include "ResAngleInfo.pypp.hpp"

#include "ResBondInfo.pypp.hpp"

#include "ResDihedralInfo.pypp.hpp"

#include "ResInternalInfo_Angle_.pypp.hpp"

#include "ResInternalInfo_Bond_.pypp.hpp"

#include "ResInternalInfo_Dihedral_.pypp.hpp"

#include "SwitchFuncBase.pypp.hpp"

#include "SwitchingFunction.pypp.hpp"

#include "Tip4PFF.pypp.hpp"

#include "UreyBradleyDB.pypp.hpp"

#include "UsePassedAngles.pypp.hpp"

#include "UsePassedBonds.pypp.hpp"

#include "UsePassedDihedrals.pypp.hpp"

#include "UsePassedInternals_MolAngleInfo_.pypp.hpp"

#include "UsePassedInternals_MolBondInfo_.pypp.hpp"

#include "UsePassedInternals_MolDihedralInfo_.pypp.hpp"

#include "assign_angles.pypp.hpp"

#include "assign_bonds.pypp.hpp"

#include "assign_dihedrals.pypp.hpp"

#include "assign_internals_MolAngleInfo_.pypp.hpp"

#include "assign_internals_MolBondInfo_.pypp.hpp"

#include "assign_internals_MolDihedralInfo_.pypp.hpp"

#include "siremm_containers.h"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(_MM){
    register_SireMM_containers();

    register_AtomTable_ChargeParameter__class();

    register_AtomTable_LJParameter__class();

    register_AngleAtomID_class();

    register_InternalDB_MolAngleInfo__class();

    register_AngleDB_class();

    register_InternalGenerator_MolAngleInfo__class();

    register_AngleGenerator_class();

    register_AngleResID_class();

    register_InternalTableBase_MolAngleInfo__class();

    register_InternalTable_MolAngleInfo_Expression__class();

    register_AngleTable_Expression__class();

    register_AngleTable_class();

    register_CombiningRuleBase_class();

    register_ArithmeticCombiningRules_class();

    register_AtomicCharges_class();

    register_AtomicLJs_class();

    register_BondAtomID_class();

    register_InternalDB_MolBondInfo__class();

    register_BondDB_class();

    register_InternalGenerator_MolBondInfo__class();

    register_BondGenerator_class();

    register_BondResID_class();

    register_InternalTableBase_MolBondInfo__class();

    register_InternalTable_MolBondInfo_Expression__class();

    register_BondTable_Expression__class();

    register_BondTable_class();

    register_CLJFF_class();

    register_CLJPair_class();

    register_ChargeParameter_class();

    register_LJParameter_class();

    register_CLJParameter_class();

    register_ChargeDB_class();

    register_ChargeTable_class();

    register_CombiningRules_class();

    register_DihedralAtomID_class();

    register_InternalDB_MolDihedralInfo__class();

    register_DihedralDB_class();

    register_InternalGenerator_MolDihedralInfo__class();

    register_DihedralGenerator_class();

    register_DihedralResID_class();

    register_InternalTableBase_MolDihedralInfo__class();

    register_InternalTable_MolDihedralInfo_Expression__class();

    register_DihedralTable_Expression__class();

    register_DihedralTable_class();

    register_GeometricCombiningRules_class();

    register_SwitchFuncBase_class();

    register_HarmonicSwitchingFunction_class();

    register_InterCLJFF_class();

    register_LJFF_class();

    register_InterLJFF_class();

    register_LJDB_class();

    register_LJPair_class();

    register_LJTable_class();

    register_InternalInfo_Angle__class();

    register_MolInternalInfo_Angle__class();

    register_MolAngleInfo_class();

    register_InternalInfo_Bond__class();

    register_MolInternalInfo_Bond__class();

    register_MolBondInfo_class();

    register_InternalInfo_Dihedral__class();

    register_MolInternalInfo_Dihedral__class();

    register_MolDihedralInfo_class();

    register_NoCutoff_class();

    register_ResInternalInfo_Angle__class();

    register_ResAngleInfo_class();

    register_ResInternalInfo_Bond__class();

    register_ResBondInfo_class();

    register_ResInternalInfo_Dihedral__class();

    register_ResDihedralInfo_class();

    register_SwitchingFunction_class();

    register_Tip4PFF_class();

    register_UreyBradleyDB_class();

    register_UsePassedInternals_MolAngleInfo__class();

    register_UsePassedAngles_class();

    register_UsePassedInternals_MolBondInfo__class();

    register_UsePassedBonds_class();

    register_UsePassedInternals_MolDihedralInfo__class();

    register_UsePassedDihedrals_class();

    register_assign_internals_MolAngleInfo__class();

    register_assign_angles_class();

    register_assign_internals_MolBondInfo__class();

    register_assign_bonds_class();

    register_assign_internals_MolDihedralInfo__class();

    bp::implicitly_convertible< QVector< QVector<SireMM::ChargeParameter> >, SireMM::AtomicCharges >();

    bp::implicitly_convertible< QVector< SireMM::ChargeParameter >, SireMM::AtomicCharges >();

    bp::implicitly_convertible< QVector< QVector<SireMM::LJParameter> >, SireMM::AtomicLJs >();

    bp::implicitly_convertible< QVector< SireMM::LJParameter >, SireMM::AtomicLJs >();

    bp::implicitly_convertible< double, SireMM::ChargeParameter >();

    bp::implicitly_convertible< const SireMM::CombiningRuleBase&, SireMM::CombiningRules >();

    bp::implicitly_convertible< SireMM::SwitchFuncBase, SireMM::SwitchingFunction >();

    bp::implicitly_convertible< QSet<SireMol::Angle>, SireMM::UsePassedAngles >();

    bp::implicitly_convertible< const SireMM::AngleGeneratorBase&, SireMM::assign_angles >();

    bp::implicitly_convertible< SireMM::UsePassedAngles, SireMM::assign_angles >();

    bp::implicitly_convertible< QSet<SireMol::Bond>, SireMM::UsePassedBonds >();

    bp::implicitly_convertible< const SireMM::BondGeneratorBase&, SireMM::assign_bonds >();

    bp::implicitly_convertible< SireMM::UsePassedBonds, SireMM::assign_bonds >();

    bp::implicitly_convertible< QSet<SireMol::Dihedral>, SireMM::UsePassedDihedrals >();

    bp::implicitly_convertible< const SireMM::DihedralGeneratorBase&, SireMM::assign_dihedrals >();

    bp::implicitly_convertible< SireMM::UsePassedDihedrals, SireMM::assign_dihedrals >();

    register_assign_dihedrals_class();
}

