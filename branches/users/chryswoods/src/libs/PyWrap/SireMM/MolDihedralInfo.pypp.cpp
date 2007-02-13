// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MolDihedralInfo.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

void register_MolDihedralInfo_class(){

    bp::class_< SireMM::MolDihedralInfo, bp::bases< SireMM::MolInternalInfo<SireMol::Dihedral> > >( "MolDihedralInfo" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") )) )    
        .def( bp::init< SireMol::Molecule const &, SireMM::InternalGenerator<SireMM::MolDihedralInfo> const & >(( bp::arg("mol"), bp::arg("dihedralgenerator") )) )    
        .def( 
            "addDihedral"
            , (::SireMol::GroupIndexID ( ::SireMM::MolDihedralInfo::* )( ::SireMol::Dihedral const & ) )( &::SireMM::MolDihedralInfo::addDihedral )
            , ( bp::arg("dihedral") ) )    
        .def( 
            "dihedrals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Dihedral> ( ::SireMM::MolDihedralInfo::* )(  ) const)( &::SireMM::MolDihedralInfo::dihedrals ) )    
        .def( 
            "dihedrals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Dihedral> ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::dihedrals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "dihedrals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Dihedral> ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::dihedrals )
            , ( bp::arg("res0"), bp::arg("res1") ) )    
        .def( 
            "dihedrals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Dihedral> ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum,::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::dihedrals )
            , ( bp::arg("res0"), bp::arg("res1"), bp::arg("res2") ) )    
        .def( 
            "dihedrals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Dihedral> ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum,::SireMol::ResNum,::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::dihedrals )
            , ( bp::arg("res0"), bp::arg("res1"), bp::arg("res2"), bp::arg("res3") ) )    
        .def( 
            "interDihedrals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Dihedral> ( ::SireMM::MolDihedralInfo::* )(  ) const)( &::SireMM::MolDihedralInfo::interDihedrals ) )    
        .def( 
            "interDihedrals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Dihedral> ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::interDihedrals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "intraDihedrals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Dihedral> ( ::SireMM::MolDihedralInfo::* )(  ) const)( &::SireMM::MolDihedralInfo::intraDihedrals ) )    
        .def( 
            "intraDihedrals"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Dihedral> ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::intraDihedrals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nDihedrals"
            , (int ( ::SireMM::MolDihedralInfo::* )(  ) const)( &::SireMM::MolDihedralInfo::nDihedrals ) )    
        .def( 
            "nDihedrals"
            , (int ( ::SireMM::MolDihedralInfo::* )( ::SireMol::GroupID ) const)( &::SireMM::MolDihedralInfo::nDihedrals )
            , ( bp::arg("group") ) )    
        .def( 
            "nDihedrals"
            , (int ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::nDihedrals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nDihedrals"
            , (int ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::nDihedrals )
            , ( bp::arg("res0"), bp::arg("res1") ) )    
        .def( 
            "nDihedrals"
            , (int ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum,::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::nDihedrals )
            , ( bp::arg("res0"), bp::arg("res1"), bp::arg("res2") ) )    
        .def( 
            "nDihedrals"
            , (int ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum,::SireMol::ResNum,::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::nDihedrals )
            , ( bp::arg("res0"), bp::arg("res1"), bp::arg("res2"), bp::arg("res3") ) )    
        .def( 
            "nInterDihedrals"
            , (int ( ::SireMM::MolDihedralInfo::* )(  ) const)( &::SireMM::MolDihedralInfo::nInterDihedrals ) )    
        .def( 
            "nInterDihedrals"
            , (int ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::nInterDihedrals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nIntraDihedrals"
            , (int ( ::SireMM::MolDihedralInfo::* )(  ) const)( &::SireMM::MolDihedralInfo::nIntraDihedrals ) )    
        .def( 
            "nIntraDihedrals"
            , (int ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::nIntraDihedrals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeDihedral"
            , (void ( ::SireMM::MolDihedralInfo::* )( ::SireMol::Dihedral const & ) )( &::SireMM::MolDihedralInfo::removeDihedral )
            , ( bp::arg("dihedral") ) )    
        .def( 
            "removeDihedral"
            , (void ( ::SireMM::MolDihedralInfo::* )( ::SireMol::GroupIndexID const & ) )( &::SireMM::MolDihedralInfo::removeDihedral )
            , ( bp::arg("id") ) )    
        .def( 
            "residue"
            , (::SireMM::ResDihedralInfo ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::residue )
            , ( bp::arg("resnum") ) )    
        .def( 
            "residuesDihedraled"
            , (bool ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::residuesDihedraled )
            , ( bp::arg("res0"), bp::arg("res1") ) )    
        .def( 
            "residuesDihedraled"
            , (bool ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum,::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::residuesDihedraled )
            , ( bp::arg("res0"), bp::arg("res1"), bp::arg("res2") ) )    
        .def( 
            "residuesDihedraled"
            , (bool ( ::SireMM::MolDihedralInfo::* )( ::SireMol::ResNum,::SireMol::ResNum,::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolDihedralInfo::residuesDihedraled )
            , ( bp::arg("res0"), bp::arg("res1"), bp::arg("res2"), bp::arg("res3") ) );

}
