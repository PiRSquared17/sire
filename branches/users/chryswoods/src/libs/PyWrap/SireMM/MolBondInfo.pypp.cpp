// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MolBondInfo.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

SireMM::MolBondInfo __copy__(const SireMM::MolBondInfo &other){ return SireMM::MolBondInfo(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::MolBondInfo&){ return "SireMM::MolBondInfo";}

void register_MolBondInfo_class(){

    bp::class_< SireMM::MolBondInfo, bp::bases< SireMM::MolInternalInfo<SireMol::Bond> > >( "MolBondInfo" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") )) )    
        .def( bp::init< SireMol::Molecule const &, SireMM::InternalGenerator<SireMM::MolBondInfo> const & >(( bp::arg("mol"), bp::arg("bondgenerator") )) )    
        .def( 
            "addBond"
            , (::SireMol::GroupIndexID ( ::SireMM::MolBondInfo::* )( ::SireMol::Bond const & ) )( &::SireMM::MolBondInfo::addBond )
            , ( bp::arg("bond") ) )    
        .def( 
            "bonds"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::MolBondInfo::* )(  ) const)( &::SireMM::MolBondInfo::bonds ) )    
        .def( 
            "bonds"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::MolBondInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolBondInfo::bonds )
            , ( bp::arg("resnum") ) )    
        .def( 
            "bonds"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::MolBondInfo::* )( ::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolBondInfo::bonds )
            , ( bp::arg("res0"), bp::arg("res1") ) )    
        .def( 
            "getBond"
            , (::SireMol::Bond ( ::SireMM::MolBondInfo::* )( ::SireMol::GroupIndexID const & ) const)( &::SireMM::MolBondInfo::getBond )
            , ( bp::arg("id") ) )    
        .def( 
            "interBonds"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::MolBondInfo::* )(  ) const)( &::SireMM::MolBondInfo::interBonds ) )    
        .def( 
            "interBonds"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::MolBondInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolBondInfo::interBonds )
            , ( bp::arg("resnum") ) )    
        .def( 
            "intraBonds"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::MolBondInfo::* )(  ) const)( &::SireMM::MolBondInfo::intraBonds ) )    
        .def( 
            "intraBonds"
            , (::SireMM::detail::InternalGroupIterator<SireMol::Bond> ( ::SireMM::MolBondInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolBondInfo::intraBonds )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nBonds"
            , (int ( ::SireMM::MolBondInfo::* )(  ) const)( &::SireMM::MolBondInfo::nBonds ) )    
        .def( 
            "nBonds"
            , (int ( ::SireMM::MolBondInfo::* )( ::SireMol::GroupID ) const)( &::SireMM::MolBondInfo::nBonds )
            , ( bp::arg("group") ) )    
        .def( 
            "nBonds"
            , (int ( ::SireMM::MolBondInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolBondInfo::nBonds )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nBonds"
            , (int ( ::SireMM::MolBondInfo::* )( ::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolBondInfo::nBonds )
            , ( bp::arg("res0"), bp::arg("res1") ) )    
        .def( 
            "nInterBonds"
            , (int ( ::SireMM::MolBondInfo::* )(  ) const)( &::SireMM::MolBondInfo::nInterBonds ) )    
        .def( 
            "nInterBonds"
            , (int ( ::SireMM::MolBondInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolBondInfo::nInterBonds )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nIntraBonds"
            , (int ( ::SireMM::MolBondInfo::* )(  ) const)( &::SireMM::MolBondInfo::nIntraBonds ) )    
        .def( 
            "nIntraBonds"
            , (int ( ::SireMM::MolBondInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolBondInfo::nIntraBonds )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeBond"
            , (void ( ::SireMM::MolBondInfo::* )( ::SireMol::Bond const & ) )( &::SireMM::MolBondInfo::removeBond )
            , ( bp::arg("bond") ) )    
        .def( 
            "removeBond"
            , (void ( ::SireMM::MolBondInfo::* )( ::SireMol::GroupIndexID const & ) )( &::SireMM::MolBondInfo::removeBond )
            , ( bp::arg("id") ) )    
        .def( 
            "residue"
            , (::SireMM::ResBondInfo ( ::SireMM::MolBondInfo::* )( ::SireMol::ResNum ) const)( &::SireMM::MolBondInfo::residue )
            , ( bp::arg("resnum") ) )    
        .def( 
            "residuesBonded"
            , (bool ( ::SireMM::MolBondInfo::* )( ::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMM::MolBondInfo::residuesBonded )
            , ( bp::arg("res0"), bp::arg("res1") ) )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMM::MolBondInfo >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMM::MolBondInfo >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
