// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MolInternalInfo_Angle_.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

SireMM::MolInternalInfo<SireMol::Angle> __copy__(const SireMM::MolInternalInfo<SireMol::Angle> &other){ return SireMM::MolInternalInfo<SireMol::Angle>(other); }

const char* pvt_get_name(const SireMM::MolInternalInfo<SireMol::Angle>&){ return "SireMM::MolInternalInfo<SireMol::Angle>";}

void register_MolInternalInfo_Angle__class(){

    bp::class_< SireMM::MolInternalInfo<SireMol::Angle>, bp::bases< SireMM::detail::InternalInfo<SireMol::Angle> > >( "MolInternalInfo_Angle_" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") )) )    
        .def( 
            "addInternal"
            , (::SireMol::GroupIndexID ( ::SireMM::MolInternalInfo<SireMol::Angle>::* )( ::SireMol::Angle const & ) )( &::SireMM::MolInternalInfo<SireMol::Angle>::addInternal )
            , ( bp::arg("internal") ) )    
        .def( 
            "removeInterInternals"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Angle>::* )( ::SireMol::ResNum ) )( &::SireMM::MolInternalInfo<SireMol::Angle>::removeInterInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeInternal"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Angle>::* )( ::SireMol::Angle const & ) )( &::SireMM::MolInternalInfo<SireMol::Angle>::removeInternal )
            , ( bp::arg("internal") ) )    
        .def( 
            "removeInternal"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Angle>::* )( ::SireMol::GroupIndexID const & ) )( &::SireMM::MolInternalInfo<SireMol::Angle>::removeInternal )
            , ( bp::arg("id") ) )    
        .def( 
            "removeInternal"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Angle>::* )( ::SireMol::GroupID ) )( &::SireMM::MolInternalInfo<SireMol::Angle>::removeInternal )
            , ( bp::arg("id") ) )    
        .def( 
            "removeInternals"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Angle>::* )(  ) )( &::SireMM::MolInternalInfo<SireMol::Angle>::removeInternals ) )    
        .def( 
            "removeInternals"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Angle>::* )( ::SireMol::GroupID ) )( &::SireMM::MolInternalInfo<SireMol::Angle>::removeInternals )
            , ( bp::arg("id") ) )    
        .def( 
            "removeInternals"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Angle>::* )( ::SireMol::ResNum ) )( &::SireMM::MolInternalInfo<SireMol::Angle>::removeInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeIntraInternals"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Angle>::* )( ::SireMol::ResNum ) )( &::SireMM::MolInternalInfo<SireMol::Angle>::removeIntraInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "residue"
            , (::SireMM::ResInternalInfo<SireMol::Angle> ( ::SireMM::MolInternalInfo<SireMol::Angle>::* )( ::SireMol::ResNum ) const)( &::SireMM::MolInternalInfo<SireMol::Angle>::residue )
            , ( bp::arg("resnum") ) )    
        .def( "__copy__", &__copy__)    
        .def( "__str__", &pvt_get_name);

}
