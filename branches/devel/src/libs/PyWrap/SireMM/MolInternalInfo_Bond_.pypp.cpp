// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MolInternalInfo_Bond_.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

void register_MolInternalInfo_Bond__class(){

    bp::class_< SireMM::MolInternalInfo<SireMol::Bond>, bp::bases< SireMM::detail::InternalInfo<SireMol::Bond> > >( "MolInternalInfo_Bond_" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") )) )    
        .def( 
            "addInternal"
            , (::SireMol::GroupIndexID ( ::SireMM::MolInternalInfo<SireMol::Bond>::* )( ::SireMol::Bond const & ) )( &::SireMM::MolInternalInfo<SireMol::Bond>::addInternal )
            , ( bp::arg("internal") ) )    
        .def( 
            "removeInterInternals"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) )( &::SireMM::MolInternalInfo<SireMol::Bond>::removeInterInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeInternal"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Bond>::* )( ::SireMol::Bond const & ) )( &::SireMM::MolInternalInfo<SireMol::Bond>::removeInternal )
            , ( bp::arg("internal") ) )    
        .def( 
            "removeInternal"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Bond>::* )( ::SireMol::GroupIndexID const & ) )( &::SireMM::MolInternalInfo<SireMol::Bond>::removeInternal )
            , ( bp::arg("id") ) )    
        .def( 
            "removeInternal"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Bond>::* )( ::SireMol::GroupID ) )( &::SireMM::MolInternalInfo<SireMol::Bond>::removeInternal )
            , ( bp::arg("id") ) )    
        .def( 
            "removeInternals"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Bond>::* )(  ) )( &::SireMM::MolInternalInfo<SireMol::Bond>::removeInternals ) )    
        .def( 
            "removeInternals"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Bond>::* )( ::SireMol::GroupID ) )( &::SireMM::MolInternalInfo<SireMol::Bond>::removeInternals )
            , ( bp::arg("id") ) )    
        .def( 
            "removeInternals"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) )( &::SireMM::MolInternalInfo<SireMol::Bond>::removeInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeIntraInternals"
            , (void ( ::SireMM::MolInternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) )( &::SireMM::MolInternalInfo<SireMol::Bond>::removeIntraInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "residue"
            , (::SireMM::ResInternalInfo<SireMol::Bond> ( ::SireMM::MolInternalInfo<SireMol::Bond>::* )( ::SireMol::ResNum ) const)( &::SireMM::MolInternalInfo<SireMol::Bond>::residue )
            , ( bp::arg("resnum") ) );

}
