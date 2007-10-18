// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "InternalGenerator_MolBondInfo_.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

const char* pvt_get_name(const SireMM::InternalGenerator<SireMM::MolBondInfo>&){ return "SireMM::InternalGenerator<SireMM::MolBondInfo>";}

void register_InternalGenerator_MolBondInfo__class(){

    bp::class_< SireMM::InternalGenerator<SireMM::MolBondInfo>, boost::noncopyable >( "InternalGenerator_MolBondInfo_", bp::no_init )    
        .def( 
            "generate"
            , (void ( ::SireMM::InternalGenerator<SireMM::MolBondInfo>::* )( ::SireMol::Molecule const &,::SireMM::MolBondInfo & ) const)( &::SireMM::InternalGenerator<SireMM::MolBondInfo>::generate )
            , ( bp::arg("molecule"), bp::arg("internalinfo") ) )    
        .def( 
            "generate"
            , (::SireMM::MolBondInfo ( ::SireMM::InternalGenerator<SireMM::MolBondInfo>::* )( ::SireMol::Molecule const & ) const)( &::SireMM::InternalGenerator<SireMM::MolBondInfo>::generate )
            , ( bp::arg("molecule") ) )    
        .def( 
            "what"
            , (char const * ( ::SireMM::InternalGenerator<SireMM::MolBondInfo>::* )(  ) const)( &::SireMM::InternalGenerator<SireMM::MolBondInfo>::what ) )    
        .def( "__str__", &pvt_get_name);

}
