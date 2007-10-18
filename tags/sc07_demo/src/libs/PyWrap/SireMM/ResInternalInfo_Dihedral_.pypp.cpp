// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ResInternalInfo_Dihedral_.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

SireMM::ResInternalInfo<SireMol::Dihedral> __copy__(const SireMM::ResInternalInfo<SireMol::Dihedral> &other){ return SireMM::ResInternalInfo<SireMol::Dihedral>(other); }

const char* pvt_get_name(const SireMM::ResInternalInfo<SireMol::Dihedral>&){ return "SireMM::ResInternalInfo<SireMol::Dihedral>";}

void register_ResInternalInfo_Dihedral__class(){

    bp::class_< SireMM::ResInternalInfo<SireMol::Dihedral>, bp::bases< SireMM::detail::InternalInfo<SireMol::Dihedral> > >( "ResInternalInfo_Dihedral_" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMM::MolInternalInfo<SireMol::Dihedral> const &, SireMol::ResNum >(( bp::arg("molinfo"), bp::arg("resnum") )) )    
        .def( bp::init< SireMM::detail::InternalInfo<SireMol::Dihedral> const & >(( bp::arg("other") )) )    
        .def( "__copy__", &__copy__)    
        .def( "__str__", &pvt_get_name);

}
