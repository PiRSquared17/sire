// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ResInternalInfo_Angle_.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

void register_ResInternalInfo_Angle__class(){

    bp::class_< SireMM::ResInternalInfo<SireMol::Angle>, bp::bases< SireMM::detail::InternalInfo<SireMol::Angle> > >( "ResInternalInfo_Angle_" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMM::detail::InternalInfo<SireMol::Angle> const & >(( bp::arg("other") )) )    
        .def( bp::init< SireMM::MolInternalInfo<SireMol::Angle> const &, SireMol::ResNum >(( bp::arg("molinfo"), bp::arg("resnum") )) );

}
