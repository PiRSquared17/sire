// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "AtomInfo.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_AtomInfo_class(){

    bp::class_< SireMol::AtomInfo, bp::bases< SireMol::AtomIndex, SireMol::Element > >( "AtomInfo" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::AtomIndex const &, SireMol::Element const & >(( bp::arg("atomindex"), bp::arg("element") )) )    
        .def( bp::init< SireMol::AtomInfo const &, SireMol::ResNum >(( bp::arg("other"), bp::arg("resnum") )) )    
        .def( bp::init< SireMol::AtomNum, SireMol::AtomIndex const &, SireMol::Element const & >(( bp::arg("atmnum"), bp::arg("atomindex"), bp::arg("element") )) )    
        .def( bp::init< SireMol::AtomNum, SireMol::AtomInfo const & >(( bp::arg("atmnum"), bp::arg("other") )) )    
        .def( bp::init< SireMol::AtomNum, SireMol::AtomInfo const &, SireMol::ResNum >(( bp::arg("atmnum"), bp::arg("other"), bp::arg("resnum") )) )    
        .def( bp::init< QString const & >(( bp::arg("name") )) )    
        .def( bp::init< QString const &, QString const & >(( bp::arg("name"), bp::arg("element") )) )    
        .def( bp::init< QString const &, SireMol::ResNum >(( bp::arg("name"), bp::arg("resnum") )) )    
        .def( bp::init< QString const &, SireMol::ResNum, QString const & >(( bp::arg("name"), bp::arg("resnum"), bp::arg("element") )) )    
        .def( bp::init< SireMol::AtomIndex const & >(( bp::arg("atm") )) )    
        .def( bp::init< QString const &, SireMol::AtomInfo const & >(( bp::arg("name"), bp::arg("other") )) )    
        .def( bp::init< SireMol::AtomNum, QString const & >(( bp::arg("num"), bp::arg("name") )) )    
        .def( bp::init< SireMol::AtomNum, QString const &, QString const & >(( bp::arg("num"), bp::arg("name"), bp::arg("element") )) )    
        .def( 
            "atomNum"
            , &::SireMol::AtomInfo::atomNum )    
        .def( 
            "name"
            , &::SireMol::AtomInfo::name )    
        .def( 
            "number"
            , &::SireMol::AtomInfo::number )    
        .def( bp::self != bp::other< SireMol::AtomIndex >() )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::other< SireMol::AtomIndex >() )    
        .def( bp::self == bp::self )    
        .def( 
            "toString"
            , &::SireMol::AtomInfo::toString );

}
