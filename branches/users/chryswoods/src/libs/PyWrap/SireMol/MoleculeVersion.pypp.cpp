// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MoleculeVersion.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMol/moleculedata.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"
#include "SireVol/space.h"

namespace bp = boost::python;

SireMol::MoleculeVersion __copy__(const SireMol::MoleculeVersion &other){ return SireMol::MoleculeVersion(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_MoleculeVersion_class(){

    bp::class_< SireMol::MoleculeVersion >( "MoleculeVersion" )    
        .def( bp::init< >() )    
        .def( bp::init< quint32, quint32 >(( bp::arg("major"), bp::arg("minor") )) )    
        .def( 
            "incrementMajor"
            , &::SireMol::MoleculeVersion::incrementMajor )    
        .def( 
            "incrementMinor"
            , &::SireMol::MoleculeVersion::incrementMinor )    
        .def( 
            "major"
            , &::SireMol::MoleculeVersion::major )    
        .def( 
            "majorDifference"
            , &::SireMol::MoleculeVersion::majorDifference
            , ( bp::arg("other") ) )    
        .def( 
            "minor"
            , &::SireMol::MoleculeVersion::minor )    
        .def( 
            "minorDifference"
            , &::SireMol::MoleculeVersion::minorDifference
            , ( bp::arg("other") ) )    
        .def( bp::self != bp::self )    
        .def( bp::self < bp::self )    
        .def( bp::self <= bp::self )    
        .def( bp::self == bp::self )    
        .def( bp::self > bp::self )    
        .def( bp::self >= bp::self )    
        .def( 
            "sameVersion"
            , &::SireMol::MoleculeVersion::sameVersion
            , ( bp::arg("other") ) )    
        .def( 
            "toString"
            , &::SireMol::MoleculeVersion::toString )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::MoleculeVersion >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::MoleculeVersion >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireMol::MoleculeVersion > );

}
