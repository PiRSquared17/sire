// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Version.pypp.hpp"
#include "boost/python.hpp"
#include "sirebase_headers.h"

namespace bp = boost::python;

void register_Version_class(){

    bp::class_< SireBase::Version >( "Version", bp::init< bp::optional< quint32, quint32 > >(( bp::arg("major")=(unsigned int)(0), bp::arg("minor")=(unsigned int)(0) )) )    
        .def( 
            "major"
            , &::SireBase::Version::major )    
        .def( 
            "minor"
            , &::SireBase::Version::minor )    
        .def( bp::self != bp::self )    
        .def( bp::self < bp::self )    
        .def( bp::self <= bp::self )    
        .def( bp::self == bp::self )    
        .def( bp::self > bp::self )    
        .def( bp::self >= bp::self )    
        .def( 
            "sameMajorVersion"
            , &::SireBase::Version::sameMajorVersion
            , ( bp::arg("other") ) )    
        .def( 
            "sameMinorVersion"
            , &::SireBase::Version::sameMinorVersion
            , ( bp::arg("other") ) )    
        .def( 
            "setMajor"
            , &::SireBase::Version::setMajor
            , ( bp::arg("major") ) )    
        .def( 
            "setMinor"
            , &::SireBase::Version::setMinor
            , ( bp::arg("minor") ) )    
        .def( 
            "toString"
            , &::SireBase::Version::toString );

}
