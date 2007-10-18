// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Filters.pypp.hpp"
#include "boost/python.hpp"
#include "sireqt_headers.h"

namespace bp = boost::python;

const char* pvt_get_name(const QFlags<QDir::Filter>&){ return "QFlags<QDir::Filter>";}

void register_Filters_class(){

    bp::class_< QFlags<QDir::Filter> >( "Filters", bp::init< QDir::Filter >(( bp::arg("f") )) )    
        .def( bp::init< bp::optional< void * * > >(( bp::arg("arg0")=bp::object() )) )    
        .def( bp::init< QFlag >(( bp::arg("f") )) )    
        .def( "__int__", &QFlags<QDir::Filter>::operator int  )    
        .def( !bp::self )    
        .def( bp::self & bp::other< int >() )    
        .def( bp::self & bp::other< uint >() )    
        .def( bp::self & bp::other< QDir::Filter >() )    
        .def( bp::self &= bp::other< int >() )    
        .def( bp::self &= bp::other< uint >() )    
        .def( bp::self ^ bp::self )    
        .def( bp::self ^ bp::other< QDir::Filter >() )    
        .def( bp::self ^= bp::self )    
        .def( bp::self ^= bp::other< QDir::Filter >() )    
        .def( bp::self | bp::self )    
        .def( bp::self | bp::other< QDir::Filter >() )    
        .def( bp::self |= bp::self )    
        .def( bp::self |= bp::other< QDir::Filter >() )    
        .def( ~bp::self )    
        .def( "__str__", &pvt_get_name);

}
