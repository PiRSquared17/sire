// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MajVersion.pypp.hpp"
#include "boost/python.hpp"
#include "sirebase_headers.h"

namespace bp = boost::python;

void register_MajVersion_class(){

    bp::class_< SireBase::MajVersion >( "MajVersion" )    
        .def( bp::init< >() )    
        .def( bp::init< SireBase::Incremint * >(( bp::arg("incremint") )) )    
        .def( 
            "assertSameVersion"
            , (void ( ::SireBase::MajVersion::* )( ::SireBase::MajVersion const & ) const)( &::SireBase::MajVersion::assertSameVersion )
            , ( bp::arg("other") ) )    
        .def( 
            "increment"
            , (void ( ::SireBase::MajVersion::* )(  ) )( &::SireBase::MajVersion::increment ) )    
        .def( "__int__", &SireBase::MajVersion::operator ::quint32  )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self );

}
