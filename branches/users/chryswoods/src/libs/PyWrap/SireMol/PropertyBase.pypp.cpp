// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "PropertyBase.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_PropertyBase_class(){

    bp::class_< SireMol::PropertyBase, boost::noncopyable >( "PropertyBase", bp::no_init )    
        .def( 
            "assertCompatibleWith"
            , (void ( ::SireMol::PropertyBase::* )( ::SireMol::Molecule const & ) const)( &::SireMol::PropertyBase::assertCompatibleWith )
            , ( bp::arg("molecule") ) )    
        .def( 
            "isCompatibleWith"
            , (bool ( ::SireMol::PropertyBase::* )( ::SireMol::Molecule const & ) const)( &::SireMol::PropertyBase::isCompatibleWith )
            , ( bp::arg("molecule") ) )    
        .def( 
            "null_property"
            , (::SireMol::Property (*)(  ))( &::SireMol::PropertyBase::null_property ) )    
        .def( 
            "what"
            , (char const * ( ::SireMol::PropertyBase::* )(  ) const)( &::SireMol::PropertyBase::what ) )    
        .staticmethod( "null_property" );

}
