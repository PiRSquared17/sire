// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Matrix4.pypp.hpp"
#include "boost/python.hpp"
#include "siremaths_headers.h"

namespace bp = boost::python;

void register_Matrix4_class(){

    bp::class_< SireMaths::Matrix4 >( "Matrix4" )    
        .def( bp::init< >() )    
        .def( bp::init< double >(( bp::arg("diagonal_value") )) )    
        .def( bp::init< SireMaths::Matrix const & >(( bp::arg("m") )) )    
        .def( 
            "column0"
            , (::SireMaths::Quaternion ( ::SireMaths::Matrix4::* )(  ) const)( &::SireMaths::Matrix4::column0 ) )    
        .def( 
            "getPrincipalAxes"
            , (::SireMaths::Matrix4 ( ::SireMaths::Matrix4::* )(  ) const)( &::SireMaths::Matrix4::getPrincipalAxes ) )    
        .def( 
            "identity"
            , (::SireMaths::Matrix4 (*)(  ))( &::SireMaths::Matrix4::identity ) )    
        .def( bp::self != bp::self )    
        .def( bp::self += bp::self )    
        .def( bp::self -= bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "toString"
            , (::QString ( ::SireMaths::Matrix4::* )(  ) const)( &::SireMaths::Matrix4::toString ) )    
        .def( 
            "transpose"
            , (::SireMaths::Matrix4 ( ::SireMaths::Matrix4::* )(  ) const)( &::SireMaths::Matrix4::transpose ) )    
        .def( 
            "zero"
            , (::SireMaths::Matrix4 (*)(  ))( &::SireMaths::Matrix4::zero ) )    
        .staticmethod( "identity" )    
        .staticmethod( "zero" );

}
