// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Vector.pypp.hpp"
#include "boost/python.hpp"
#include "siremaths_headers.h"

namespace bp = boost::python;

SireMaths::Vector __copy__(const SireMaths::Vector &other){ return SireMaths::Vector(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_Vector_class(){

    bp::class_< SireMaths::Vector >( "Vector", bp::init< bp::optional< double > >(( bp::arg("val")=0.0 )) )    
        .def( bp::init< double, double, double >(( bp::arg("xpos"), bp::arg("ypos"), bp::arg("zpos") )) )    
        .def( bp::init< boost::tuples::tuple<double, double, double, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type> const & >(( bp::arg("pos") )) )    
        .def( bp::init< QString const & >(( bp::arg("str") )) )    
        .def( 
            "angle"
            , (::SireMaths::Angle (*)( ::SireMaths::Vector const &,::SireMaths::Vector const & ))( &::SireMaths::Vector::angle )
            , ( bp::arg("v0"), bp::arg("v1") ) )    
        .def( 
            "angle"
            , (::SireMaths::Angle (*)( ::SireMaths::Vector const &,::SireMaths::Vector const &,::SireMaths::Vector const & ))( &::SireMaths::Vector::angle )
            , ( bp::arg("v0"), bp::arg("v1"), bp::arg("v2") ) )    
        .def( 
            "at"
            , (double const & ( ::SireMaths::Vector::* )( unsigned int ) const)( &::SireMaths::Vector::at )
            , ( bp::arg("i") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "b"
            , (double ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::b ) )    
        .def( 
            "bearing"
            , (::SireMaths::Angle ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::bearing ) )    
        .def( 
            "bearingXY"
            , (::SireMaths::Angle ( ::SireMaths::Vector::* )( ::SireMaths::Vector const & ) const)( &::SireMaths::Vector::bearingXY )
            , ( bp::arg("v") ) )    
        .def( 
            "bearingXZ"
            , (::SireMaths::Angle ( ::SireMaths::Vector::* )( ::SireMaths::Vector const & ) const)( &::SireMaths::Vector::bearingXZ )
            , ( bp::arg("v") ) )    
        .def( 
            "bearingYZ"
            , (::SireMaths::Angle ( ::SireMaths::Vector::* )( ::SireMaths::Vector const & ) const)( &::SireMaths::Vector::bearingYZ )
            , ( bp::arg("v") ) )    
        .def( 
            "count"
            , (unsigned int ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::count ) )    
        .def( 
            "cross"
            , (::SireMaths::Vector (*)( ::SireMaths::Vector const &,::SireMaths::Vector const & ))( &::SireMaths::Vector::cross )
            , ( bp::arg("v0"), bp::arg("v1") ) )    
        .def( 
            "dihedral"
            , (::SireMaths::Angle (*)( ::SireMaths::Vector const &,::SireMaths::Vector const &,::SireMaths::Vector const &,::SireMaths::Vector const & ))( &::SireMaths::Vector::dihedral )
            , ( bp::arg("v0"), bp::arg("v1"), bp::arg("v2"), bp::arg("v3") ) )    
        .def( 
            "distance"
            , (double (*)( ::SireMaths::Vector const &,::SireMaths::Vector const & ))( &::SireMaths::Vector::distance )
            , ( bp::arg("v1"), bp::arg("v2") ) )    
        .def( 
            "distance2"
            , (double (*)( ::SireMaths::Vector const &,::SireMaths::Vector const & ))( &::SireMaths::Vector::distance2 )
            , ( bp::arg("v1"), bp::arg("v2") ) )    
        .def( 
            "dot"
            , (double (*)( ::SireMaths::Vector const &,::SireMaths::Vector const & ))( &::SireMaths::Vector::dot )
            , ( bp::arg("v0"), bp::arg("v1") ) )    
        .def( 
            "fromString"
            , (::SireMaths::Vector (*)( ::QString const & ))( &::SireMaths::Vector::fromString )
            , ( bp::arg("str") ) )    
        .def( 
            "g"
            , (double ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::g ) )    
        .def( 
            "generate"
            , (::SireMaths::Vector (*)( double,::SireMaths::Vector const &,::SireMaths::Angle const &,::SireMaths::Vector const &,::SireMaths::Angle const &,::SireMaths::Vector const & ))( &::SireMaths::Vector::generate )
            , ( bp::arg("dst"), bp::arg("v1"), bp::arg("ang"), bp::arg("v2"), bp::arg("dih"), bp::arg("v3") ) )    
        .def( 
            "invDistance"
            , (double (*)( ::SireMaths::Vector const &,::SireMaths::Vector const & ))( &::SireMaths::Vector::invDistance )
            , ( bp::arg("v1"), bp::arg("v2") ) )    
        .def( 
            "invDistance2"
            , (double (*)( ::SireMaths::Vector const &,::SireMaths::Vector const & ))( &::SireMaths::Vector::invDistance2 )
            , ( bp::arg("v1"), bp::arg("v2") ) )    
        .def( 
            "invLength"
            , (double ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::invLength ) )    
        .def( 
            "invLength2"
            , (double ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::invLength2 ) )    
        .def( 
            "isZero"
            , (bool ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::isZero ) )    
        .def( 
            "length"
            , (double ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::length ) )    
        .def( 
            "length2"
            , (double ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::length2 ) )    
        .def( 
            "manhattanLength"
            , (double ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::manhattanLength ) )    
        .def( 
            "max"
            , (::SireMaths::Vector ( ::SireMaths::Vector::* )( ::SireMaths::Vector const & ) const)( &::SireMaths::Vector::max )
            , ( bp::arg("other") ) )    
        .def( 
            "metricTensor"
            , (::SireMaths::Matrix ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::metricTensor ) )    
        .def( 
            "min"
            , (::SireMaths::Vector ( ::SireMaths::Vector::* )( ::SireMaths::Vector const & ) const)( &::SireMaths::Vector::min )
            , ( bp::arg("other") ) )    
        .def( 
            "normalise"
            , (::SireMaths::Vector ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::normalise ) )    
        .def( bp::self != bp::self )    
        .def( bp::self *= bp::other< double >() )    
        .def( bp::self += bp::self )    
        .def( -bp::self )    
        .def( bp::self -= bp::self )    
        .def( bp::self /= bp::other< double >() )    
        .def( bp::self == bp::self )    
        .def( 
            "__getitem__"
            , (double const & ( ::SireMaths::Vector::* )( unsigned int ) const)( &::SireMaths::Vector::operator[] )
            , ( bp::arg("i") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "r"
            , (double ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::r ) )    
        .def( 
            "set"
            , (void ( ::SireMaths::Vector::* )( double,double,double ) )( &::SireMaths::Vector::set )
            , ( bp::arg("x"), bp::arg("y"), bp::arg("z") ) )    
        .def( 
            "set"
            , (void ( ::SireMaths::Vector::* )( int,double const & ) )( &::SireMaths::Vector::set )
            , ( bp::arg("i"), bp::arg("v") ) )    
        .def( 
            "setB"
            , (void ( ::SireMaths::Vector::* )( double ) )( &::SireMaths::Vector::setB )
            , ( bp::arg("z") ) )    
        .def( 
            "setG"
            , (void ( ::SireMaths::Vector::* )( double ) )( &::SireMaths::Vector::setG )
            , ( bp::arg("y") ) )    
        .def( 
            "setMax"
            , (void ( ::SireMaths::Vector::* )( ::SireMaths::Vector const & ) )( &::SireMaths::Vector::setMax )
            , ( bp::arg("other") ) )    
        .def( 
            "setMin"
            , (void ( ::SireMaths::Vector::* )( ::SireMaths::Vector const & ) )( &::SireMaths::Vector::setMin )
            , ( bp::arg("other") ) )    
        .def( 
            "setR"
            , (void ( ::SireMaths::Vector::* )( double ) )( &::SireMaths::Vector::setR )
            , ( bp::arg("x") ) )    
        .def( 
            "setX"
            , (void ( ::SireMaths::Vector::* )( double ) )( &::SireMaths::Vector::setX )
            , ( bp::arg("x") ) )    
        .def( 
            "setY"
            , (void ( ::SireMaths::Vector::* )( double ) )( &::SireMaths::Vector::setY )
            , ( bp::arg("y") ) )    
        .def( 
            "setZ"
            , (void ( ::SireMaths::Vector::* )( double ) )( &::SireMaths::Vector::setZ )
            , ( bp::arg("z") ) )    
        .def( 
            "toString"
            , (::QString ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::toString ) )    
        .def( 
            "x"
            , (double ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::x ) )    
        .def( 
            "y"
            , (double ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::y ) )    
        .def( 
            "z"
            , (double ( ::SireMaths::Vector::* )(  ) const)( &::SireMaths::Vector::z ) )    
        .staticmethod( "angle" )    
        .staticmethod( "cross" )    
        .staticmethod( "dihedral" )    
        .staticmethod( "distance" )    
        .staticmethod( "distance2" )    
        .staticmethod( "dot" )    
        .staticmethod( "fromString" )    
        .staticmethod( "generate" )    
        .staticmethod( "invDistance" )    
        .staticmethod( "invDistance2" )    
        .def(self + self)    
        .def(self - self)    
        .def(self * other<double>())    
        .def(self / other<double>())    
        .def(self += self)    
        .def(self -= self)    
        .def(self *= other<double>())    
        .def(self /= other<double>())    
        .def(-self)    
        .def(self == self)    
        .def(self != self)    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMaths::Vector >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMaths::Vector >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireMaths::Vector > );

}
