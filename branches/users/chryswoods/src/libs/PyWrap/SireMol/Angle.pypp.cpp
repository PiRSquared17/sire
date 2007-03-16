// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Angle.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

SireMol::Angle __copy__(const SireMol::Angle &other){ return SireMol::Angle(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_Angle_class(){

    bp::class_< SireMol::Angle >( "Angle" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::AtomIndex const &, SireMol::AtomIndex const &, SireMol::AtomIndex const & >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2") )) )    
        .def( bp::init< boost::tuples::tuple<SireMol::AtomIndex, SireMol::AtomIndex, SireMol::AtomIndex, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type> const & >(( bp::arg("tuple") )) )    
        .def( bp::init< QString const &, QString const &, QString const &, SireMol::ResNum >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("resnum") )) )    
        .def( bp::init< QString const &, SireMol::ResNum, QString const &, SireMol::ResNum, QString const &, SireMol::ResNum >(( bp::arg("atm0"), bp::arg("res0"), bp::arg("atm1"), bp::arg("res1"), bp::arg("atm2"), bp::arg("res2") )) )    
        .def( 
            "at"
            , &::SireMol::Angle::at
            , ( bp::arg("i") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "atom0"
            , &::SireMol::Angle::atom0
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "atom1"
            , &::SireMol::Angle::atom1
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "atom2"
            , &::SireMol::Angle::atom2
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Angle::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::Angle::contains )
            , ( bp::arg("atom") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Angle::* )( ::QString const & ) const)( &::SireMol::Angle::contains )
            , ( bp::arg("atmnam") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Angle::* )( ::SireMol::ResNum ) const)( &::SireMol::Angle::contains )
            , ( bp::arg("resnum") ) )    
        .def( 
            "count"
            , &::SireMol::Angle::count )    
        .def( 
            "interResidue"
            , &::SireMol::Angle::interResidue )    
        .def( 
            "intraResidue"
            , &::SireMol::Angle::intraResidue )    
        .def( bp::self != bp::self )    
        .def( bp::self < bp::self )    
        .def( bp::self <= bp::self )    
        .def( bp::self == bp::self )    
        .def( bp::self > bp::self )    
        .def( bp::self >= bp::self )    
        .def( 
            "__getitem__"
            , &::SireMol::Angle::operator[]
            , ( bp::arg("i") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "resNum0"
            , &::SireMol::Angle::resNum0 )    
        .def( 
            "resNum1"
            , &::SireMol::Angle::resNum1 )    
        .def( 
            "resNum2"
            , &::SireMol::Angle::resNum2 )    
        .def( 
            "size"
            , &::SireMol::Angle::size )    
        .def( 
            "toString"
            , &::SireMol::Angle::toString )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::Angle >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::Angle >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireMol::Angle > );

}
