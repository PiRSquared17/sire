// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Dihedral.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

SireMol::Dihedral __copy__(const SireMol::Dihedral &other){ return SireMol::Dihedral(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_Dihedral_class(){

    bp::class_< SireMol::Dihedral >( "Dihedral" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::AtomIndex const &, SireMol::AtomIndex const &, SireMol::AtomIndex const &, SireMol::AtomIndex const & >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") )) )    
        .def( bp::init< boost::tuples::tuple<SireMol::AtomIndex, SireMol::AtomIndex, SireMol::AtomIndex, SireMol::AtomIndex, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type> const & >(( bp::arg("tuple") )) )    
        .def( bp::init< QString const &, QString const &, QString const &, QString const &, SireMol::ResNum >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3"), bp::arg("resnum") )) )    
        .def( bp::init< QString const &, SireMol::ResNum, QString const &, SireMol::ResNum, QString const &, SireMol::ResNum, QString const &, SireMol::ResNum >(( bp::arg("atm0"), bp::arg("res0"), bp::arg("atm1"), bp::arg("res1"), bp::arg("atm2"), bp::arg("res2"), bp::arg("atm3"), bp::arg("res3") )) )    
        .def( 
            "at"
            , &::SireMol::Dihedral::at
            , ( bp::arg("i") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "atom0"
            , &::SireMol::Dihedral::atom0
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "atom1"
            , &::SireMol::Dihedral::atom1
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "atom2"
            , &::SireMol::Dihedral::atom2
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "atom3"
            , &::SireMol::Dihedral::atom3
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Dihedral::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::Dihedral::contains )
            , ( bp::arg("atom") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Dihedral::* )( ::QString const & ) const)( &::SireMol::Dihedral::contains )
            , ( bp::arg("atmnam") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Dihedral::* )( ::SireMol::ResNum ) const)( &::SireMol::Dihedral::contains )
            , ( bp::arg("resnum") ) )    
        .def( 
            "count"
            , &::SireMol::Dihedral::count )    
        .def( 
            "interResidue"
            , &::SireMol::Dihedral::interResidue )    
        .def( 
            "intraResidue"
            , &::SireMol::Dihedral::intraResidue )    
        .def( bp::self != bp::self )    
        .def( bp::self < bp::self )    
        .def( bp::self <= bp::self )    
        .def( bp::self == bp::self )    
        .def( bp::self > bp::self )    
        .def( bp::self >= bp::self )    
        .def( 
            "__getitem__"
            , &::SireMol::Dihedral::operator[]
            , ( bp::arg("i") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "resNum0"
            , &::SireMol::Dihedral::resNum0 )    
        .def( 
            "resNum1"
            , &::SireMol::Dihedral::resNum1 )    
        .def( 
            "resNum2"
            , &::SireMol::Dihedral::resNum2 )    
        .def( 
            "resNum3"
            , &::SireMol::Dihedral::resNum3 )    
        .def( 
            "size"
            , &::SireMol::Dihedral::size )    
        .def( 
            "toString"
            , &::SireMol::Dihedral::toString )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::Dihedral >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::Dihedral >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireMol::Dihedral > );

}
