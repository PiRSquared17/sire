// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Bond.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMol/moleculedata.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

SireMol::Bond __copy__(const SireMol::Bond &other){ return SireMol::Bond(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_Bond_class(){

    bp::class_< SireMol::Bond >( "Bond" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::AtomIndex const &, SireMol::AtomIndex const & >(( bp::arg("atom0"), bp::arg("atom1") )) )    
        .def( bp::init< boost::tuples::tuple<SireMol::AtomIndex, SireMol::AtomIndex, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type> const & >(( bp::arg("atoms") )) )    
        .def( bp::init< QString const &, QString const &, SireMol::ResNum >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("resnum") )) )    
        .def( bp::init< QString const &, SireMol::ResNum, QString const &, SireMol::ResNum >(( bp::arg("nam0"), bp::arg("res0"), bp::arg("nam1"), bp::arg("res1") )) )    
        .def( 
            "at"
            , (::SireMol::AtomIndex const & ( ::SireMol::Bond::* )( int ) const)( &::SireMol::Bond::at )
            , ( bp::arg("i") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "at"
            , (::SireMol::AtomIndex const & ( ::SireMol::Bond::* )( ::SireMol::ResNum ) const)( &::SireMol::Bond::at )
            , ( bp::arg("resnum") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "atom"
            , &::SireMol::Bond::atom
            , ( bp::arg("resnum") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "atom0"
            , &::SireMol::Bond::atom0
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "atom1"
            , &::SireMol::Bond::atom1
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Bond::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::Bond::contains )
            , ( bp::arg("atom") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Bond::* )( ::QString const & ) const)( &::SireMol::Bond::contains )
            , ( bp::arg("atmnam") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::Bond::* )( ::SireMol::ResNum ) const)( &::SireMol::Bond::contains )
            , ( bp::arg("resnum") ) )    
        .def( 
            "count"
            , &::SireMol::Bond::count )    
        .def( 
            "interResidue"
            , &::SireMol::Bond::interResidue )    
        .def( 
            "intraResidue"
            , &::SireMol::Bond::intraResidue )    
        .def( bp::self != bp::self )    
        .def( bp::self < bp::self )    
        .def( bp::self <= bp::self )    
        .def( bp::self == bp::self )    
        .def( bp::self > bp::self )    
        .def( bp::self >= bp::self )    
        .def( 
            "__getitem__"
            , (::SireMol::AtomIndex const & ( ::SireMol::Bond::* )( int ) const)( &::SireMol::Bond::operator[] )
            , ( bp::arg("i") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "__getitem__"
            , (::SireMol::AtomIndex const & ( ::SireMol::Bond::* )( ::SireMol::ResNum ) const)( &::SireMol::Bond::operator[] )
            , ( bp::arg("resnum") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "other"
            , (::SireMol::AtomIndex const & ( ::SireMol::Bond::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::Bond::other )
            , ( bp::arg("atm") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "other"
            , (::SireMol::AtomIndex const & ( ::SireMol::Bond::* )( ::SireMol::ResNum ) const)( &::SireMol::Bond::other )
            , ( bp::arg("resnum") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "otherRes"
            , &::SireMol::Bond::otherRes
            , ( bp::arg("resnum") ) )    
        .def( 
            "renumber"
            , &::SireMol::Bond::renumber
            , ( bp::arg("oldnum"), bp::arg("newnum") ) )    
        .def( 
            "resNum0"
            , &::SireMol::Bond::resNum0 )    
        .def( 
            "resNum1"
            , &::SireMol::Bond::resNum1 )    
        .def( 
            "size"
            , &::SireMol::Bond::size )    
        .def( 
            "toString"
            , &::SireMol::Bond::toString )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::Bond >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::Bond >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireMol::Bond > );

}
