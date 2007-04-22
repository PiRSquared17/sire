// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MoleculeBonds.pypp.hpp"
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

SireMol::MoleculeBonds __copy__(const SireMol::MoleculeBonds &other){ return SireMol::MoleculeBonds(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_MoleculeBonds_class(){

    bp::class_< SireMol::MoleculeBonds >( "MoleculeBonds" )    
        .def( bp::init< >() )    
        .def( 
            "add"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::AtomIndex const &,::SireMol::AtomIndex const & ) )( &::SireMol::MoleculeBonds::add )
            , ( bp::arg("atom0"), bp::arg("atom1") ) )    
        .def( 
            "add"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::ResNum,::QString const &,::QString const & ) )( &::SireMol::MoleculeBonds::add )
            , ( bp::arg("resnum"), bp::arg("atom0"), bp::arg("atom1") ) )    
        .def( 
            "add"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::Bond const & ) )( &::SireMol::MoleculeBonds::add )
            , ( bp::arg("bond") ) )    
        .def( 
            "atomsBondedTo"
            , &::SireMol::MoleculeBonds::atomsBondedTo
            , ( bp::arg("atom") ) )    
        .def( 
            "bonded"
            , (bool ( ::SireMol::MoleculeBonds::* )( ::SireMol::ResNum,::SireMol::ResNum ) const)( &::SireMol::MoleculeBonds::bonded )
            , ( bp::arg("resnum0"), bp::arg("resnum1") ) )    
        .def( 
            "bonded"
            , (bool ( ::SireMol::MoleculeBonds::* )( ::SireMol::AtomIndex const &,::SireMol::AtomIndex const & ) const)( &::SireMol::MoleculeBonds::bonded )
            , ( bp::arg("atom0"), bp::arg("atom1") ) )    
        .def( 
            "bondedResidues"
            , &::SireMol::MoleculeBonds::bondedResidues
            , ( bp::arg("resnum") ) )    
        .def( 
            "bonds"
            , (::QList<SireMol::Bond> ( ::SireMol::MoleculeBonds::* )(  ) const)( &::SireMol::MoleculeBonds::bonds ) )    
        .def( 
            "bonds"
            , (::QList<SireMol::Bond> ( ::SireMol::MoleculeBonds::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::MoleculeBonds::bonds )
            , ( bp::arg("atom") ) )    
        .def( 
            "bonds"
            , (::QList<SireMol::Bond> ( ::SireMol::MoleculeBonds::* )( ::SireMol::ResNum ) const)( &::SireMol::MoleculeBonds::bonds )
            , ( bp::arg("resnum") ) )    
        .def( 
            "clear"
            , &::SireMol::MoleculeBonds::clear )    
        .def( 
            "contains"
            , (bool ( ::SireMol::MoleculeBonds::* )( ::SireMol::ResNum ) const)( &::SireMol::MoleculeBonds::contains )
            , ( bp::arg("resnum") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::MoleculeBonds::* )( ::SireMol::AtomIndex const & ) const)( &::SireMol::MoleculeBonds::contains )
            , ( bp::arg("atom") ) )    
        .def( 
            "contains"
            , (bool ( ::SireMol::MoleculeBonds::* )( ::SireMol::Bond const & ) const)( &::SireMol::MoleculeBonds::contains )
            , ( bp::arg("bond") ) )    
        .def( 
            "finalise"
            , &::SireMol::MoleculeBonds::finalise )    
        .def( 
            "isEmpty"
            , &::SireMol::MoleculeBonds::isEmpty )    
        .def( 
            "nBonds"
            , &::SireMol::MoleculeBonds::nBonds )    
        .def( 
            "nInterBonds"
            , &::SireMol::MoleculeBonds::nInterBonds )    
        .def( 
            "nIntraBonds"
            , &::SireMol::MoleculeBonds::nIntraBonds )    
        .def( 
            "nResidues"
            , &::SireMol::MoleculeBonds::nResidues )    
        .def( bp::self != bp::self )    
        .def( bp::self == bp::self )    
        .def( 
            "remove"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::AtomIndex const &,::SireMol::AtomIndex const & ) )( &::SireMol::MoleculeBonds::remove )
            , ( bp::arg("atom0"), bp::arg("atom1") ) )    
        .def( 
            "remove"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::ResNum,::QString const &,::QString const & ) )( &::SireMol::MoleculeBonds::remove )
            , ( bp::arg("resnum"), bp::arg("atom0"), bp::arg("atom1") ) )    
        .def( 
            "remove"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::Bond const & ) )( &::SireMol::MoleculeBonds::remove )
            , ( bp::arg("bond") ) )    
        .def( 
            "remove"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::AtomIndex const & ) )( &::SireMol::MoleculeBonds::remove )
            , ( bp::arg("atom") ) )    
        .def( 
            "remove"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::ResNum ) )( &::SireMol::MoleculeBonds::remove )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeAll"
            , &::SireMol::MoleculeBonds::removeAll )    
        .def( 
            "removeAllInter"
            , &::SireMol::MoleculeBonds::removeAllInter )    
        .def( 
            "removeAllIntra"
            , &::SireMol::MoleculeBonds::removeAllIntra )    
        .def( 
            "removeInter"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::ResNum ) )( &::SireMol::MoleculeBonds::removeInter )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeInter"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::AtomIndex const & ) )( &::SireMol::MoleculeBonds::removeInter )
            , ( bp::arg("atom") ) )    
        .def( 
            "removeIntra"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::ResNum ) )( &::SireMol::MoleculeBonds::removeIntra )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeIntra"
            , (void ( ::SireMol::MoleculeBonds::* )( ::SireMol::AtomIndex const & ) )( &::SireMol::MoleculeBonds::removeIntra )
            , ( bp::arg("atom") ) )    
        .def( 
            "renumber"
            , &::SireMol::MoleculeBonds::renumber
            , ( bp::arg("oldnum"), bp::arg("newnum") ) )    
        .def( 
            "resNums"
            , &::SireMol::MoleculeBonds::resNums )    
        .def( 
            "resNumsBondedTo"
            , &::SireMol::MoleculeBonds::resNumsBondedTo
            , ( bp::arg("resnum") ) )    
        .def( 
            "residue"
            , &::SireMol::MoleculeBonds::residue
            , ( bp::arg("resnum") ) )    
        .def( 
            "toString"
            , &::SireMol::MoleculeBonds::toString )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::MoleculeBonds >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::MoleculeBonds >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &SirePy::__str__< ::SireMol::MoleculeBonds > );

}
