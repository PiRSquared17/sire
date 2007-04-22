// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MoleculeView.pypp.hpp"
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

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMol::MoleculeView&){ return "SireMol::MoleculeView";}

void register_MoleculeView_class(){

    bp::class_< SireMol::MoleculeView, bp::bases< SireMol::MolDataView >, boost::noncopyable >( "MoleculeView", bp::no_init )    
        .def( 
            "assertSameMajorVersion"
            , &::SireMol::MoleculeView::assertSameMajorVersion
            , ( bp::arg("other") ) )    
        .def( 
            "assertSameMolecule"
            , &::SireMol::MoleculeView::assertSameMolecule
            , ( bp::arg("other") ) )    
        .def( 
            "assertSameVersion"
            , &::SireMol::MoleculeView::assertSameVersion
            , ( bp::arg("other") ) )    
        .def( 
            "extract"
            , (::SireMol::PropertyExtractor ( ::SireMol::MoleculeView::* )(  ) const)( &::SireMol::MoleculeView::extract ) )    
        .def( 
            "extract"
            , (::SireMol::PropertyExtractor ( ::SireMol::MoleculeView::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::MoleculeView::extract )
            , ( bp::arg("selected_atoms") ) )    
        .def( 
            "molecule"
            , &::SireMol::MoleculeView::molecule )    
        .def( 
            "move"
            , (::SireMol::MoleculeMover ( ::SireMol::MoleculeView::* )(  ) const)( &::SireMol::MoleculeView::move ) )    
        .def( 
            "move"
            , (::SireMol::MoleculeMover ( ::SireMol::MoleculeView::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::MoleculeView::move )
            , ( bp::arg("selection") ) )    
        .def( 
            "selectedAtoms"
            , &::SireMol::MoleculeView::selectedAtoms
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "selection"
            , (::SireMol::AtomSelector ( ::SireMol::MoleculeView::* )(  ) const)( &::SireMol::MoleculeView::selection ) )    
        .def( 
            "selection"
            , (::SireMol::AtomSelector ( ::SireMol::MoleculeView::* )( ::SireMol::SelectionFromMol const & ) const)( &::SireMol::MoleculeView::selection )
            , ( bp::arg("selected_atoms") ) )    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::MoleculeView >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::MoleculeView >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
