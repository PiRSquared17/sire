// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "AtomTable.pypp.hpp"
#include "boost/python.hpp"
#include "siredb_headers.h"
#include "SireCAS/function.h"
#include "SireCAS/expression.h"

namespace bp = boost::python;

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireDB::AtomTable&){ return "SireDB::AtomTable";}

void register_AtomTable_class(){

    bp::class_< SireDB::AtomTable, bp::bases< SireDB::TableBase >, boost::noncopyable >( "AtomTable", bp::no_init )    
        .def( 
            "assertHaveParameter"
            , (::SireMol::CGAtomID const & ( ::SireDB::AtomTable::* )( ::SireMol::AtomID ) const)( &::SireDB::AtomTable::assertHaveParameter )
            , ( bp::arg("atomid") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "assertHaveParameter"
            , (::SireMol::CGAtomID const & ( ::SireDB::AtomTable::* )( ::SireMol::CGAtomID const & ) const)( &::SireDB::AtomTable::assertHaveParameter )
            , ( bp::arg("cgatomid") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "assertHaveParameter"
            , (::SireMol::CGAtomID const & ( ::SireDB::AtomTable::* )( ::SireMol::CGNumAtomID const & ) const)( &::SireDB::AtomTable::assertHaveParameter )
            , ( bp::arg("cgatomid") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "assertHaveParameter"
            , (::SireMol::CGAtomID const & ( ::SireDB::AtomTable::* )( ::SireMol::ResNumAtomID const & ) const)( &::SireDB::AtomTable::assertHaveParameter )
            , ( bp::arg("resatomid") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "assertHaveParameter"
            , (::SireMol::CGAtomID const & ( ::SireDB::AtomTable::* )( ::SireMol::ResIDAtomID const & ) const)( &::SireDB::AtomTable::assertHaveParameter )
            , ( bp::arg("resatomid") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "assertHaveParameter"
            , (::SireMol::CGAtomID const & ( ::SireDB::AtomTable::* )( ::SireMol::AtomIndex const & ) const)( &::SireDB::AtomTable::assertHaveParameter )
            , ( bp::arg("atomindex") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "assignedParameter"
            , (bool ( ::SireDB::AtomTable::* )( ::SireMol::CGAtomID const & ) const)( &::SireDB::AtomTable::assignedParameter )
            , ( bp::arg("cgatomid") ) )    
        .def( 
            "assignedParameter"
            , (bool ( ::SireDB::AtomTable::* )( ::SireMol::CGNumAtomID const & ) const)( &::SireDB::AtomTable::assignedParameter )
            , ( bp::arg("cgatomid") ) )    
        .def( 
            "assignedParameter"
            , (bool ( ::SireDB::AtomTable::* )( ::SireMol::ResNumAtomID const & ) const)( &::SireDB::AtomTable::assignedParameter )
            , ( bp::arg("resatomid") ) )    
        .def( 
            "assignedParameter"
            , (bool ( ::SireDB::AtomTable::* )( ::SireMol::ResIDAtomID const & ) const)( &::SireDB::AtomTable::assignedParameter )
            , ( bp::arg("resatomid") ) )    
        .def( 
            "assignedParameter"
            , (bool ( ::SireDB::AtomTable::* )( ::SireMol::AtomIndex const & ) const)( &::SireDB::AtomTable::assignedParameter )
            , ( bp::arg("atom") ) )    
        .def( 
            "assignedParameter"
            , (bool ( ::SireDB::AtomTable::* )( ::SireMol::AtomID ) const)( &::SireDB::AtomTable::assignedParameter )
            , ( bp::arg("atomid") ) )    
        .def( 
            "clear"
            , (void ( ::SireDB::AtomTable::* )(  ) )( &::SireDB::AtomTable::clear ) )    
        .def( 
            "clear"
            , (void ( ::SireDB::AtomTable::* )( ::SireMol::CGAtomID const & ) )( &::SireDB::AtomTable::clear )
            , ( bp::arg("cgatomid") ) )    
        .def( 
            "clear"
            , (void ( ::SireDB::AtomTable::* )( ::SireMol::CGNumAtomID const & ) )( &::SireDB::AtomTable::clear )
            , ( bp::arg("cgatomid") ) )    
        .def( 
            "clear"
            , (void ( ::SireDB::AtomTable::* )( ::SireMol::ResNumAtomID const & ) )( &::SireDB::AtomTable::clear )
            , ( bp::arg("resatomid") ) )    
        .def( 
            "clear"
            , (void ( ::SireDB::AtomTable::* )( ::SireMol::ResIDAtomID const & ) )( &::SireDB::AtomTable::clear )
            , ( bp::arg("resatomid") ) )    
        .def( 
            "clear"
            , (void ( ::SireDB::AtomTable::* )( ::SireMol::AtomIndex const & ) )( &::SireDB::AtomTable::clear )
            , ( bp::arg("atom") ) )    
        .def( 
            "clear"
            , (void ( ::SireDB::AtomTable::* )( ::SireMol::AtomID ) )( &::SireDB::AtomTable::clear )
            , ( bp::arg("atomid") ) )    
        .def( 
            "clear"
            , (void ( ::SireDB::AtomTable::* )( ::SireMol::ResNum ) )( &::SireDB::AtomTable::clear )
            , ( bp::arg("resnum") ) )    
        .def( 
            "clear"
            , (void ( ::SireDB::AtomTable::* )( ::SireMol::ResID ) )( &::SireDB::AtomTable::clear )
            , ( bp::arg("resid") ) )    
        .def( 
            "clear"
            , (void ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupNum ) )( &::SireDB::AtomTable::clear )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "clear"
            , (void ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupID ) )( &::SireDB::AtomTable::clear )
            , ( bp::arg("cgid") ) )    
        .def( 
            "hasMissingParameters"
            , (bool ( ::SireDB::AtomTable::* )(  ) const)( &::SireDB::AtomTable::hasMissingParameters ) )    
        .def( 
            "hasMissingParameters"
            , (bool ( ::SireDB::AtomTable::* )( ::SireMol::ResNum ) const)( &::SireDB::AtomTable::hasMissingParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "hasMissingParameters"
            , (bool ( ::SireDB::AtomTable::* )( ::SireMol::ResID ) const)( &::SireDB::AtomTable::hasMissingParameters )
            , ( bp::arg("resid") ) )    
        .def( 
            "hasMissingParameters"
            , (bool ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupID ) const)( &::SireDB::AtomTable::hasMissingParameters )
            , ( bp::arg("cgid") ) )    
        .def( 
            "hasMissingParameters"
            , (bool ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupNum ) const)( &::SireDB::AtomTable::hasMissingParameters )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "isEmpty"
            , (bool ( ::SireDB::AtomTable::* )(  ) const)( &::SireDB::AtomTable::isEmpty ) )    
        .def( 
            "missingParameters"
            , (::QSet<SireMol::AtomIndex> ( ::SireDB::AtomTable::* )(  ) const)( &::SireDB::AtomTable::missingParameters ) )    
        .def( 
            "missingParameters"
            , (::QSet<SireMol::AtomIndex> ( ::SireDB::AtomTable::* )( ::SireMol::ResNum ) const)( &::SireDB::AtomTable::missingParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "missingParameters"
            , (::QSet<SireMol::AtomIndex> ( ::SireDB::AtomTable::* )( ::SireMol::ResID ) const)( &::SireDB::AtomTable::missingParameters )
            , ( bp::arg("resid") ) )    
        .def( 
            "missingParameters"
            , (::QSet<SireMol::AtomIndex> ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupID ) const)( &::SireDB::AtomTable::missingParameters )
            , ( bp::arg("cgid") ) )    
        .def( 
            "missingParameters"
            , (::QSet<SireMol::AtomIndex> ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupNum ) const)( &::SireDB::AtomTable::missingParameters )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "nAssignedParameters"
            , (int ( ::SireDB::AtomTable::* )(  ) const)( &::SireDB::AtomTable::nAssignedParameters ) )    
        .def( 
            "nAssignedParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::ResNum ) const)( &::SireDB::AtomTable::nAssignedParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nAssignedParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::ResID ) const)( &::SireDB::AtomTable::nAssignedParameters )
            , ( bp::arg("resid") ) )    
        .def( 
            "nAssignedParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupID ) const)( &::SireDB::AtomTable::nAssignedParameters )
            , ( bp::arg("cgid") ) )    
        .def( 
            "nAssignedParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupNum ) const)( &::SireDB::AtomTable::nAssignedParameters )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "nMissingParameters"
            , (int ( ::SireDB::AtomTable::* )(  ) const)( &::SireDB::AtomTable::nMissingParameters ) )    
        .def( 
            "nMissingParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::ResNum ) const)( &::SireDB::AtomTable::nMissingParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nMissingParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::ResID ) const)( &::SireDB::AtomTable::nMissingParameters )
            , ( bp::arg("resid") ) )    
        .def( 
            "nMissingParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupID ) const)( &::SireDB::AtomTable::nMissingParameters )
            , ( bp::arg("cgid") ) )    
        .def( 
            "nMissingParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupNum ) const)( &::SireDB::AtomTable::nMissingParameters )
            , ( bp::arg("cgnum") ) )    
        .def( 
            "nParameters"
            , (int ( ::SireDB::AtomTable::* )(  ) const)( &::SireDB::AtomTable::nParameters ) )    
        .def( 
            "nParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::ResNum ) const)( &::SireDB::AtomTable::nParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::ResID ) const)( &::SireDB::AtomTable::nParameters )
            , ( bp::arg("resid") ) )    
        .def( 
            "nParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupID ) const)( &::SireDB::AtomTable::nParameters )
            , ( bp::arg("cgid") ) )    
        .def( 
            "nParameters"
            , (int ( ::SireDB::AtomTable::* )( ::SireMol::CutGroupNum ) const)( &::SireDB::AtomTable::nParameters )
            , ( bp::arg("cgnum") ) )    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireDB::AtomTable >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireDB::AtomTable >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
