// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "InternalTable_MolDihedralInfo_Expression_.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>&){ return "SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>";}

void register_InternalTable_MolDihedralInfo_Expression__class(){

    bp::class_< SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>, bp::bases< SireMM::InternalTableBase<SireMM::MolDihedralInfo> >, boost::noncopyable >( "InternalTable_MolDihedralInfo_Expression_", bp::no_init )    
        .def( 
            "addInternal"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const & ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::addInternal )
            , ( bp::arg("internal") ) )    
        .def( 
            "addInternal"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const &,::SireCAS::Expression const & ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::addInternal )
            , ( bp::arg("internal"), bp::arg("param") ) )    
        .def( 
            "assignedParameter"
            , (bool ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::assignedParameter )
            , ( bp::arg("internal") ) )    
        .def( 
            "assignedParameter"
            , (bool ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupIndexID const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::assignedParameter )
            , ( bp::arg("id") ) )    
        .def( 
            "at"
            , (::SireCAS::Expression const & ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::at )
            , ( bp::arg("internal") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "at"
            , (::SireCAS::Expression const & ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupIndexID const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::at )
            , ( bp::arg("id") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "clear"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::clear ) )    
        .def( 
            "clear"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::clear )
            , ( bp::arg("resnum") ) )    
        .def( 
            "clear"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupID ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::clear )
            , ( bp::arg("groupid") ) )    
        .def( 
            "clear"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const & ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::clear )
            , ( bp::arg("internal") ) )    
        .def( 
            "clearInter"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::clearInter ) )    
        .def( 
            "clearInter"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::clearInter )
            , ( bp::arg("resnum") ) )    
        .def( 
            "clearIntra"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::clearIntra ) )    
        .def( 
            "clearIntra"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::clearIntra )
            , ( bp::arg("resnum") ) )    
        .def( 
            "getParameter"
            , (::SireCAS::Expression const & ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::getParameter )
            , ( bp::arg("internal") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "getParameter"
            , (::SireCAS::Expression const & ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupIndexID const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::getParameter )
            , ( bp::arg("id") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "hasMissingInterParameters"
            , (bool ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::hasMissingInterParameters ) )    
        .def( 
            "hasMissingInterParameters"
            , (bool ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::hasMissingInterParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "hasMissingIntraParameters"
            , (bool ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::hasMissingIntraParameters ) )    
        .def( 
            "hasMissingIntraParameters"
            , (bool ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::hasMissingIntraParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "hasMissingParameters"
            , (bool ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::hasMissingParameters ) )    
        .def( 
            "hasMissingParameters"
            , (bool ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::hasMissingParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "hasMissingParameters"
            , (bool ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupID ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::hasMissingParameters )
            , ( bp::arg("id") ) )    
        .def( 
            "interParameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::interParameters ) )    
        .def( 
            "interParameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::interParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "interParametersByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireCAS::Expression> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::interParametersByResidue ) )    
        .def( 
            "intraParameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::intraParameters ) )    
        .def( 
            "intraParameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::intraParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "intraParametersByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireCAS::Expression> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::intraParametersByResidue ) )    
        .def( 
            "missingInterParameters"
            , (::QSet<SireMol::Dihedral> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::missingInterParameters ) )    
        .def( 
            "missingInterParameters"
            , (::QSet<SireMol::Dihedral> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::missingInterParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "missingIntraParameters"
            , (::QSet<SireMol::Dihedral> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::missingIntraParameters ) )    
        .def( 
            "missingIntraParameters"
            , (::QSet<SireMol::Dihedral> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::missingIntraParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "missingParameters"
            , (::QSet<SireMol::Dihedral> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::missingParameters ) )    
        .def( 
            "missingParameters"
            , (::QSet<SireMol::Dihedral> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::missingParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "missingParameters"
            , (::QSet<SireMol::Dihedral> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::QSet<SireMol::GroupID> const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::missingParameters )
            , ( bp::arg("groupids") ) )    
        .def( 
            "nInterParameters"
            , (int ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::nInterParameters ) )    
        .def( 
            "nInterParameters"
            , (int ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::nInterParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nIntraParameters"
            , (int ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::nIntraParameters ) )    
        .def( 
            "nIntraParameters"
            , (int ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::nIntraParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nParameters"
            , (int ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::nParameters ) )    
        .def( 
            "nParameters"
            , (int ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::nParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nParameters"
            , (int ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupID ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::nParameters )
            , ( bp::arg("groupid") ) )    
        .def( 
            "__getitem__"
            , (::SireCAS::Expression const & ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::operator[] )
            , ( bp::arg("internal") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "__getitem__"
            , (::SireCAS::Expression const & ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupIndexID const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::operator[] )
            , ( bp::arg("id") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "__getitem__"
            , (::SireCAS::Expression & ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const & ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::operator[] )
            , ( bp::arg("internal") )
            , bp::return_internal_reference< 1 >() )    
        .def( 
            "__getitem__"
            , (::SireCAS::Expression & ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupIndexID const & ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::operator[] )
            , ( bp::arg("id") )
            , bp::return_internal_reference< 1 >() )    
        .def( 
            "parameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::parameters ) )    
        .def( 
            "parameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::parameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "parameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupID ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::parameters )
            , ( bp::arg("groupid") ) )    
        .def( 
            "parametersByGroup"
            , (::SireMol::GroupedVector<SireMol::GroupIndexID,SireCAS::Expression> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::parametersByGroup ) )    
        .def( 
            "parametersByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireCAS::Expression> ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::parametersByResidue ) )    
        .def( 
            "removeInterInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeInterInternals ) )    
        .def( 
            "removeInterInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeInterInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeInternal"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const & ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeInternal )
            , ( bp::arg("internal") ) )    
        .def( 
            "removeInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeInternals ) )    
        .def( 
            "removeInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupID ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeInternals )
            , ( bp::arg("groupid") ) )    
        .def( 
            "removeInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::QSet<SireMol::GroupID> const & ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeInternals )
            , ( bp::arg("groupids") ) )    
        .def( 
            "removeIntraInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeIntraInternals ) )    
        .def( 
            "removeIntraInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeIntraInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeMissingInterInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeMissingInterInternals ) )    
        .def( 
            "removeMissingInterInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeMissingInterInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeMissingInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeMissingInternals ) )    
        .def( 
            "removeMissingInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeMissingInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeMissingInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupID ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeMissingInternals )
            , ( bp::arg("groupid") ) )    
        .def( 
            "removeMissingInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::QSet<SireMol::GroupID> const & ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeMissingInternals )
            , ( bp::arg("groupids") ) )    
        .def( 
            "removeMissingIntraInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )(  ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeMissingIntraInternals ) )    
        .def( 
            "removeMissingIntraInternals"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::removeMissingIntraInternals )
            , ( bp::arg("resnum") ) )    
        .def( 
            "setParameter"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const &,::SireCAS::Expression const & ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::setParameter )
            , ( bp::arg("internal"), bp::arg("param") ) )    
        .def( 
            "setParameter"
            , (void ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupIndexID const &,::SireCAS::Expression const & ) )( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::setParameter )
            , ( bp::arg("id"), bp::arg("param") ) )    
        .def( 
            "value"
            , (::SireCAS::Expression ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::value )
            , ( bp::arg("internal") ) )    
        .def( 
            "value"
            , (::SireCAS::Expression ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::Dihedral const &,::SireCAS::Expression const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::value )
            , ( bp::arg("internal"), bp::arg("defaultParam") ) )    
        .def( 
            "value"
            , (::SireCAS::Expression ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupIndexID const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::value )
            , ( bp::arg("id") ) )    
        .def( 
            "value"
            , (::SireCAS::Expression ( ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::* )( ::SireMol::GroupIndexID const &,::SireCAS::Expression const & ) const)( &::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression>::value )
            , ( bp::arg("id"), bp::arg("defaultParam") ) )    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMM::InternalTable<SireMM::MolDihedralInfo, SireCAS::Expression> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
