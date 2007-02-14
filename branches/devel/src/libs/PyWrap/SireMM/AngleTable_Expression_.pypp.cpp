// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "AngleTable_Expression_.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

void register_AngleTable_Expression__class(){

    bp::class_< SireMM::AngleTableT<SireCAS::Expression>, bp::bases< SireMM::InternalTable<SireMM::MolAngleInfo, SireCAS::Expression> >, boost::noncopyable >( "AngleTable_Expression_", bp::no_init )    
        .def( 
            "addAngle"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::Angle const & ) )( &::SireMM::AngleTableT<SireCAS::Expression>::addAngle )
            , ( bp::arg("angle") ) )    
        .def( 
            "addAngle"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::Angle const &,::SireCAS::Expression const & ) )( &::SireMM::AngleTableT<SireCAS::Expression>::addAngle )
            , ( bp::arg("angle"), bp::arg("param") ) )    
        .def( 
            "angleParameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::angleParameters ) )    
        .def( 
            "angleParameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::angleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "angleParameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::GroupID ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::angleParameters )
            , ( bp::arg("groupid") ) )    
        .def( 
            "angleParametersByGroup"
            , (::SireMol::GroupedVector<SireMol::GroupIndexID,SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::angleParametersByGroup ) )    
        .def( 
            "angleParametersByGroup"
            , (::SireMol::GroupedVector<SireMol::GroupIndexID,SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::QSet<SireMol::GroupID> const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::angleParametersByGroup )
            , ( bp::arg("gids") ) )    
        .def( 
            "angleParametersByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::angleParametersByResidue ) )    
        .def( 
            "angleParametersByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::QSet<SireMol::ResNum> const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::angleParametersByResidue )
            , ( bp::arg("resnums") ) )    
        .def( 
            "angles"
            , (::QVector<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::angles ) )    
        .def( 
            "angles"
            , (::QVector<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::angles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "angles"
            , (::QVector<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::GroupID ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::angles )
            , ( bp::arg("groupid") ) )    
        .def( 
            "anglesByGroup"
            , (::SireMol::GroupedVector<SireMol::GroupIndexID,SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::anglesByGroup ) )    
        .def( 
            "anglesByGroup"
            , (::SireMol::GroupedVector<SireMol::GroupIndexID,SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::QSet<SireMol::GroupID> const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::anglesByGroup )
            , ( bp::arg("gids") ) )    
        .def( 
            "anglesByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::anglesByResidue ) )    
        .def( 
            "anglesByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::QSet<SireMol::ResNum> const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::anglesByResidue )
            , ( bp::arg("resnums") ) )    
        .def( 
            "assignedAngleParameter"
            , (bool ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::Angle const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::assignedAngleParameter )
            , ( bp::arg("angle") ) )    
        .def( 
            "assignedAngleParameter"
            , (bool ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::GroupIndexID const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::assignedAngleParameter )
            , ( bp::arg("id") ) )    
        .def( 
            "clearAngle"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::Angle const & ) )( &::SireMM::AngleTableT<SireCAS::Expression>::clearAngle )
            , ( bp::arg("angle") ) )    
        .def( 
            "clearAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) )( &::SireMM::AngleTableT<SireCAS::Expression>::clearAngles ) )    
        .def( 
            "clearAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::AngleTableT<SireCAS::Expression>::clearAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "clearAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::GroupID ) )( &::SireMM::AngleTableT<SireCAS::Expression>::clearAngles )
            , ( bp::arg("groupid") ) )    
        .def( 
            "clearInterAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) )( &::SireMM::AngleTableT<SireCAS::Expression>::clearInterAngles ) )    
        .def( 
            "clearInterAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::AngleTableT<SireCAS::Expression>::clearInterAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "clearIntraAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) )( &::SireMM::AngleTableT<SireCAS::Expression>::clearIntraAngles ) )    
        .def( 
            "clearIntraAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::AngleTableT<SireCAS::Expression>::clearIntraAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "getAngleParameter"
            , (::SireCAS::Expression const & ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::Angle const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::getAngleParameter )
            , ( bp::arg("angle") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "getAngleParameter"
            , (::SireCAS::Expression const & ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::GroupIndexID const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::getAngleParameter )
            , ( bp::arg("id") )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "hasMissingAngleParameters"
            , (bool ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::hasMissingAngleParameters ) )    
        .def( 
            "hasMissingAngleParameters"
            , (bool ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::hasMissingAngleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "hasMissingAngleParameters"
            , (bool ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::GroupID ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::hasMissingAngleParameters )
            , ( bp::arg("id") ) )    
        .def( 
            "hasMissingInterAngleParameters"
            , (bool ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::hasMissingInterAngleParameters ) )    
        .def( 
            "hasMissingInterAngleParameters"
            , (bool ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::hasMissingInterAngleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "hasMissingIntraAngleParameters"
            , (bool ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::hasMissingIntraAngleParameters ) )    
        .def( 
            "hasMissingIntraAngleParameters"
            , (bool ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::hasMissingIntraAngleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "interAngleParameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::interAngleParameters ) )    
        .def( 
            "interAngleParameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::interAngleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "interAngleParametersByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::interAngleParametersByResidue ) )    
        .def( 
            "interAngleParametersByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::QSet<SireMol::ResNum> const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::interAngleParametersByResidue )
            , ( bp::arg("resnums") ) )    
        .def( 
            "interAngles"
            , (::QVector<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::interAngles ) )    
        .def( 
            "interAngles"
            , (::QVector<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::interAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "interAnglesByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::interAnglesByResidue ) )    
        .def( 
            "interAnglesByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::QSet<SireMol::ResNum> const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::interAnglesByResidue )
            , ( bp::arg("resnums") ) )    
        .def( 
            "intraAngleParameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::intraAngleParameters ) )    
        .def( 
            "intraAngleParameters"
            , (::QVector<SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::intraAngleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "intraAngleParametersByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::intraAngleParametersByResidue ) )    
        .def( 
            "intraAngleParametersByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireCAS::Expression> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::QSet<SireMol::ResNum> const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::intraAngleParametersByResidue )
            , ( bp::arg("resnums") ) )    
        .def( 
            "intraAngles"
            , (::QVector<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::intraAngles ) )    
        .def( 
            "intraAngles"
            , (::QVector<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::intraAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "intraAnglesByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::intraAnglesByResidue ) )    
        .def( 
            "intraAnglesByResidue"
            , (::SireMol::GroupedVector<SireMol::ResNumIndexID,SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::QSet<SireMol::ResNum> const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::intraAnglesByResidue )
            , ( bp::arg("resnums") ) )    
        .def( 
            "missingAngleParameters"
            , (::QSet<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::missingAngleParameters ) )    
        .def( 
            "missingAngleParameters"
            , (::QSet<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::missingAngleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "missingAngleParameters"
            , (::QSet<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::QSet<SireMol::GroupID> const & ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::missingAngleParameters )
            , ( bp::arg("groupids") ) )    
        .def( 
            "missingInterAngleParameters"
            , (::QSet<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::missingInterAngleParameters ) )    
        .def( 
            "missingInterAngleParameters"
            , (::QSet<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::missingInterAngleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "missingIntraAngleParameters"
            , (::QSet<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::missingIntraAngleParameters ) )    
        .def( 
            "missingIntraAngleParameters"
            , (::QSet<SireMol::Angle> ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::missingIntraAngleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nAngleParameters"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nAngleParameters ) )    
        .def( 
            "nAngleParameters"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nAngleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nAngleParameters"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::GroupID ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nAngleParameters )
            , ( bp::arg("groupid") ) )    
        .def( 
            "nAngles"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nAngles ) )    
        .def( 
            "nAngles"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nAngles"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::GroupID ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nAngles )
            , ( bp::arg("groupid") ) )    
        .def( 
            "nInterAngleParameters"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nInterAngleParameters ) )    
        .def( 
            "nInterAngleParameters"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nInterAngleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nInterAngles"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nInterAngles ) )    
        .def( 
            "nInterAngles"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nInterAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nIntraAngleParameters"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nIntraAngleParameters ) )    
        .def( 
            "nIntraAngleParameters"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nIntraAngleParameters )
            , ( bp::arg("resnum") ) )    
        .def( 
            "nIntraAngles"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nIntraAngles ) )    
        .def( 
            "nIntraAngles"
            , (int ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) const)( &::SireMM::AngleTableT<SireCAS::Expression>::nIntraAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeAngle"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::Angle const & ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeAngle )
            , ( bp::arg("angle") ) )    
        .def( 
            "removeAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeAngles ) )    
        .def( 
            "removeAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::GroupID ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeAngles )
            , ( bp::arg("groupid") ) )    
        .def( 
            "removeAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::QSet<SireMol::GroupID> const & ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeAngles )
            , ( bp::arg("groupids") ) )    
        .def( 
            "removeInterAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeInterAngles ) )    
        .def( 
            "removeInterAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeInterAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeIntraAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeIntraAngles ) )    
        .def( 
            "removeIntraAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeIntraAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeMissingAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeMissingAngles ) )    
        .def( 
            "removeMissingAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeMissingAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeMissingAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::GroupID ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeMissingAngles )
            , ( bp::arg("groupid") ) )    
        .def( 
            "removeMissingAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::QSet<SireMol::GroupID> const & ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeMissingAngles )
            , ( bp::arg("groupids") ) )    
        .def( 
            "removeMissingInterAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeMissingInterAngles ) )    
        .def( 
            "removeMissingInterAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeMissingInterAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "removeMissingIntraAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )(  ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeMissingIntraAngles ) )    
        .def( 
            "removeMissingIntraAngles"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::ResNum ) )( &::SireMM::AngleTableT<SireCAS::Expression>::removeMissingIntraAngles )
            , ( bp::arg("resnum") ) )    
        .def( 
            "setAngleParameter"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::Angle const &,::SireCAS::Expression const & ) )( &::SireMM::AngleTableT<SireCAS::Expression>::setAngleParameter )
            , ( bp::arg("angle"), bp::arg("param") ) )    
        .def( 
            "setAngleParameter"
            , (void ( ::SireMM::AngleTableT<SireCAS::Expression>::* )( ::SireMol::GroupIndexID const &,::SireCAS::Expression const & ) )( &::SireMM::AngleTableT<SireCAS::Expression>::setAngleParameter )
            , ( bp::arg("id"), bp::arg("param") ) );

}
