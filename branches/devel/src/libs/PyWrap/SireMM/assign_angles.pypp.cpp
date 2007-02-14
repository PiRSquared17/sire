// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "assign_angles.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

#include "SireDB/using_database.h"

using namespace SireMM;

void register_assign_angles_class(){

    bp::class_< SireMM::assign_angles, bp::bases< SireMM::assign_internals<SireMM::MolAngleInfo> > >( "assign_angles", bp::init< bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("db0"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("db0"), bp::arg("db1"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("db5"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("db5"), bp::arg("db6"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("db5"), bp::arg("db6"), bp::arg("db7"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("db5"), bp::arg("db6"), bp::arg("db7"), bp::arg("db8"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("db5"), bp::arg("db6"), bp::arg("db7"), bp::arg("db8"), bp::arg("db9"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< QSet<SireMol::Angle> const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("angles"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< QSet<SireMol::Angle> const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("angles"), bp::arg("db0"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< QSet<SireMol::Angle> const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("angles"), bp::arg("db0"), bp::arg("db1"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< QSet<SireMol::Angle> const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("angles"), bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< QSet<SireMol::Angle> const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("angles"), bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< QSet<SireMol::Angle> const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("angles"), bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< QSet<SireMol::Angle> const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("angles"), bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("db5"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< QSet<SireMol::Angle> const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("angles"), bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("db5"), bp::arg("db6"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< QSet<SireMol::Angle> const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("angles"), bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("db5"), bp::arg("db6"), bp::arg("db7"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< QSet<SireMol::Angle> const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("angles"), bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("db5"), bp::arg("db6"), bp::arg("db7"), bp::arg("db8"), bp::arg("generator")=AngleGenerator() )) )    
        .def( bp::init< QSet<SireMol::Angle> const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, SireDB::using_database const &, bp::optional< SireMM::AngleGeneratorBase const & > >(( bp::arg("angles"), bp::arg("db0"), bp::arg("db1"), bp::arg("db2"), bp::arg("db3"), bp::arg("db4"), bp::arg("db5"), bp::arg("db6"), bp::arg("db7"), bp::arg("db8"), bp::arg("db9"), bp::arg("generator")=AngleGenerator() )) )    
        .def( 
            "assignParameters"
            , (void ( ::SireMM::assign_angles::* )( ::SireMol::Molecule const &,::SireDB::ParameterTable &,::SireDB::ParameterDB &,::SireDB::MatchMRData const & ) const)( &::SireMM::assign_angles::assignParameters )
            , ( bp::arg("molecule"), bp::arg("param_table"), bp::arg("database"), bp::arg("matchmr")=::SireDB::MatchMRData( ) ) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireMM::assign_angles::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireMM::assign_angles::* )(  ) const)( &::SireMM::assign_angles::what ) )    
        .staticmethod( "typeName" );

}
