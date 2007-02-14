// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "BondDB.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

void register_BondDB_class(){

    bp::class_< SireMM::BondDB, bp::bases< SireDB::Term12DB, SireMM::InternalDB<SireMM::MolBondInfo> > >( "BondDB" )    
        .def( bp::init< >() )    
        .def( 
            "addBond"
            , (void ( ::SireMM::BondDB::* )( ::QString const &,::SireCAS::Expression const & ) )( &::SireMM::BondDB::addBond )
            , ( bp::arg("userid"), bp::arg("bondfunc") ) )    
        .def( 
            "assignParameter"
            , (bool ( ::SireMM::BondDB::* )( ::SireMol::Bond const &,::SireDB::RelateIDMap const &,::SireDB::ParameterTable & ) )( &::SireMM::BondDB::assignParameter )
            , ( bp::arg("bond"), bp::arg("relateids"), bp::arg("param_table") ) )    
        .def( 
            "getBond"
            , (::SireCAS::Expression ( ::SireMM::BondDB::* )( ::QString const &,bool * ) )( &::SireMM::BondDB::getBond )
            , ( bp::arg("userid"), bp::arg("foundbond")=bp::object() ) )    
        .def( 
            "getBond"
            , (::SireCAS::Expression ( ::SireMM::BondDB::* )( ::SireDB::RelateID,bool * ) )( &::SireMM::BondDB::getBond )
            , ( bp::arg("relateid"), bp::arg("foundbond")=bp::object() ) )    
        .def( 
            "getBond"
            , (::SireCAS::Expression ( ::SireMM::BondDB::* )( ::SireDB::RelateIDMap const &,bool * ) )( &::SireMM::BondDB::getBond )
            , ( bp::arg("relateids"), bp::arg("foundfunc")=bp::object() ) )    
        .def( 
            "r"
            , (::SireCAS::Symbol const & ( ::SireMM::BondDB::* )(  ) const)( &::SireMM::BondDB::r )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "relateBond"
            , (void ( ::SireMM::BondDB::* )( ::SireDB::AssertMatch<2> const &,::QString const & ) )( &::SireMM::BondDB::relateBond )
            , ( bp::arg("matchbond"), bp::arg("userid") ) )    
        .def( 
            "relateBond"
            , (void ( ::SireMM::BondDB::* )( ::SireDB::AssertMatch<2> const &,::SireCAS::Expression const & ) )( &::SireMM::BondDB::relateBond )
            , ( bp::arg("matchbond"), bp::arg("bondfunc") ) )    
        .def( 
            "relateBond"
            , (void ( ::SireMM::BondDB::* )( ::SireDB::RelateID,::QString const & ) )( &::SireMM::BondDB::relateBond )
            , ( bp::arg("relateid"), bp::arg("userid") ) )    
        .def( 
            "relateBond"
            , (void ( ::SireMM::BondDB::* )( ::SireDB::RelateID,::SireCAS::Expression const & ) )( &::SireMM::BondDB::relateBond )
            , ( bp::arg("relateid"), bp::arg("bondfunc") ) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireMM::BondDB::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireMM::BondDB::* )(  ) const)( &::SireMM::BondDB::what ) )    
        .staticmethod( "typeName" );

}
