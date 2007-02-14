// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ArithmeticCombiningRules.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

void register_ArithmeticCombiningRules_class(){

    bp::class_< SireMM::ArithmeticCombiningRules, bp::bases< SireMM::CombiningRuleBase > >( "ArithmeticCombiningRules" )    
        .def( bp::init< >() )    
        .def( 
            "combine"
            , (void ( ::SireMM::ArithmeticCombiningRules::* )( ::QVector<SireMM::ChargeParameter> const &,::QVector<SireMM::LJParameter> const &,::QVector<SireMM::ChargeParameter> const &,::QVector<SireMM::LJParameter> const &,::SireMM::CLJPairMatrix & ) const)( &::SireMM::ArithmeticCombiningRules::combine )
            , ( bp::arg("chg0"), bp::arg("lj0"), bp::arg("chg1"), bp::arg("lj1"), bp::arg("cljmatrix") ) )    
        .def( 
            "combine"
            , (void ( ::SireMM::ArithmeticCombiningRules::* )( ::QVector<SireMM::ChargeParameter> const &,::QVector<SireMM::LJParameter> const &,::SireMM::CLJPairMatrix & ) const)( &::SireMM::ArithmeticCombiningRules::combine )
            , ( bp::arg("chgs"), bp::arg("ljs"), bp::arg("cljmatrix") ) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireMM::ArithmeticCombiningRules::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireMM::ArithmeticCombiningRules::* )(  ) const)( &::SireMM::ArithmeticCombiningRules::what ) )    
        .staticmethod( "typeName" );

}
