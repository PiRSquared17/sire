// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ArithmeticCombiningRules.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

SireMM::ArithmeticCombiningRules __copy__(const SireMM::ArithmeticCombiningRules &other){ return SireMM::ArithmeticCombiningRules(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::ArithmeticCombiningRules&){ return "SireMM::ArithmeticCombiningRules";}

void register_ArithmeticCombiningRules_class(){

    bp::class_< SireMM::ArithmeticCombiningRules, bp::bases< SireMM::CombiningRuleBase > >( "ArithmeticCombiningRules" )    
        .def( bp::init< >() )    
        .def( 
            "combine"
            , (void ( ::SireMM::ArithmeticCombiningRules::* )( ::QVector<SireMM::CLJParameter> const &,::QVector<SireMM::CLJParameter> const &,::SireMM::CLJPairMatrix & ) const)( &::SireMM::ArithmeticCombiningRules::combine )
            , ( bp::arg("clj0"), bp::arg("clj1"), bp::arg("cljmatrix") ) )    
        .def( 
            "combine"
            , (void ( ::SireMM::ArithmeticCombiningRules::* )( ::QVector<SireMM::CLJParameter> const &,::SireMM::CLJPairMatrix & ) const)( &::SireMM::ArithmeticCombiningRules::combine )
            , ( bp::arg("cljs"), bp::arg("cljmatrix") ) )    
        .def( 
            "combine"
            , (void ( ::SireMM::ArithmeticCombiningRules::* )( ::QVector<SireMM::ChargeParameter> const &,::QVector<SireMM::LJParameter> const &,::QVector<SireMM::ChargeParameter> const &,::QVector<SireMM::LJParameter> const &,::SireMM::CLJPairMatrix & ) const)( &::SireMM::ArithmeticCombiningRules::combine )
            , ( bp::arg("chg0"), bp::arg("lj0"), bp::arg("chg1"), bp::arg("lj1"), bp::arg("cljmatrix") ) )    
        .def( 
            "combine"
            , (void ( ::SireMM::ArithmeticCombiningRules::* )( ::QVector<SireMM::ChargeParameter> const &,::QVector<SireMM::LJParameter> const &,::SireMM::CLJPairMatrix & ) const)( &::SireMM::ArithmeticCombiningRules::combine )
            , ( bp::arg("chgs"), bp::arg("ljs"), bp::arg("cljmatrix") ) )    
        .def( 
            "combine"
            , (void ( ::SireMM::ArithmeticCombiningRules::* )( ::QVector<SireMM::LJParameter> const &,::QVector<SireMM::LJParameter> const &,::SireMM::LJPairMatrix & ) const)( &::SireMM::ArithmeticCombiningRules::combine )
            , ( bp::arg("lj0"), bp::arg("lj1"), bp::arg("ljmatrix") ) )    
        .def( 
            "combine"
            , (void ( ::SireMM::ArithmeticCombiningRules::* )( ::QVector<SireMM::LJParameter> const &,::SireMM::LJPairMatrix & ) const)( &::SireMM::ArithmeticCombiningRules::combine )
            , ( bp::arg("ljs"), bp::arg("ljmatrix") ) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireMM::ArithmeticCombiningRules::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireMM::ArithmeticCombiningRules::* )(  ) const)( &::SireMM::ArithmeticCombiningRules::what ) )    
        .staticmethod( "typeName" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMM::ArithmeticCombiningRules >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMM::ArithmeticCombiningRules >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
