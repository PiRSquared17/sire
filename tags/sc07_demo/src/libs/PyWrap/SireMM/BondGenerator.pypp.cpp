// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "BondGenerator.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

SireMM::BondGenerator __copy__(const SireMM::BondGenerator &other){ return SireMM::BondGenerator(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::BondGenerator&){ return "SireMM::BondGenerator";}

void register_BondGenerator_class(){

    bp::class_< SireMM::BondGenerator, bp::bases< SireMM::InternalGenerator<SireMM::MolBondInfo> > >( "BondGenerator" )    
        .def( bp::init< >() )    
        .def( 
            "generate"
            , (void ( ::SireMM::BondGenerator::* )( ::SireMol::Molecule const &,::SireMM::MolBondInfo & ) const)( &::SireMM::BondGenerator::generate )
            , ( bp::arg("mol"), bp::arg("bondinfo") ) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireMM::BondGenerator::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireMM::BondGenerator::* )(  ) const)( &::SireMM::BondGenerator::what ) )    
        .staticmethod( "typeName" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMM::BondGenerator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMM::BondGenerator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
