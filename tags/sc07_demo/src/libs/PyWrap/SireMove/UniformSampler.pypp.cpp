// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "UniformSampler.pypp.hpp"
#include "boost/python.hpp"
#include "siremove_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/moleculegroup.h"
#include "SireSystem/simsystem.h"

namespace bp = boost::python;

SireMove::UniformSampler __copy__(const SireMove::UniformSampler &other){ return SireMove::UniformSampler(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMove::UniformSampler&){ return "SireMove::UniformSampler";}

void register_UniformSampler_class(){

    bp::class_< SireMove::UniformSampler, bp::bases< SireMove::SamplerBase > >( "UniformSampler", bp::init< bp::optional< SireMaths::RanGenerator const & > >(( bp::arg("rangenerator")=::SireMaths::RanGenerator( ) )) )    
        .def( bp::init< SireMol::MoleculeGroup const &, bp::optional< SireMaths::RanGenerator const & > >(( bp::arg("group"), bp::arg("rangenerator")=::SireMaths::RanGenerator( ) )) )    
        .def( 
            "assertCompatibleWith"
            , (void ( ::SireMove::UniformSampler::* )( ::SireSystem::QuerySystem const & ) const)( &::SireMove::UniformSampler::assertCompatibleWith )
            , ( bp::arg("system") ) )    
        .def( 
            "probabilityOf"
            , (double ( ::SireMove::UniformSampler::* )( ::SireMol::PartialMolecule const &,::SireSystem::QuerySystem const & ) )( &::SireMove::UniformSampler::probabilityOf )
            , ( bp::arg("molecule"), bp::arg("system") ) )    
        .def( 
            "sample"
            , (::boost::tuples::tuple<SireMol::PartialMolecule,double,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type> ( ::SireMove::UniformSampler::* )( ::SireSystem::QuerySystem const & ) )( &::SireMove::UniformSampler::sample )
            , ( bp::arg("system") ) )    
        .def( 
            "setGroup"
            , (void ( ::SireMove::UniformSampler::* )( ::SireMol::MoleculeGroup const & ) )( &::SireMove::UniformSampler::setGroup )
            , ( bp::arg("group") ) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireMove::UniformSampler::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireMove::UniformSampler::* )(  ) const)( &::SireMove::UniformSampler::what ) )    
        .staticmethod( "typeName" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMove::UniformSampler >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMove::UniformSampler >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
