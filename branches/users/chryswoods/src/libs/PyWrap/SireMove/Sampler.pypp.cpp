// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Sampler.pypp.hpp"
#include "boost/python.hpp"
#include "siremove_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/moleculegroup.h"
#include "SireSystem/simsystem.h"

namespace bp = boost::python;

SireMove::Sampler __copy__(const SireMove::Sampler &other){ return SireMove::Sampler(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMove::Sampler&){ return "SireMove::Sampler";}

void register_Sampler_class(){

    bp::class_< SireMove::Sampler >( "Sampler" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMove::SamplerBase const & >(( bp::arg("sampler") )) )    
        .def( 
            "generator"
            , (::SireMaths::RanGenerator const & ( ::SireMove::Sampler::* )(  ) const)( &::SireMove::Sampler::generator )
            , bp::return_value_policy< bp::copy_const_reference >() )    
        .def( 
            "probabilityOf"
            , (double ( ::SireMove::Sampler::* )( ::SireMol::PartialMolecule const &,::SireSystem::QuerySystem const & ) )( &::SireMove::Sampler::probabilityOf )
            , ( bp::arg("molecule"), bp::arg("system") ) )    
        .def( 
            "sample"
            , (::boost::tuples::tuple<SireMol::PartialMolecule,double,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type,boost::tuples::null_type> ( ::SireMove::Sampler::* )( ::SireSystem::QuerySystem const & ) )( &::SireMove::Sampler::sample )
            , ( bp::arg("system") ) )    
        .def( 
            "setGenerator"
            , (void ( ::SireMove::Sampler::* )( ::SireMaths::RanGenerator const & ) )( &::SireMove::Sampler::setGenerator )
            , ( bp::arg("generator") ) )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMove::Sampler >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMove::Sampler >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
