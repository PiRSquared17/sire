// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MapAsCutGroups.pypp.hpp"
#include "boost/python.hpp"
#include "siremove_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/moleculegroup.h"
#include "SireSystem/simsystem.h"

namespace bp = boost::python;

SireMove::MapAsCutGroups __copy__(const SireMove::MapAsCutGroups &other){ return SireMove::MapAsCutGroups(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMove::MapAsCutGroups&){ return "SireMove::MapAsCutGroups";}

void register_MapAsCutGroups_class(){

    bp::class_< SireMove::MapAsCutGroups, bp::bases< SireMove::MolMappingFunctionBase > >( "MapAsCutGroups" )    
        .def( bp::init< >() )    
        .def( bp::init< SireMol::MoleculeGroup const & >(( bp::arg("molgroup") )) )    
        .def( bp::init< SireMol::MoleculeGroup const &, SireVol::MappingFunction const & >(( bp::arg("molgroup"), bp::arg("mappingfunction") )) )    
        .def( bp::init< QList<SireMol::MoleculeGroup> const & >(( bp::arg("molgroups") )) )    
        .def( bp::init< QList<SireMol::MoleculeGroup> const &, SireVol::MappingFunction const & >(( bp::arg("molgroups"), bp::arg("mappingfunction") )) )    
        .def( bp::init< SireMol::MoleculeGroups const & >(( bp::arg("molgroups") )) )    
        .def( bp::init< SireMol::MoleculeGroups const &, SireVol::MappingFunction const & >(( bp::arg("molgroups"), bp::arg("mappingfunction") )) )    
        .def( 
            "assertCompatibleWith"
            , (void ( ::SireMove::MapAsCutGroups::* )( ::SireSystem::QuerySystem & ) const)( &::SireMove::MapAsCutGroups::assertCompatibleWith )
            , ( bp::arg("system") ) )    
        .def( 
            "map"
            , (::SireMol::Molecules ( ::SireMove::MapAsCutGroups::* )( ::SireSystem::QuerySystem &,::SireVol::Space const & ) const)( &::SireMove::MapAsCutGroups::map )
            , ( bp::arg("system"), bp::arg("new_space") ) )    
        .def( 
            "typeName"
            , (char const * (*)(  ))( &::SireMove::MapAsCutGroups::typeName ) )    
        .def( 
            "what"
            , (char const * ( ::SireMove::MapAsCutGroups::* )(  ) const)( &::SireMove::MapAsCutGroups::what ) )    
        .staticmethod( "typeName" )    
        .def( "__copy__", &__copy__)    
        .def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMove::MapAsCutGroups >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMove::MapAsCutGroups >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )    
        .def( "__str__", &pvt_get_name);

}
