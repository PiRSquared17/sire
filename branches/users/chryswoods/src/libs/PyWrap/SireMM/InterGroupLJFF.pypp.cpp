// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "InterGroupLJFF.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

SireMM::InterGroupLJFF __copy__(const SireMM::InterGroupLJFF &other){ return SireMM::InterGroupLJFF(other); }

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_InterGroupLJFF_class(){

    { //::SireMM::InterGroupLJFF
        typedef bp::class_< SireMM::InterGroupLJFF, bp::bases< SireMM::LJFF > > InterGroupLJFF_exposer_t;
        InterGroupLJFF_exposer_t InterGroupLJFF_exposer = InterGroupLJFF_exposer_t( "InterGroupLJFF" );
        bp::scope InterGroupLJFF_scope( InterGroupLJFF_exposer );
        bp::class_< SireMM::InterGroupLJFF::Components, bp::bases< SireMM::LJFF::Components > >( "Components" )    
            .def( bp::init< >() )    
            .def( bp::init< SireFF::FFBase const &, SireCAS::Symbols const & >(( bp::arg("ffbase"), bp::arg("symbols") )) );
        bp::class_< SireMM::InterGroupLJFF::Groups, bp::bases< SireMM::LJFF::Groups > >( "Groups" )    
            .def( bp::init< >() )    
            .def( 
                "A"
                , (::SireFF::FFBase::Group ( ::SireMM::InterGroupLJFF::Groups::* )(  ) const)( &::SireMM::InterGroupLJFF::Groups::A ) )    
            .def( 
                "B"
                , (::SireFF::FFBase::Group ( ::SireMM::InterGroupLJFF::Groups::* )(  ) const)( &::SireMM::InterGroupLJFF::Groups::B ) );
        bp::class_< SireMM::InterGroupLJFF::Parameters, bp::bases< SireMM::LJFF::Parameters > >( "Parameters" )    
            .def( bp::init< >() );
        InterGroupLJFF_exposer.def( bp::init< >() );
        InterGroupLJFF_exposer.def( bp::init< SireVol::Space const &, SireMM::SwitchingFunction const & >(( bp::arg("space"), bp::arg("switchfunc") )) );
        { //::SireMM::InterGroupLJFF::add
        
            typedef bool ( ::SireMM::InterGroupLJFF::*add_function_type )( ::SireMol::PartialMolecule const &,::SireFF::ParameterMap const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "add"
                , add_function_type( &::SireMM::InterGroupLJFF::add )
                , ( bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) );
        
        }
        { //::SireMM::InterGroupLJFF::add
        
            typedef bool ( ::SireMM::InterGroupLJFF::*add_function_type )( ::SireMol::Molecules const &,::SireFF::ParameterMap const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "add"
                , add_function_type( &::SireMM::InterGroupLJFF::add )
                , ( bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) );
        
        }
        { //::SireMM::InterGroupLJFF::addTo
        
            typedef bool ( ::SireMM::InterGroupLJFF::*addTo_function_type )( ::SireFF::FFBase::Group const &,::SireMol::PartialMolecule const &,::SireFF::ParameterMap const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "addTo"
                , addTo_function_type( &::SireMM::InterGroupLJFF::addTo )
                , ( bp::arg("group"), bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) );
        
        }
        { //::SireMM::InterGroupLJFF::addTo
        
            typedef bool ( ::SireMM::InterGroupLJFF::*addTo_function_type )( ::SireFF::FFBase::Group const &,::SireMol::Molecules const &,::SireFF::ParameterMap const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "addTo"
                , addTo_function_type( &::SireMM::InterGroupLJFF::addTo )
                , ( bp::arg("group"), bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) );
        
        }
        { //::SireMM::InterGroupLJFF::addToA
        
            typedef bool ( ::SireMM::InterGroupLJFF::*addToA_function_type )( ::SireMol::PartialMolecule const &,::SireFF::ParameterMap const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "addToA"
                , addToA_function_type( &::SireMM::InterGroupLJFF::addToA )
                , ( bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) );
        
        }
        { //::SireMM::InterGroupLJFF::addToA
        
            typedef bool ( ::SireMM::InterGroupLJFF::*addToA_function_type )( ::SireMol::Molecules const &,::SireFF::ParameterMap const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "addToA"
                , addToA_function_type( &::SireMM::InterGroupLJFF::addToA )
                , ( bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) );
        
        }
        { //::SireMM::InterGroupLJFF::addToB
        
            typedef bool ( ::SireMM::InterGroupLJFF::*addToB_function_type )( ::SireMol::PartialMolecule const &,::SireFF::ParameterMap const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "addToB"
                , addToB_function_type( &::SireMM::InterGroupLJFF::addToB )
                , ( bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) );
        
        }
        { //::SireMM::InterGroupLJFF::addToB
        
            typedef bool ( ::SireMM::InterGroupLJFF::*addToB_function_type )( ::SireMol::Molecules const &,::SireFF::ParameterMap const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "addToB"
                , addToB_function_type( &::SireMM::InterGroupLJFF::addToB )
                , ( bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) );
        
        }
        { //::SireMM::InterGroupLJFF::change
        
            typedef bool ( ::SireMM::InterGroupLJFF::*change_function_type )( ::SireMol::PartialMolecule const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "change"
                , change_function_type( &::SireMM::InterGroupLJFF::change )
                , ( bp::arg("molecule") ) );
        
        }
        { //::SireMM::InterGroupLJFF::change
        
            typedef bool ( ::SireMM::InterGroupLJFF::*change_function_type )( ::SireMol::Molecules const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "change"
                , change_function_type( &::SireMM::InterGroupLJFF::change )
                , ( bp::arg("molecules") ) );
        
        }
        { //::SireMM::InterGroupLJFF::contains
        
            typedef bool ( ::SireMM::InterGroupLJFF::*contains_function_type )( ::SireMol::PartialMolecule const & ) const;
            
            InterGroupLJFF_exposer.def( 
                "contains"
                , contains_function_type( &::SireMM::InterGroupLJFF::contains )
                , ( bp::arg("molecule") ) );
        
        }
        { //::SireMM::InterGroupLJFF::contains
        
            typedef bool ( ::SireMM::InterGroupLJFF::*contains_function_type )( ::SireMol::PartialMolecule const &,::SireFF::FFBase::Group const & ) const;
            
            InterGroupLJFF_exposer.def( 
                "contains"
                , contains_function_type( &::SireMM::InterGroupLJFF::contains )
                , ( bp::arg("molecule"), bp::arg("group") ) );
        
        }
        { //::SireMM::InterGroupLJFF::contents
        
            typedef ::SireMol::Molecules ( ::SireMM::InterGroupLJFF::*contents_function_type )(  ) const;
            
            InterGroupLJFF_exposer.def( 
                "contents"
                , contents_function_type( &::SireMM::InterGroupLJFF::contents ) );
        
        }
        { //::SireMM::InterGroupLJFF::contents
        
            typedef ::SireMol::Molecules ( ::SireMM::InterGroupLJFF::*contents_function_type )( ::SireFF::FFBase::Group const & ) const;
            
            InterGroupLJFF_exposer.def( 
                "contents"
                , contents_function_type( &::SireMM::InterGroupLJFF::contents )
                , ( bp::arg("group") ) );
        
        }
        { //::SireMM::InterGroupLJFF::groups
        
            typedef ::SireMM::InterGroupLJFF::Groups const & ( ::SireMM::InterGroupLJFF::*groups_function_type )(  ) const;
            
            InterGroupLJFF_exposer.def( 
                "groups"
                , groups_function_type( &::SireMM::InterGroupLJFF::groups )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMM::InterGroupLJFF::groupsReferringTo
        
            typedef ::QSet<SireFF::FFBase::Group> ( ::SireMM::InterGroupLJFF::*groupsReferringTo_function_type )( ::SireMol::MoleculeID ) const;
            
            InterGroupLJFF_exposer.def( 
                "groupsReferringTo"
                , groupsReferringTo_function_type( &::SireMM::InterGroupLJFF::groupsReferringTo )
                , ( bp::arg("molid") ) );
        
        }
        { //::SireMM::InterGroupLJFF::molecule
        
            typedef ::SireMol::PartialMolecule ( ::SireMM::InterGroupLJFF::*molecule_function_type )( ::SireMol::MoleculeID ) const;
            
            InterGroupLJFF_exposer.def( 
                "molecule"
                , molecule_function_type( &::SireMM::InterGroupLJFF::molecule )
                , ( bp::arg("molid") ) );
        
        }
        { //::SireMM::InterGroupLJFF::molecule
        
            typedef ::SireMol::PartialMolecule ( ::SireMM::InterGroupLJFF::*molecule_function_type )( ::SireMol::MoleculeID,::SireFF::FFBase::Group const & ) const;
            
            InterGroupLJFF_exposer.def( 
                "molecule"
                , molecule_function_type( &::SireMM::InterGroupLJFF::molecule )
                , ( bp::arg("molid"), bp::arg("group") ) );
        
        }
        { //::SireMM::InterGroupLJFF::moleculeIDs
        
            typedef ::QSet<SireMol::MoleculeID> ( ::SireMM::InterGroupLJFF::*moleculeIDs_function_type )(  ) const;
            
            InterGroupLJFF_exposer.def( 
                "moleculeIDs"
                , moleculeIDs_function_type( &::SireMM::InterGroupLJFF::moleculeIDs ) );
        
        }
        { //::SireMM::InterGroupLJFF::moleculeIDs
        
            typedef ::QSet<SireMol::MoleculeID> ( ::SireMM::InterGroupLJFF::*moleculeIDs_function_type )( ::SireFF::FFBase::Group const & ) const;
            
            InterGroupLJFF_exposer.def( 
                "moleculeIDs"
                , moleculeIDs_function_type( &::SireMM::InterGroupLJFF::moleculeIDs )
                , ( bp::arg("group") ) );
        
        }
        { //::SireMM::InterGroupLJFF::mustNowRecalculateFromScratch
        
            typedef void ( ::SireMM::InterGroupLJFF::*mustNowRecalculateFromScratch_function_type )(  ) ;
            
            InterGroupLJFF_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_type( &::SireMM::InterGroupLJFF::mustNowRecalculateFromScratch ) );
        
        }
        { //::SireMM::InterGroupLJFF::refersTo
        
            typedef bool ( ::SireMM::InterGroupLJFF::*refersTo_function_type )( ::SireMol::MoleculeID ) const;
            
            InterGroupLJFF_exposer.def( 
                "refersTo"
                , refersTo_function_type( &::SireMM::InterGroupLJFF::refersTo )
                , ( bp::arg("molid") ) );
        
        }
        { //::SireMM::InterGroupLJFF::refersTo
        
            typedef bool ( ::SireMM::InterGroupLJFF::*refersTo_function_type )( ::SireMol::MoleculeID,::SireFF::FFBase::Group const & ) const;
            
            InterGroupLJFF_exposer.def( 
                "refersTo"
                , refersTo_function_type( &::SireMM::InterGroupLJFF::refersTo )
                , ( bp::arg("molid"), bp::arg("group") ) );
        
        }
        { //::SireMM::InterGroupLJFF::remove
        
            typedef bool ( ::SireMM::InterGroupLJFF::*remove_function_type )( ::SireMol::PartialMolecule const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "remove"
                , remove_function_type( &::SireMM::InterGroupLJFF::remove )
                , ( bp::arg("molecule") ) );
        
        }
        { //::SireMM::InterGroupLJFF::remove
        
            typedef bool ( ::SireMM::InterGroupLJFF::*remove_function_type )( ::SireMol::Molecules const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "remove"
                , remove_function_type( &::SireMM::InterGroupLJFF::remove )
                , ( bp::arg("molecules") ) );
        
        }
        { //::SireMM::InterGroupLJFF::removeFrom
        
            typedef bool ( ::SireMM::InterGroupLJFF::*removeFrom_function_type )( ::SireFF::FFBase::Group const &,::SireMol::PartialMolecule const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "removeFrom"
                , removeFrom_function_type( &::SireMM::InterGroupLJFF::removeFrom )
                , ( bp::arg("group"), bp::arg("molecule") ) );
        
        }
        { //::SireMM::InterGroupLJFF::removeFrom
        
            typedef bool ( ::SireMM::InterGroupLJFF::*removeFrom_function_type )( ::SireFF::FFBase::Group const &,::SireMol::Molecules const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "removeFrom"
                , removeFrom_function_type( &::SireMM::InterGroupLJFF::removeFrom )
                , ( bp::arg("group"), bp::arg("molecules") ) );
        
        }
        { //::SireMM::InterGroupLJFF::removeFromA
        
            typedef bool ( ::SireMM::InterGroupLJFF::*removeFromA_function_type )( ::SireMol::PartialMolecule const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "removeFromA"
                , removeFromA_function_type( &::SireMM::InterGroupLJFF::removeFromA )
                , ( bp::arg("molecule") ) );
        
        }
        { //::SireMM::InterGroupLJFF::removeFromA
        
            typedef bool ( ::SireMM::InterGroupLJFF::*removeFromA_function_type )( ::SireMol::Molecules const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "removeFromA"
                , removeFromA_function_type( &::SireMM::InterGroupLJFF::removeFromA )
                , ( bp::arg("molecules") ) );
        
        }
        { //::SireMM::InterGroupLJFF::removeFromB
        
            typedef bool ( ::SireMM::InterGroupLJFF::*removeFromB_function_type )( ::SireMol::PartialMolecule const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "removeFromB"
                , removeFromB_function_type( &::SireMM::InterGroupLJFF::removeFromB )
                , ( bp::arg("molecule") ) );
        
        }
        { //::SireMM::InterGroupLJFF::removeFromB
        
            typedef bool ( ::SireMM::InterGroupLJFF::*removeFromB_function_type )( ::SireMol::Molecules const & ) ;
            
            InterGroupLJFF_exposer.def( 
                "removeFromB"
                , removeFromB_function_type( &::SireMM::InterGroupLJFF::removeFromB )
                , ( bp::arg("molecules") ) );
        
        }
        { //::SireMM::InterGroupLJFF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            
            InterGroupLJFF_exposer.def( 
                "typeName"
                , typeName_function_type( &::SireMM::InterGroupLJFF::typeName ) );
        
        }
        { //::SireMM::InterGroupLJFF::what
        
            typedef char const * ( ::SireMM::InterGroupLJFF::*what_function_type )(  ) const;
            
            InterGroupLJFF_exposer.def( 
                "what"
                , what_function_type( &::SireMM::InterGroupLJFF::what ) );
        
        }
        InterGroupLJFF_exposer.staticmethod( "typeName" );
        InterGroupLJFF_exposer.def( "__copy__", &__copy__);
        InterGroupLJFF_exposer.def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMM::InterGroupLJFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InterGroupLJFF_exposer.def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMM::InterGroupLJFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InterGroupLJFF_exposer.def( "__str__", &SirePy::__str__< ::SireMM::InterGroupLJFF > );
    }

}
