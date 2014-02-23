// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CenterOfMass.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/evaluator.h"

#include "SireMol/mgidx.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleculegroups.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/aabox.h"

#include "SireVol/errors.h"

#include "energytable.h"

#include "forcetable.h"

#include "point.h"

#include <boost/tuple/tuple.hpp>

#include "point.h"

SireFF::CenterOfMass __copy__(const SireFF::CenterOfMass &other){ return SireFF::CenterOfMass(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_CenterOfMass_class(){

    { //::SireFF::CenterOfMass
        typedef bp::class_< SireFF::CenterOfMass, bp::bases< SireFF::Point, SireBase::Property > > CenterOfMass_exposer_t;
        CenterOfMass_exposer_t CenterOfMass_exposer = CenterOfMass_exposer_t( "CenterOfMass" );
        bp::scope CenterOfMass_scope( CenterOfMass_exposer );
        CenterOfMass_exposer.def( bp::init< >() );
        CenterOfMass_exposer.def( bp::init< SireMol::MoleculeView const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() )) );
        CenterOfMass_exposer.def( bp::init< SireMol::Molecules const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecules"), bp::arg("map")=SireBase::PropertyMap() )) );
        CenterOfMass_exposer.def( bp::init< SireFF::CenterOfMass const & >(( bp::arg("other") )) );
        { //::SireFF::CenterOfMass::addEnergy
        
            typedef bool ( ::SireFF::CenterOfMass::*addEnergy_function_type )( ::SireFF::MolEnergyTable &,double const & ) const;
            addEnergy_function_type addEnergy_function_value( &::SireFF::CenterOfMass::addEnergy );
            
            CenterOfMass_exposer.def( 
                "addEnergy"
                , addEnergy_function_value
                , ( bp::arg("molenergies"), bp::arg("energy") ) );
        
        }
        { //::SireFF::CenterOfMass::addEnergy
        
            typedef bool ( ::SireFF::CenterOfMass::*addEnergy_function_type )( ::SireFF::EnergyTable &,double const & ) const;
            addEnergy_function_type addEnergy_function_value( &::SireFF::CenterOfMass::addEnergy );
            
            CenterOfMass_exposer.def( 
                "addEnergy"
                , addEnergy_function_value
                , ( bp::arg("energies"), bp::arg("energy") ) );
        
        }
        { //::SireFF::CenterOfMass::addForce
        
            typedef bool ( ::SireFF::CenterOfMass::*addForce_function_type )( ::SireFF::MolForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::CenterOfMass::addForce );
            
            CenterOfMass_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("molforces"), bp::arg("force") ) );
        
        }
        { //::SireFF::CenterOfMass::addForce
        
            typedef bool ( ::SireFF::CenterOfMass::*addForce_function_type )( ::SireFF::ForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::CenterOfMass::addForce );
            
            CenterOfMass_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("forces"), bp::arg("force") ) );
        
        }
        { //::SireFF::CenterOfMass::contains
        
            typedef bool ( ::SireFF::CenterOfMass::*contains_function_type )( ::SireMol::MolNum ) const;
            contains_function_type contains_function_value( &::SireFF::CenterOfMass::contains );
            
            CenterOfMass_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molnum") ) );
        
        }
        { //::SireFF::CenterOfMass::contains
        
            typedef bool ( ::SireFF::CenterOfMass::*contains_function_type )( ::SireMol::MolID const & ) const;
            contains_function_type contains_function_value( &::SireFF::CenterOfMass::contains );
            
            CenterOfMass_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molid") ) );
        
        }
        { //::SireFF::CenterOfMass::isExtraMoleculePoint
        
            typedef bool ( ::SireFF::CenterOfMass::*isExtraMoleculePoint_function_type )(  ) const;
            isExtraMoleculePoint_function_type isExtraMoleculePoint_function_value( &::SireFF::CenterOfMass::isExtraMoleculePoint );
            
            CenterOfMass_exposer.def( 
                "isExtraMoleculePoint"
                , isExtraMoleculePoint_function_value );
        
        }
        { //::SireFF::CenterOfMass::isInterMoleculePoint
        
            typedef bool ( ::SireFF::CenterOfMass::*isInterMoleculePoint_function_type )(  ) const;
            isInterMoleculePoint_function_type isInterMoleculePoint_function_value( &::SireFF::CenterOfMass::isInterMoleculePoint );
            
            CenterOfMass_exposer.def( 
                "isInterMoleculePoint"
                , isInterMoleculePoint_function_value );
        
        }
        { //::SireFF::CenterOfMass::isIntraMoleculePoint
        
            typedef bool ( ::SireFF::CenterOfMass::*isIntraMoleculePoint_function_type )(  ) const;
            isIntraMoleculePoint_function_type isIntraMoleculePoint_function_value( &::SireFF::CenterOfMass::isIntraMoleculePoint );
            
            CenterOfMass_exposer.def( 
                "isIntraMoleculePoint"
                , isIntraMoleculePoint_function_value );
        
        }
        { //::SireFF::CenterOfMass::molecules
        
            typedef ::SireMol::Molecules ( ::SireFF::CenterOfMass::*molecules_function_type )(  ) const;
            molecules_function_type molecules_function_value( &::SireFF::CenterOfMass::molecules );
            
            CenterOfMass_exposer.def( 
                "molecules"
                , molecules_function_value );
        
        }
        { //::SireFF::CenterOfMass::nMolecules
        
            typedef int ( ::SireFF::CenterOfMass::*nMolecules_function_type )(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireFF::CenterOfMass::nMolecules );
            
            CenterOfMass_exposer.def( 
                "nMolecules"
                , nMolecules_function_value );
        
        }
        CenterOfMass_exposer.def( bp::self != bp::self );
        { //::SireFF::CenterOfMass::operator=
        
            typedef ::SireFF::CenterOfMass & ( ::SireFF::CenterOfMass::*assign_function_type )( ::SireFF::CenterOfMass const & ) ;
            assign_function_type assign_function_value( &::SireFF::CenterOfMass::operator= );
            
            CenterOfMass_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        CenterOfMass_exposer.def( bp::self == bp::self );
        { //::SireFF::CenterOfMass::setSpace
        
            typedef void ( ::SireFF::CenterOfMass::*setSpace_function_type )( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireFF::CenterOfMass::setSpace );
            
            CenterOfMass_exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("space") ) );
        
        }
        { //::SireFF::CenterOfMass::toString
        
            typedef ::QString ( ::SireFF::CenterOfMass::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireFF::CenterOfMass::toString );
            
            CenterOfMass_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireFF::CenterOfMass::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::CenterOfMass::typeName );
            
            CenterOfMass_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireFF::CenterOfMass::update
        
            typedef bool ( ::SireFF::CenterOfMass::*update_function_type )( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireFF::CenterOfMass::update );
            
            CenterOfMass_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") ) );
        
        }
        { //::SireFF::CenterOfMass::update
        
            typedef bool ( ::SireFF::CenterOfMass::*update_function_type )( ::SireMol::Molecules const & ) ;
            update_function_type update_function_value( &::SireFF::CenterOfMass::update );
            
            CenterOfMass_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molecules") ) );
        
        }
        { //::SireFF::CenterOfMass::update
        
            typedef bool ( ::SireFF::CenterOfMass::*update_function_type )( ::SireMol::MoleculeGroup const & ) ;
            update_function_type update_function_value( &::SireFF::CenterOfMass::update );
            
            CenterOfMass_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireFF::CenterOfMass::update
        
            typedef bool ( ::SireFF::CenterOfMass::*update_function_type )( ::SireMol::MolGroupsBase const & ) ;
            update_function_type update_function_value( &::SireFF::CenterOfMass::update );
            
            CenterOfMass_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroups") ) );
        
        }
        { //::SireFF::CenterOfMass::usesMoleculesIn
        
            typedef bool ( ::SireFF::CenterOfMass::*usesMoleculesIn_function_type )( ::SireFF::ForceTable const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::CenterOfMass::usesMoleculesIn );
            
            CenterOfMass_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("forcetable") ) );
        
        }
        { //::SireFF::CenterOfMass::usesMoleculesIn
        
            typedef bool ( ::SireFF::CenterOfMass::*usesMoleculesIn_function_type )( ::SireFF::EnergyTable const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::CenterOfMass::usesMoleculesIn );
            
            CenterOfMass_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("energytable") ) );
        
        }
        { //::SireFF::CenterOfMass::usesMoleculesIn
        
            typedef bool ( ::SireFF::CenterOfMass::*usesMoleculesIn_function_type )( ::SireMol::Molecules const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::CenterOfMass::usesMoleculesIn );
            
            CenterOfMass_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molecules") ) );
        
        }
        { //::SireFF::CenterOfMass::usesMoleculesIn
        
            typedef bool ( ::SireFF::CenterOfMass::*usesMoleculesIn_function_type )( ::SireMol::MoleculeGroup const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::CenterOfMass::usesMoleculesIn );
            
            CenterOfMass_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireFF::CenterOfMass::usesMoleculesIn
        
            typedef bool ( ::SireFF::CenterOfMass::*usesMoleculesIn_function_type )( ::SireMol::MolGroupsBase const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::CenterOfMass::usesMoleculesIn );
            
            CenterOfMass_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroups") ) );
        
        }
        { //::SireFF::CenterOfMass::wouldUpdate
        
            typedef bool ( ::SireFF::CenterOfMass::*wouldUpdate_function_type )( ::SireMol::MoleculeData const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::CenterOfMass::wouldUpdate );
            
            CenterOfMass_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("moldata") ) );
        
        }
        { //::SireFF::CenterOfMass::wouldUpdate
        
            typedef bool ( ::SireFF::CenterOfMass::*wouldUpdate_function_type )( ::SireMol::Molecules const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::CenterOfMass::wouldUpdate );
            
            CenterOfMass_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molecules") ) );
        
        }
        { //::SireFF::CenterOfMass::wouldUpdate
        
            typedef bool ( ::SireFF::CenterOfMass::*wouldUpdate_function_type )( ::SireMol::MoleculeGroup const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::CenterOfMass::wouldUpdate );
            
            CenterOfMass_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireFF::CenterOfMass::wouldUpdate
        
            typedef bool ( ::SireFF::CenterOfMass::*wouldUpdate_function_type )( ::SireMol::MolGroupsBase const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::CenterOfMass::wouldUpdate );
            
            CenterOfMass_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroups") ) );
        
        }
        CenterOfMass_exposer.staticmethod( "typeName" );
        CenterOfMass_exposer.def( "__copy__", &__copy__);
        CenterOfMass_exposer.def( "__deepcopy__", &__copy__);
        CenterOfMass_exposer.def( "clone", &__copy__);
        CenterOfMass_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::CenterOfMass >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CenterOfMass_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::CenterOfMass >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CenterOfMass_exposer.def( "__str__", &__str__< ::SireFF::CenterOfMass > );
        CenterOfMass_exposer.def( "__repr__", &__str__< ::SireFF::CenterOfMass > );
    }

}
