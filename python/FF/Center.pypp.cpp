// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Center.pypp.hpp"

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

#include "forcetable.h"

#include "point.h"

#include <boost/tuple/tuple.hpp>

#include "point.h"

SireFF::Center __copy__(const SireFF::Center &other){ return SireFF::Center(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Center_class(){

    { //::SireFF::Center
        typedef bp::class_< SireFF::Center, bp::bases< SireFF::Point, SireBase::Property > > Center_exposer_t;
        Center_exposer_t Center_exposer = Center_exposer_t( "Center", bp::init< >() );
        bp::scope Center_scope( Center_exposer );
        Center_exposer.def( bp::init< SireMol::MoleculeView const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() )) );
        Center_exposer.def( bp::init< SireMol::Molecules const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecules"), bp::arg("map")=SireBase::PropertyMap() )) );
        Center_exposer.def( bp::init< SireFF::Center const & >(( bp::arg("other") )) );
        { //::SireFF::Center::addForce
        
            typedef bool ( ::SireFF::Center::*addForce_function_type )( ::SireFF::MolForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::Center::addForce );
            
            Center_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("molforces"), bp::arg("force") ) );
        
        }
        { //::SireFF::Center::addForce
        
            typedef bool ( ::SireFF::Center::*addForce_function_type )( ::SireFF::ForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::Center::addForce );
            
            Center_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("forces"), bp::arg("force") ) );
        
        }
        { //::SireFF::Center::contains
        
            typedef bool ( ::SireFF::Center::*contains_function_type )( ::SireMol::MolNum ) const;
            contains_function_type contains_function_value( &::SireFF::Center::contains );
            
            Center_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molnum") ) );
        
        }
        { //::SireFF::Center::contains
        
            typedef bool ( ::SireFF::Center::*contains_function_type )( ::SireMol::MolID const & ) const;
            contains_function_type contains_function_value( &::SireFF::Center::contains );
            
            Center_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molid") ) );
        
        }
        { //::SireFF::Center::isExtraMoleculePoint
        
            typedef bool ( ::SireFF::Center::*isExtraMoleculePoint_function_type )(  ) const;
            isExtraMoleculePoint_function_type isExtraMoleculePoint_function_value( &::SireFF::Center::isExtraMoleculePoint );
            
            Center_exposer.def( 
                "isExtraMoleculePoint"
                , isExtraMoleculePoint_function_value );
        
        }
        { //::SireFF::Center::isInterMoleculePoint
        
            typedef bool ( ::SireFF::Center::*isInterMoleculePoint_function_type )(  ) const;
            isInterMoleculePoint_function_type isInterMoleculePoint_function_value( &::SireFF::Center::isInterMoleculePoint );
            
            Center_exposer.def( 
                "isInterMoleculePoint"
                , isInterMoleculePoint_function_value );
        
        }
        { //::SireFF::Center::isIntraMoleculePoint
        
            typedef bool ( ::SireFF::Center::*isIntraMoleculePoint_function_type )(  ) const;
            isIntraMoleculePoint_function_type isIntraMoleculePoint_function_value( &::SireFF::Center::isIntraMoleculePoint );
            
            Center_exposer.def( 
                "isIntraMoleculePoint"
                , isIntraMoleculePoint_function_value );
        
        }
        { //::SireFF::Center::molecules
        
            typedef ::SireMol::Molecules ( ::SireFF::Center::*molecules_function_type )(  ) const;
            molecules_function_type molecules_function_value( &::SireFF::Center::molecules );
            
            Center_exposer.def( 
                "molecules"
                , molecules_function_value );
        
        }
        { //::SireFF::Center::nMolecules
        
            typedef int ( ::SireFF::Center::*nMolecules_function_type )(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireFF::Center::nMolecules );
            
            Center_exposer.def( 
                "nMolecules"
                , nMolecules_function_value );
        
        }
        Center_exposer.def( bp::self != bp::self );
        { //::SireFF::Center::operator=
        
            typedef ::SireFF::Center & ( ::SireFF::Center::*assign_function_type )( ::SireFF::Center const & ) ;
            assign_function_type assign_function_value( &::SireFF::Center::operator= );
            
            Center_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Center_exposer.def( bp::self == bp::self );
        { //::SireFF::Center::setSpace
        
            typedef void ( ::SireFF::Center::*setSpace_function_type )( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireFF::Center::setSpace );
            
            Center_exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("space") ) );
        
        }
        { //::SireFF::Center::toString
        
            typedef ::QString ( ::SireFF::Center::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireFF::Center::toString );
            
            Center_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireFF::Center::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Center::typeName );
            
            Center_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireFF::Center::update
        
            typedef bool ( ::SireFF::Center::*update_function_type )( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireFF::Center::update );
            
            Center_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") ) );
        
        }
        { //::SireFF::Center::update
        
            typedef bool ( ::SireFF::Center::*update_function_type )( ::SireMol::Molecules const & ) ;
            update_function_type update_function_value( &::SireFF::Center::update );
            
            Center_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molecules") ) );
        
        }
        { //::SireFF::Center::update
        
            typedef bool ( ::SireFF::Center::*update_function_type )( ::SireMol::MoleculeGroup const & ) ;
            update_function_type update_function_value( &::SireFF::Center::update );
            
            Center_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireFF::Center::update
        
            typedef bool ( ::SireFF::Center::*update_function_type )( ::SireMol::MolGroupsBase const & ) ;
            update_function_type update_function_value( &::SireFF::Center::update );
            
            Center_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroups") ) );
        
        }
        { //::SireFF::Center::usesMoleculesIn
        
            typedef bool ( ::SireFF::Center::*usesMoleculesIn_function_type )( ::SireFF::ForceTable const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Center::usesMoleculesIn );
            
            Center_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("forcetable") ) );
        
        }
        { //::SireFF::Center::usesMoleculesIn
        
            typedef bool ( ::SireFF::Center::*usesMoleculesIn_function_type )( ::SireMol::Molecules const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Center::usesMoleculesIn );
            
            Center_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molecules") ) );
        
        }
        { //::SireFF::Center::usesMoleculesIn
        
            typedef bool ( ::SireFF::Center::*usesMoleculesIn_function_type )( ::SireMol::MoleculeGroup const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Center::usesMoleculesIn );
            
            Center_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireFF::Center::usesMoleculesIn
        
            typedef bool ( ::SireFF::Center::*usesMoleculesIn_function_type )( ::SireMol::MolGroupsBase const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Center::usesMoleculesIn );
            
            Center_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroups") ) );
        
        }
        { //::SireFF::Center::wouldUpdate
        
            typedef bool ( ::SireFF::Center::*wouldUpdate_function_type )( ::SireMol::MoleculeData const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Center::wouldUpdate );
            
            Center_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("moldata") ) );
        
        }
        { //::SireFF::Center::wouldUpdate
        
            typedef bool ( ::SireFF::Center::*wouldUpdate_function_type )( ::SireMol::Molecules const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Center::wouldUpdate );
            
            Center_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molecules") ) );
        
        }
        { //::SireFF::Center::wouldUpdate
        
            typedef bool ( ::SireFF::Center::*wouldUpdate_function_type )( ::SireMol::MoleculeGroup const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Center::wouldUpdate );
            
            Center_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireFF::Center::wouldUpdate
        
            typedef bool ( ::SireFF::Center::*wouldUpdate_function_type )( ::SireMol::MolGroupsBase const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Center::wouldUpdate );
            
            Center_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroups") ) );
        
        }
        Center_exposer.staticmethod( "typeName" );
        Center_exposer.def( "__copy__", &__copy__);
        Center_exposer.def( "__deepcopy__", &__copy__);
        Center_exposer.def( "clone", &__copy__);
        Center_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::Center >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Center_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::Center >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Center_exposer.def( "__str__", &__str__< ::SireFF::Center > );
        Center_exposer.def( "__repr__", &__str__< ::SireFF::Center > );
    }

}
