// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "RestraintFF.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireCAS/errors.h"

#include "SireError/errors.h"

#include "SireFF/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculedata.h"

#include "SireMol/molecules.h"

#include "SireMol/partialmolecule.h"

#include "SireMol/viewsofmol.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "restraintff.h"

#include "restraintff.h"

SireMM::RestraintFF __copy__(const SireMM::RestraintFF &other){ return SireMM::RestraintFF(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_RestraintFF_class(){

    { //::SireMM::RestraintFF
        typedef bp::class_< SireMM::RestraintFF, bp::bases< SireFF::FF3D, SireFF::G1FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > RestraintFF_exposer_t;
        RestraintFF_exposer_t RestraintFF_exposer = RestraintFF_exposer_t( "RestraintFF", bp::init< >() );
        bp::scope RestraintFF_scope( RestraintFF_exposer );
        RestraintFF_exposer.def( bp::init< QString const & >(( bp::arg("name") )) );
        RestraintFF_exposer.def( bp::init< SireMM::RestraintFF const & >(( bp::arg("other") )) );
        { //::SireMM::RestraintFF::add
        
            typedef bool ( ::SireMM::RestraintFF::*add_function_type )( ::SireMM::Restraint3D const & ) ;
            add_function_type add_function_value( &::SireMM::RestraintFF::add );
            
            RestraintFF_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("restraint") ) );
        
        }
        { //::SireMM::RestraintFF::builtinSymbols
        
            typedef ::SireCAS::Symbols ( ::SireMM::RestraintFF::*builtinSymbols_function_type )(  ) const;
            builtinSymbols_function_type builtinSymbols_function_value( &::SireMM::RestraintFF::builtinSymbols );
            
            RestraintFF_exposer.def( 
                "builtinSymbols"
                , builtinSymbols_function_value );
        
        }
        { //::SireMM::RestraintFF::components
        
            typedef ::SireMM::RestraintComponent const & ( ::SireMM::RestraintFF::*components_function_type )(  ) const;
            components_function_type components_function_value( &::SireMM::RestraintFF::components );
            
            RestraintFF_exposer.def( 
                "components"
                , components_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::RestraintFF::contains
        
            typedef bool ( ::SireMM::RestraintFF::*contains_function_type )( ::SireMM::Restraint3D const & ) const;
            contains_function_type contains_function_value( &::SireMM::RestraintFF::contains );
            
            RestraintFF_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("restraint") ) );
        
        }
        { //::SireMM::RestraintFF::containsProperty
        
            typedef bool ( ::SireMM::RestraintFF::*containsProperty_function_type )( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireMM::RestraintFF::containsProperty );
            
            RestraintFF_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") ) );
        
        }
        { //::SireMM::RestraintFF::differentiate
        
            typedef ::SireMM::RestraintFF ( ::SireMM::RestraintFF::*differentiate_function_type )( ::SireCAS::Symbol const & ) const;
            differentiate_function_type differentiate_function_value( &::SireMM::RestraintFF::differentiate );
            
            RestraintFF_exposer.def( 
                "differentiate"
                , differentiate_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireMM::RestraintFF::field
        
            typedef void ( ::SireMM::RestraintFF::*field_function_type )( ::SireFF::FieldTable &,double ) ;
            field_function_type field_function_value( &::SireMM::RestraintFF::field );
            
            RestraintFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireMM::RestraintFF::field
        
            typedef void ( ::SireMM::RestraintFF::*field_function_type )( ::SireFF::FieldTable &,::SireCAS::Symbol const &,double ) ;
            field_function_type field_function_value( &::SireMM::RestraintFF::field );
            
            RestraintFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireMM::RestraintFF::field
        
            typedef void ( ::SireMM::RestraintFF::*field_function_type )( ::SireFF::FieldTable &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireMM::RestraintFF::field );
            
            RestraintFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireMM::RestraintFF::field
        
            typedef void ( ::SireMM::RestraintFF::*field_function_type )( ::SireFF::FieldTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireMM::RestraintFF::field );
            
            RestraintFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireMM::RestraintFF::force
        
            typedef void ( ::SireMM::RestraintFF::*force_function_type )( ::SireFF::ForceTable &,double ) ;
            force_function_type force_function_value( &::SireMM::RestraintFF::force );
            
            RestraintFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 ) );
        
        }
        { //::SireMM::RestraintFF::force
        
            typedef void ( ::SireMM::RestraintFF::*force_function_type )( ::SireFF::ForceTable &,::SireCAS::Symbol const &,double ) ;
            force_function_type force_function_value( &::SireMM::RestraintFF::force );
            
            RestraintFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("symbol"), bp::arg("scale_force")=1 ) );
        
        }
        { //::SireMM::RestraintFF::getValue
        
            typedef double ( ::SireMM::RestraintFF::*getValue_function_type )( ::SireCAS::Symbol const & ) const;
            getValue_function_type getValue_function_value( &::SireMM::RestraintFF::getValue );
            
            RestraintFF_exposer.def( 
                "getValue"
                , getValue_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireMM::RestraintFF::hasValue
        
            typedef bool ( ::SireMM::RestraintFF::*hasValue_function_type )( ::SireCAS::Symbol const & ) const;
            hasValue_function_type hasValue_function_value( &::SireMM::RestraintFF::hasValue );
            
            RestraintFF_exposer.def( 
                "hasValue"
                , hasValue_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireMM::RestraintFF::mustNowRecalculateFromScratch
        
            typedef void ( ::SireMM::RestraintFF::*mustNowRecalculateFromScratch_function_type )(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireMM::RestraintFF::mustNowRecalculateFromScratch );
            
            RestraintFF_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value );
        
        }
        { //::SireMM::RestraintFF::nRestraints
        
            typedef int ( ::SireMM::RestraintFF::*nRestraints_function_type )(  ) const;
            nRestraints_function_type nRestraints_function_value( &::SireMM::RestraintFF::nRestraints );
            
            RestraintFF_exposer.def( 
                "nRestraints"
                , nRestraints_function_value );
        
        }
        RestraintFF_exposer.def( bp::self != bp::self );
        { //::SireMM::RestraintFF::operator=
        
            typedef ::SireMM::RestraintFF & ( ::SireMM::RestraintFF::*assign_function_type )( ::SireMM::RestraintFF const & ) ;
            assign_function_type assign_function_value( &::SireMM::RestraintFF::operator= );
            
            RestraintFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        RestraintFF_exposer.def( bp::self == bp::self );
        { //::SireMM::RestraintFF::potential
        
            typedef void ( ::SireMM::RestraintFF::*potential_function_type )( ::SireFF::PotentialTable &,double ) ;
            potential_function_type potential_function_value( &::SireMM::RestraintFF::potential );
            
            RestraintFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireMM::RestraintFF::potential
        
            typedef void ( ::SireMM::RestraintFF::*potential_function_type )( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,double ) ;
            potential_function_type potential_function_value( &::SireMM::RestraintFF::potential );
            
            RestraintFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireMM::RestraintFF::potential
        
            typedef void ( ::SireMM::RestraintFF::*potential_function_type )( ::SireFF::PotentialTable &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireMM::RestraintFF::potential );
            
            RestraintFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("probe"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireMM::RestraintFF::potential
        
            typedef void ( ::SireMM::RestraintFF::*potential_function_type )( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireMM::RestraintFF::potential );
            
            RestraintFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireMM::RestraintFF::properties
        
            typedef ::SireBase::Properties const & ( ::SireMM::RestraintFF::*properties_function_type )(  ) const;
            properties_function_type properties_function_value( &::SireMM::RestraintFF::properties );
            
            RestraintFF_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMM::RestraintFF::property
        
            typedef ::SireBase::Property const & ( ::SireMM::RestraintFF::*property_function_type )( ::QString const & ) const;
            property_function_type property_function_value( &::SireMM::RestraintFF::property );
            
            RestraintFF_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::RestraintFF::remove
        
            typedef bool ( ::SireMM::RestraintFF::*remove_function_type )( ::SireMM::Restraint3D const & ) ;
            remove_function_type remove_function_value( &::SireMM::RestraintFF::remove );
            
            RestraintFF_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("restraint") ) );
        
        }
        { //::SireMM::RestraintFF::removeRestraintAt
        
            typedef void ( ::SireMM::RestraintFF::*removeRestraintAt_function_type )( int ) ;
            removeRestraintAt_function_type removeRestraintAt_function_value( &::SireMM::RestraintFF::removeRestraintAt );
            
            RestraintFF_exposer.def( 
                "removeRestraintAt"
                , removeRestraintAt_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMM::RestraintFF::restraintAt
        
            typedef ::SireMM::Restraint3D const & ( ::SireMM::RestraintFF::*restraintAt_function_type )( int ) const;
            restraintAt_function_type restraintAt_function_value( &::SireMM::RestraintFF::restraintAt );
            
            RestraintFF_exposer.def( 
                "restraintAt"
                , restraintAt_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::RestraintFF::restraints
        
            typedef ::QVector< SireBase::PropPtr< SireMM::Restraint3D > > ( ::SireMM::RestraintFF::*restraints_function_type )(  ) const;
            restraints_function_type restraints_function_value( &::SireMM::RestraintFF::restraints );
            
            RestraintFF_exposer.def( 
                "restraints"
                , restraints_function_value );
        
        }
        { //::SireMM::RestraintFF::setProperty
        
            typedef bool ( ::SireMM::RestraintFF::*setProperty_function_type )( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMM::RestraintFF::setProperty );
            
            RestraintFF_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("property") ) );
        
        }
        { //::SireMM::RestraintFF::setSpace
        
            typedef bool ( ::SireMM::RestraintFF::*setSpace_function_type )( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireMM::RestraintFF::setSpace );
            
            RestraintFF_exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("space") ) );
        
        }
        { //::SireMM::RestraintFF::setValue
        
            typedef bool ( ::SireMM::RestraintFF::*setValue_function_type )( ::SireCAS::Symbol const &,double ) ;
            setValue_function_type setValue_function_value( &::SireMM::RestraintFF::setValue );
            
            RestraintFF_exposer.def( 
                "setValue"
                , setValue_function_value
                , ( bp::arg("symbol"), bp::arg("value") ) );
        
        }
        { //::SireMM::RestraintFF::space
        
            typedef ::SireVol::Space const & ( ::SireMM::RestraintFF::*space_function_type )(  ) const;
            space_function_type space_function_value( &::SireMM::RestraintFF::space );
            
            RestraintFF_exposer.def( 
                "space"
                , space_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::RestraintFF::symbols
        
            typedef ::SireCAS::Symbols ( ::SireMM::RestraintFF::*symbols_function_type )(  ) const;
            symbols_function_type symbols_function_value( &::SireMM::RestraintFF::symbols );
            
            RestraintFF_exposer.def( 
                "symbols"
                , symbols_function_value );
        
        }
        { //::SireMM::RestraintFF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::RestraintFF::typeName );
            
            RestraintFF_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMM::RestraintFF::userSymbols
        
            typedef ::SireCAS::Symbols ( ::SireMM::RestraintFF::*userSymbols_function_type )(  ) const;
            userSymbols_function_type userSymbols_function_value( &::SireMM::RestraintFF::userSymbols );
            
            RestraintFF_exposer.def( 
                "userSymbols"
                , userSymbols_function_value );
        
        }
        { //::SireMM::RestraintFF::userValues
        
            typedef ::SireCAS::Values ( ::SireMM::RestraintFF::*userValues_function_type )(  ) const;
            userValues_function_type userValues_function_value( &::SireMM::RestraintFF::userValues );
            
            RestraintFF_exposer.def( 
                "userValues"
                , userValues_function_value );
        
        }
        RestraintFF_exposer.staticmethod( "typeName" );
        RestraintFF_exposer.def( "__copy__", &__copy__);
        RestraintFF_exposer.def( "__deepcopy__", &__copy__);
        RestraintFF_exposer.def( "clone", &__copy__);
        RestraintFF_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::RestraintFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RestraintFF_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::RestraintFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RestraintFF_exposer.def( "__str__", &__str__< ::SireMM::RestraintFF > );
        RestraintFF_exposer.def( "__repr__", &__str__< ::SireMM::RestraintFF > );
    }

}
