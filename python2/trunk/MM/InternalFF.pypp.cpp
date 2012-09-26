// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "InternalFF.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireBase/property.h"

#include "SireError/errors.h"

#include "SireFF/detail/atomiccoords3d.h"

#include "SireFF/errors.h"

#include "SireMaths/line.h"

#include "SireMaths/torsion.h"

#include "SireMaths/triangle.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/units.h"

#include "internalff.h"

#include "tostring.h"

#include <QDebug>

#include <cstdio>

#include "internalff.h"

SireMM::InternalFF __copy__(const SireMM::InternalFF &other){ return SireMM::InternalFF(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_InternalFF_class(){

    { //::SireMM::InternalFF
        typedef bp::class_< SireMM::InternalFF, bp::bases< SireFF::FF3D, SireFF::G1FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > InternalFF_exposer_t;
        InternalFF_exposer_t InternalFF_exposer = InternalFF_exposer_t( "InternalFF", bp::init< >() );
        bp::scope InternalFF_scope( InternalFF_exposer );
        InternalFF_exposer.def( bp::init< QString const & >(( bp::arg("name") )) );
        InternalFF_exposer.def( bp::init< SireMM::InternalFF const & >(( bp::arg("other") )) );
        { //::SireMM::InternalFF::components
        
            typedef ::SireMM::InternalComponent const & ( ::SireMM::InternalFF::*components_function_type )(  ) const;
            components_function_type components_function_value( &::SireMM::InternalFF::components );
            
            InternalFF_exposer.def( 
                "components"
                , components_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::InternalFF::containsProperty
        
            typedef bool ( ::SireMM::InternalFF::*containsProperty_function_type )( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireMM::InternalFF::containsProperty );
            
            InternalFF_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") ) );
        
        }
        { //::SireMM::InternalFF::field
        
            typedef void ( ::SireMM::InternalFF::*field_function_type )( ::SireFF::FieldTable &,double ) ;
            field_function_type field_function_value( &::SireMM::InternalFF::field );
            
            InternalFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireMM::InternalFF::field
        
            typedef void ( ::SireMM::InternalFF::*field_function_type )( ::SireFF::FieldTable &,::SireCAS::Symbol const &,double ) ;
            field_function_type field_function_value( &::SireMM::InternalFF::field );
            
            InternalFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireMM::InternalFF::field
        
            typedef void ( ::SireMM::InternalFF::*field_function_type )( ::SireFF::FieldTable &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireMM::InternalFF::field );
            
            InternalFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireMM::InternalFF::field
        
            typedef void ( ::SireMM::InternalFF::*field_function_type )( ::SireFF::FieldTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireMM::InternalFF::field );
            
            InternalFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireMM::InternalFF::force
        
            typedef void ( ::SireMM::InternalFF::*force_function_type )( ::SireFF::ForceTable &,double ) ;
            force_function_type force_function_value( &::SireMM::InternalFF::force );
            
            InternalFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 ) );
        
        }
        { //::SireMM::InternalFF::force
        
            typedef void ( ::SireMM::InternalFF::*force_function_type )( ::SireFF::ForceTable &,::SireCAS::Symbol const &,double ) ;
            force_function_type force_function_value( &::SireMM::InternalFF::force );
            
            InternalFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("symbol"), bp::arg("scale_force")=1 ) );
        
        }
        { //::SireMM::InternalFF::isStrict
        
            typedef bool ( ::SireMM::InternalFF::*isStrict_function_type )(  ) const;
            isStrict_function_type isStrict_function_value( &::SireMM::InternalFF::isStrict );
            
            InternalFF_exposer.def( 
                "isStrict"
                , isStrict_function_value );
        
        }
        { //::SireMM::InternalFF::mustNowRecalculateFromScratch
        
            typedef void ( ::SireMM::InternalFF::*mustNowRecalculateFromScratch_function_type )(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireMM::InternalFF::mustNowRecalculateFromScratch );
            
            InternalFF_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value );
        
        }
        InternalFF_exposer.def( bp::self != bp::self );
        { //::SireMM::InternalFF::operator=
        
            typedef ::SireMM::InternalFF & ( ::SireMM::InternalFF::*assign_function_type )( ::SireMM::InternalFF const & ) ;
            assign_function_type assign_function_value( &::SireMM::InternalFF::operator= );
            
            InternalFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        InternalFF_exposer.def( bp::self == bp::self );
        { //::SireMM::InternalFF::parameters
        
            typedef ::SireMM::InternalParameterNames3D const & ( ::SireMM::InternalFF::*parameters_function_type )(  ) const;
            parameters_function_type parameters_function_value( &::SireMM::InternalFF::parameters );
            
            InternalFF_exposer.def( 
                "parameters"
                , parameters_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMM::InternalFF::potential
        
            typedef void ( ::SireMM::InternalFF::*potential_function_type )( ::SireFF::PotentialTable &,double ) ;
            potential_function_type potential_function_value( &::SireMM::InternalFF::potential );
            
            InternalFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireMM::InternalFF::potential
        
            typedef void ( ::SireMM::InternalFF::*potential_function_type )( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,double ) ;
            potential_function_type potential_function_value( &::SireMM::InternalFF::potential );
            
            InternalFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireMM::InternalFF::potential
        
            typedef void ( ::SireMM::InternalFF::*potential_function_type )( ::SireFF::PotentialTable &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireMM::InternalFF::potential );
            
            InternalFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("probe"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireMM::InternalFF::potential
        
            typedef void ( ::SireMM::InternalFF::*potential_function_type )( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireMM::InternalFF::potential );
            
            InternalFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireMM::InternalFF::properties
        
            typedef ::SireBase::Properties const & ( ::SireMM::InternalFF::*properties_function_type )(  ) const;
            properties_function_type properties_function_value( &::SireMM::InternalFF::properties );
            
            InternalFF_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMM::InternalFF::property
        
            typedef ::SireBase::Property const & ( ::SireMM::InternalFF::*property_function_type )( ::QString const & ) const;
            property_function_type property_function_value( &::SireMM::InternalFF::property );
            
            InternalFF_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::InternalFF::setProperty
        
            typedef bool ( ::SireMM::InternalFF::*setProperty_function_type )( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMM::InternalFF::setProperty );
            
            InternalFF_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("property") ) );
        
        }
        { //::SireMM::InternalFF::setStrict
        
            typedef bool ( ::SireMM::InternalFF::*setStrict_function_type )( bool ) ;
            setStrict_function_type setStrict_function_value( &::SireMM::InternalFF::setStrict );
            
            InternalFF_exposer.def( 
                "setStrict"
                , setStrict_function_value
                , ( bp::arg("isstrict") ) );
        
        }
        { //::SireMM::InternalFF::symbols
        
            typedef ::SireMM::InternalSymbols const & ( ::SireMM::InternalFF::*symbols_function_type )(  ) const;
            symbols_function_type symbols_function_value( &::SireMM::InternalFF::symbols );
            
            InternalFF_exposer.def( 
                "symbols"
                , symbols_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMM::InternalFF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::InternalFF::typeName );
            
            InternalFF_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        InternalFF_exposer.staticmethod( "typeName" );
        InternalFF_exposer.def( "__copy__", &__copy__);
        InternalFF_exposer.def( "__deepcopy__", &__copy__);
        InternalFF_exposer.def( "clone", &__copy__);
        InternalFF_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::InternalFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InternalFF_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::InternalFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InternalFF_exposer.def( "__str__", &__str__< ::SireMM::InternalFF > );
        InternalFF_exposer.def( "__repr__", &__str__< ::SireMM::InternalFF > );
    }

}
