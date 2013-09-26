// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "PartialMolecule.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "chain.h"

#include "cutgroup.h"

#include "editor.hpp"

#include "evaluator.h"

#include "molecule.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include <QDebug>

#include "partialmolecule.h"

SireMol::PartialMolecule __copy__(const SireMol::PartialMolecule &other){ return SireMol::PartialMolecule(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_PartialMolecule_class(){

    { //::SireMol::PartialMolecule
        typedef bp::class_< SireMol::PartialMolecule, bp::bases< SireMol::MoleculeView, SireBase::Property > > PartialMolecule_exposer_t;
        PartialMolecule_exposer_t PartialMolecule_exposer = PartialMolecule_exposer_t( "PartialMolecule", bp::init< >() );
        bp::scope PartialMolecule_scope( PartialMolecule_exposer );
        PartialMolecule_exposer.def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molecule") )) );
        PartialMolecule_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomSelection const & >(( bp::arg("moldata"), bp::arg("atoms") )) );
        PartialMolecule_exposer.def( bp::init< SireMol::PartialMolecule const & >(( bp::arg("other") )) );
        { //::SireMol::PartialMolecule::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMol::PartialMolecule::*evaluate_function_type )(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::PartialMolecule::evaluate );
            
            PartialMolecule_exposer.def( 
                "evaluate"
                , evaluate_function_value );
        
        }
        { //::SireMol::PartialMolecule::extract
        
            typedef ::SireMol::PartialMolecule ( ::SireMol::PartialMolecule::*extract_function_type )(  ) const;
            extract_function_type extract_function_value( &::SireMol::PartialMolecule::extract );
            
            PartialMolecule_exposer.def( 
                "extract"
                , extract_function_value );
        
        }
        { //::SireMol::PartialMolecule::hasMetadata
        
            typedef bool ( ::SireMol::PartialMolecule::*hasMetadata_function_type )( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::PartialMolecule::hasMetadata );
            
            PartialMolecule_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") ) );
        
        }
        { //::SireMol::PartialMolecule::hasMetadata
        
            typedef bool ( ::SireMol::PartialMolecule::*hasMetadata_function_type )( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::PartialMolecule::hasMetadata );
            
            PartialMolecule_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") ) );
        
        }
        { //::SireMol::PartialMolecule::hasProperty
        
            typedef bool ( ::SireMol::PartialMolecule::*hasProperty_function_type )( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::PartialMolecule::hasProperty );
            
            PartialMolecule_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::PartialMolecule::isEmpty
        
            typedef bool ( ::SireMol::PartialMolecule::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::PartialMolecule::isEmpty );
            
            PartialMolecule_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::PartialMolecule::metadata
        
            typedef ::SireBase::Property const & ( ::SireMol::PartialMolecule::*metadata_function_type )( ::SireBase::PropertyName const & ) const;
            metadata_function_type metadata_function_value( &::SireMol::PartialMolecule::metadata );
            
            PartialMolecule_exposer.def( 
                "metadata"
                , metadata_function_value
                , ( bp::arg("metakey") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMol::PartialMolecule::metadata
        
            typedef ::SireBase::Property const & ( ::SireMol::PartialMolecule::*metadata_function_type )( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            metadata_function_type metadata_function_value( &::SireMol::PartialMolecule::metadata );
            
            PartialMolecule_exposer.def( 
                "metadata"
                , metadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMol::PartialMolecule::metadataKeys
        
            typedef ::QStringList ( ::SireMol::PartialMolecule::*metadataKeys_function_type )(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::PartialMolecule::metadataKeys );
            
            PartialMolecule_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value );
        
        }
        { //::SireMol::PartialMolecule::metadataKeys
        
            typedef ::QStringList ( ::SireMol::PartialMolecule::*metadataKeys_function_type )( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::PartialMolecule::metadataKeys );
            
            PartialMolecule_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::PartialMolecule::move
        
            typedef ::SireMol::Mover< SireMol::PartialMolecule > ( ::SireMol::PartialMolecule::*move_function_type )(  ) const;
            move_function_type move_function_value( &::SireMol::PartialMolecule::move );
            
            PartialMolecule_exposer.def( 
                "move"
                , move_function_value );
        
        }
        { //::SireMol::PartialMolecule::nAtoms
        
            typedef int ( ::SireMol::PartialMolecule::*nAtoms_function_type )(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::PartialMolecule::nAtoms );
            
            PartialMolecule_exposer.def( 
                "nAtoms"
                , nAtoms_function_value );
        
        }
        { //::SireMol::PartialMolecule::nChains
        
            typedef int ( ::SireMol::PartialMolecule::*nChains_function_type )(  ) const;
            nChains_function_type nChains_function_value( &::SireMol::PartialMolecule::nChains );
            
            PartialMolecule_exposer.def( 
                "nChains"
                , nChains_function_value );
        
        }
        { //::SireMol::PartialMolecule::nCutGroups
        
            typedef int ( ::SireMol::PartialMolecule::*nCutGroups_function_type )(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::PartialMolecule::nCutGroups );
            
            PartialMolecule_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value );
        
        }
        { //::SireMol::PartialMolecule::nResidues
        
            typedef int ( ::SireMol::PartialMolecule::*nResidues_function_type )(  ) const;
            nResidues_function_type nResidues_function_value( &::SireMol::PartialMolecule::nResidues );
            
            PartialMolecule_exposer.def( 
                "nResidues"
                , nResidues_function_value );
        
        }
        { //::SireMol::PartialMolecule::nSegments
        
            typedef int ( ::SireMol::PartialMolecule::*nSegments_function_type )(  ) const;
            nSegments_function_type nSegments_function_value( &::SireMol::PartialMolecule::nSegments );
            
            PartialMolecule_exposer.def( 
                "nSegments"
                , nSegments_function_value );
        
        }
        { //::SireMol::PartialMolecule::name
        
            typedef ::SireMol::MolName const & ( ::SireMol::PartialMolecule::*name_function_type )(  ) const;
            name_function_type name_function_value( &::SireMol::PartialMolecule::name );
            
            PartialMolecule_exposer.def( 
                "name"
                , name_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMol::PartialMolecule::number
        
            typedef ::SireMol::MolNum ( ::SireMol::PartialMolecule::*number_function_type )(  ) const;
            number_function_type number_function_value( &::SireMol::PartialMolecule::number );
            
            PartialMolecule_exposer.def( 
                "number"
                , number_function_value );
        
        }
        PartialMolecule_exposer.def( bp::self != bp::self );
        { //::SireMol::PartialMolecule::operator=
        
            typedef ::SireMol::PartialMolecule & ( ::SireMol::PartialMolecule::*assign_function_type )( ::SireMol::MoleculeView const & ) ;
            assign_function_type assign_function_value( &::SireMol::PartialMolecule::operator= );
            
            PartialMolecule_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::PartialMolecule::operator=
        
            typedef ::SireMol::PartialMolecule & ( ::SireMol::PartialMolecule::*assign_function_type )( ::SireMol::PartialMolecule const & ) ;
            assign_function_type assign_function_value( &::SireMol::PartialMolecule::operator= );
            
            PartialMolecule_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        PartialMolecule_exposer.def( bp::self == bp::self );
        { //::SireMol::PartialMolecule::property
        
            typedef ::SireBase::Property const & ( ::SireMol::PartialMolecule::*property_function_type )( ::SireBase::PropertyName const & ) const;
            property_function_type property_function_value( &::SireMol::PartialMolecule::property );
            
            PartialMolecule_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMol::PartialMolecule::propertyKeys
        
            typedef ::QStringList ( ::SireMol::PartialMolecule::*propertyKeys_function_type )(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::PartialMolecule::propertyKeys );
            
            PartialMolecule_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value );
        
        }
        { //::SireMol::PartialMolecule::selectedAll
        
            typedef bool ( ::SireMol::PartialMolecule::*selectedAll_function_type )(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::PartialMolecule::selectedAll );
            
            PartialMolecule_exposer.def( 
                "selectedAll"
                , selectedAll_function_value );
        
        }
        { //::SireMol::PartialMolecule::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMol::PartialMolecule::*selection_function_type )(  ) const;
            selection_function_type selection_function_value( &::SireMol::PartialMolecule::selection );
            
            PartialMolecule_exposer.def( 
                "selection"
                , selection_function_value );
        
        }
        { //::SireMol::PartialMolecule::toString
        
            typedef ::QString ( ::SireMol::PartialMolecule::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::PartialMolecule::toString );
            
            PartialMolecule_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::PartialMolecule::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::PartialMolecule::typeName );
            
            PartialMolecule_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::PartialMolecule::version
        
            typedef ::quint64 ( ::SireMol::PartialMolecule::*version_function_type )(  ) const;
            version_function_type version_function_value( &::SireMol::PartialMolecule::version );
            
            PartialMolecule_exposer.def( 
                "version"
                , version_function_value );
        
        }
        { //::SireMol::PartialMolecule::version
        
            typedef ::quint64 ( ::SireMol::PartialMolecule::*version_function_type )( ::SireBase::PropertyName const & ) const;
            version_function_type version_function_value( &::SireMol::PartialMolecule::version );
            
            PartialMolecule_exposer.def( 
                "version"
                , version_function_value
                , ( bp::arg("key") ) );
        
        }
        PartialMolecule_exposer.staticmethod( "typeName" );
        PartialMolecule_exposer.def( "__copy__", &__copy__);
        PartialMolecule_exposer.def( "__deepcopy__", &__copy__);
        PartialMolecule_exposer.def( "clone", &__copy__);
        PartialMolecule_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::PartialMolecule >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PartialMolecule_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::PartialMolecule >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PartialMolecule_exposer.def( "__str__", &__str__< ::SireMol::PartialMolecule > );
        PartialMolecule_exposer.def( "__repr__", &__str__< ::SireMol::PartialMolecule > );
    }

}
